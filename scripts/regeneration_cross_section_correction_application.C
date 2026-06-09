#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TChain.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TMath.h>

#include <string>
#include <array>
#include <unordered_map>

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TVector3.h>
#include <TF1.h>
#include <TLine.h>

#include <TFitResult.h>
#include <TFitResultPtr.h>

#include <TLegend.h>

#include "const.h"
#include "ErrorLogs.h"
#include "interference.h"
#include "inc/KLOETools.h"

void regeneration_cross_section_correction_application(TString root_file_date)
{
  gROOT->SetBatch(kTRUE);
  gErrorIgnoreLevel = kSysError;
  KLOE::setGlobalStyle();
  gStyle->SetOptStat(0);

  ErrorHandling::ErrorLogs logger("logs/");

  Utils::InitializeVariables(logger);
  KLOE::interference interf;

  TChain *chain = new TChain("h1");
  KloeTools::LoadSignalChains(chain, root_file_date);

  TTreeReader reader(chain);
  TTreeReaderValue<Double_t> time_ch_CM_signal_fit(reader, "KaonChTimeCMSignalFit");
  TTreeReaderValue<Double_t> time_ne_CM_signal_fit(reader, "KaonNeTimeCMSignalFit");
  TTreeReaderValue<Double_t> time_ch_CM_MC(reader, "KaonChTimeCMMC");
  TTreeReaderValue<Double_t> time_ne_CM_MC(reader, "KaonNeTimeCMMC");

  TTreeReaderValue<Double_t> Bx(reader, "Bx");
  TTreeReaderValue<Double_t> By(reader, "By");
  TTreeReaderValue<Double_t> Bz(reader, "Bz");

  TTreeReaderValue<Double_t> minv4gam(reader, "minv4gam");

  TTreeReaderArray<Double_t> Knerec(reader, "Knerec");
  TTreeReaderArray<Double_t> KnerecSix(reader, "KnerecSix");
  TTreeReaderArray<Double_t> Kchrec(reader, "Kchrec");

  TTreeReaderArray<Double_t> KnerecFit(reader, "KnerecFit");
  TTreeReaderArray<Double_t> KchrecFit(reader, "KchrecFit");

  TTreeReaderArray<Double_t> trk1Fit(reader, "trk1Fit");
  TTreeReaderArray<Double_t> trk2Fit(reader, "trk2Fit");

  TTreeReaderArray<Double_t> KchrecClosest(reader, "KchrecClosest");

  TTreeReaderArray<Double_t> ip(reader, "ip");

  TTreeReaderValue<Int_t> mctruth(reader, "mctruth");
  TTreeReaderValue<Int_t> mcflag(reader, "mcflag");

  // Definition of histograms
  std::unordered_map<std::string, TH1 *> h_dt;

  for (const auto &chann : KLOE::channName)
  {
    std::string channel_name = (std::string)chann.second;

    double t_min = -90, t_max = 90, t_res = 1.5, r_min = 0., r_max = 50.;
    int t_bins = floor((t_max - t_min) / t_res);
    int r_bins = 200;

    h_dt[channel_name] = new TH1F(("h_dt_" + channel_name).c_str(), ("Time difference " + channel_name + ";t_{ch} - t_{ne} [#tau_{S}];Entries").c_str(), t_bins, t_min, t_max);
  }
  //

  // Cuts
  double u0 = -92.9524, v0 = 7.1644, su = 33.7633, sv = 42.0163, rho = 0.556;
  double n_sigma_cut = 2.5;

  // Definicja lambdy
  auto ellipse_cut = [=](double u, double v)
  {
    double du = u - u0;
    double dv = v - v0;
    double inv_rho2 = 1.0 / (1.0 - rho * rho);

    double dist2 = inv_rho2 * ((du * du) / (su * su) +
                               (dv * dv) / (sv * sv) -
                               2.0 * rho * (du * dv) / (su * sv));

    return std::sqrt(std::max(0.0, dist2)) > n_sigma_cut;
  };
  // ---

  std::unordered_map<std::string, Int_t> entries_count; // Map to count entries for each channel

  Long64_t data_entries = 0, mc_entries = 0;

  std::unordered_map<std::string, TH2 *> h_corr_matrices;
  // Open correction file
  KloeTools::initialize_matrices_from_file("results/regeneration_cross_section_corr_factors.root", h_corr_matrices);

  TVector3 kch_position, kne_position;

  while (reader.Next())
  {
    // 1. Wstępna identyfikacja nazwy kanału (Czytelna i bezpieczna!)
    std::string current_channel = "Data";
    if (*mcflag == 1)
    {
      int truth = *mctruth;
      if (truth == -1 || truth == 0)
        continue;
      current_channel = KLOE::channName.at(truth).Data();
    }

    double t_ne = *time_ne_CM_signal_fit;
    double t_ch = *time_ch_CM_signal_fit;

    double t_ne_mc = *time_ne_CM_MC;
    double t_ch_mc = *time_ch_CM_MC;

    // 3. Wywołanie wydzielonej funkcji cięć (Kod jest czysty i czytelny!)
    bool global_cut = KloeTools::PassesGlobalCuts(KchrecFit, KnerecFit, Knerec, KchrecClosest, trk1Fit, trk2Fit,
                                                  *Bx, *By, *Bz, *minv4gam, KnerecSix[5], ellipse_cut);
    if (!global_cut)
      continue;

    EventKinematics ev;

    // ---

    if (*mcflag == 1 && (*mctruth == -1 || *mctruth == 0))
      continue;

    Bool_t is_mc = (*mcflag == 1);

    double signal_weight = 1.0, regeneration_weight = 1.0;

    if (current_channel == "Signal")
      signal_weight = interf.interf_function(t_ch_mc - t_ne_mc) / interf.double_exponential(t_ch_mc - t_ne_mc);

    kch_position.SetXYZ(Kchrec[6], Kchrec[7], Kchrec[8]);
    ev.R_charged = kch_position.Mag();
    ev.rho_charged = kch_position.Perp();

    kne_position.SetXYZ(Knerec[6], Knerec[7], Knerec[8]);
    ev.R_neutral = kne_position.Mag();
    ev.rho_neutral = kne_position.Perp();

    bool CylBPCharged = (ev.rho_charged > 3.60227 && ev.rho_charged < 5.4366),
         CylDCCharged = (ev.rho_charged > 23.4879 && ev.rho_charged < 26.8633),
         SphBPCharged = (ev.R_charged > 5.57729 && ev.R_charged < 15.3962),
         CylBPNeutral = (ev.rho_neutral > 3.34375 && ev.rho_neutral < 5.41215),
         CylDCNeutral = (ev.rho_neutral > 23.3831 && ev.rho_neutral < 26.8272),
         SphBPNeutral = (ev.R_neutral > 5.91582 && ev.R_neutral < 15.334);

    int n_regions = (int)CylBPCharged + (int)CylDCCharged + (int)SphBPCharged +
                (int)CylBPNeutral + (int)CylDCNeutral + (int)SphBPNeutral;

    if (current_channel == "Regeneration")
    {
      std::cout << n_regions << " " << CylBPCharged << " " << CylDCCharged << " " << SphBPCharged << " " << CylBPNeutral << " " << CylDCNeutral << " " << SphBPNeutral << std::endl;

      if (n_regions >= 1)
      {

        if (CylBPCharged || CylDCCharged)
          regeneration_weight = KloeTools::get_weight_2D(h_corr_matrices["cylindrical"], ev.rho_charged, t_ch - t_ne);

        if (SphBPCharged)
          regeneration_weight = KloeTools::get_weight_2D(h_corr_matrices["spherical"], ev.R_charged, t_ch - t_ne);

        if (CylBPNeutral || CylDCNeutral)
          regeneration_weight = KloeTools::get_weight_2D(h_corr_matrices["cylindrical"], ev.rho_neutral, t_ch - t_ne);

        if (SphBPNeutral)
          regeneration_weight = KloeTools::get_weight_2D(h_corr_matrices["spherical"], ev.R_neutral, t_ch - t_ne);
      }
    }

    if (current_channel == "Signal")
      h_dt[current_channel]->Fill(t_ch - t_ne, signal_weight);
    else if (current_channel == "Regeneration")
      h_dt[current_channel]->Fill(t_ch - t_ne, regeneration_weight);
    else if (current_channel == "Data")
      h_dt[current_channel]->Fill(t_ch - t_ne);
    else if (current_channel != "Semileptonic" && current_channel != "3pi0")
    {
      h_dt[current_channel]->Fill(t_ch - t_ne);
      h_dt["MC sum"]->Fill(t_ch - t_ne);
    }

    *mcflag == 1 ? mc_entries++ : data_entries++;
  }

  // Rescaling of Signal to original entries
  h_dt["Signal"]->Scale(h_dt["Signal"]->GetEntries() / h_dt["Signal"]->Integral(0, h_dt["Signal"]->GetNbinsX() + 1));

  // Adding Signal to MC sum
  h_dt["MC sum"]->Add(h_dt["Signal"]);
  h_dt["MC sum"]->Add(h_dt["Regeneration"]);

  double scaling_factor_lumi = 1.09;
  h_dt["MC sum"]->Scale(scaling_factor_lumi);

  // Zapis wyników do pliku ROOT
  TFile output_file("regeneration_cross_section_DC.root", "RECREATE");
  for (const auto &chann : KLOE::channName)
  {
    std::string channel_name = (std::string)chann.second;
    h_dt[channel_name]->Write();
  }

  h_corr_matrices["spherical"]->Write("h_ratio_R_dt");
  h_corr_matrices["cylindrical"]->Write("h_ratio_rho_dt");

  output_file.Close();

  // --- TUTAJ ZACZYNA SIĘ TWOJA ISTNIEJĄCA PĘTLA RYSOWANIA ---
  TCanvas *c = new TCanvas("c_dt", "c_dt", 800, 800);
  c->cd();

  h_dt["Data"]->SetLineColor(KLOE::channColor.at("Data"));
  h_dt["Data"]->GetYaxis()->SetRangeUser(0.0, h_dt["Data"]->GetMaximum() * 1.5);
  h_dt["Data"]->Draw("PE1");

  for (const auto &chann : KLOE::channName)
  {
    std::string channel_name = (std::string)chann.second;

    if (channel_name != "Data")
    {
      h_dt[channel_name]->SetLineColor(KLOE::channColor.at(channel_name));
      h_dt[channel_name]->Draw("SAME HIST");
    }
  }

  c->SaveAs("img/dt.png");
  c->Clear();

  delete c;
}
