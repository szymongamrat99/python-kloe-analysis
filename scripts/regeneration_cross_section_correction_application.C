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

gROOT->SetBatch(kTRUE);

void regeneration_cross_section_after_DC_displacement(TString root_file_date)
{
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

    double t_min = -300, t_max = 300, t_res = 1.5, r_min = 0., r_max = 50.;
    int t_bins = floor((t_max - t_min) / t_res);
    int r_bins = 200;

    h_dt[channel_name] = new TH1F(("h_dt_" + channel_name).c_str(), ("Time difference " + channel_name + ";t_{ch} - t_{ne} [#tau_{S}];Entries").c_str(), t_bins, t_min, t_max, r_bins, r_min, r_max);
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

  std::unordered_map<std::string, TH2 *> &h_corr_matrices
  // Open correction file
  initialize_matrices_from_file("results/regeneration_cross_section_corr_factors.root", h_corr_matrices);


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

    double signal_weight = 1.0;

    if (current_channel == "Signal")
          signal_weight = interf.interf_function(t_ch - t_ne) / interf.double_exponential(t_ch - t_ne);

    kch_position.SetXYZ(Kchrec[6], Kchrec[7], Kchrec[8]);
    ev.R_charged = kch_position.Mag();
    ev.rho_charged = kch_position.Perp();

    kne_position.SetXYZ(Knerec[6], Knerec[7], Knerec[8]);
    ev.R_neutral = kne_position.Mag();
    ev.rho_neutral = kne_position.Perp();

    bool CylBPCharged = (ev.rho_charged > 3.60227 && ev.rho_charged < 5.4366),
         CylDCCharged = (ev.rho_charged > 23.4879 && ev.rho_charged < 26.8633),
         SphBPCharged = (ev.R_charged > 7.57729 && ev.R_charged < 12.3962),
         CylBPNeutral = (ev.rho_neutral > 3.34375 && ev.rho_neutral < 5.41215),
         CylDCNeutral = (ev.rho_neutral > 23.3831 && ev.rho_neutral < 26.8272),
         SphBPNeutral = (ev.R_neutral > 7.91582 && ev.R_neutral < 12.334);

    
    if (current_channel == "Regeneration")
    {
      if (CylBPCharged || CylDCCharged)
        
    }
    else
    {
      
    }

    // Wypełnianie histogramów współrzędnych (dla analizy jakościowej)
    for (int i = 0; i < 3; ++i)
    {
      h_neu_coordinates[current_channel][i]->Fill(Knerec[6 + i] - (is_mc ? displacement_vector_mc[i] : displacement_vector_data[i]));
      h_ch_coordinates[current_channel][i]->Fill(Kchrec[6 + i] - (is_mc ? displacement_vector_mc[i] : displacement_vector_data[i]));
    }

    *mcflag == 1 ? mc_entries++ : data_entries++;
  }

  // Rescaling of Signal to original entries
  h_rho_charged_no_corr["Signal"]->Scale(h_rho_charged_no_corr["Signal"]->GetEntries() / h_rho_charged_no_corr["Signal"]->Integral(0, h_rho_charged_no_corr["Signal"]->GetNbinsX() + 1));
  h_rho_neutral_no_corr["Signal"]->Scale(h_rho_neutral_no_corr["Signal"]->GetEntries() / h_rho_neutral_no_corr["Signal"]->Integral(0, h_rho_neutral_no_corr["Signal"]->GetNbinsX() + 1));

  h_R_charged_no_corr["Signal"]->Scale(h_R_charged_no_corr["Signal"]->GetEntries() / h_R_charged_no_corr["Signal"]->Integral(0, h_R_charged_no_corr["Signal"]->GetNbinsX() + 1));
  h_R_neutral_no_corr["Signal"]->Scale(h_R_neutral_no_corr["Signal"]->GetEntries() / h_R_neutral_no_corr["Signal"]->Integral(0, h_R_neutral_no_corr["Signal"]->GetNbinsX() + 1));

  h_corr_rho_dt_charged["Signal"]->Scale(h_corr_rho_dt_charged["Signal"]->GetEntries() / h_corr_rho_dt_charged["Signal"]->Integral(0, h_corr_rho_dt_charged["Signal"]->GetNbinsX() + 1, 0, h_corr_rho_dt_charged["Signal"]->GetNbinsY() + 1));
  h_corr_rho_dt_neutral["Signal"]->Scale(h_corr_rho_dt_neutral["Signal"]->GetEntries() / h_corr_rho_dt_neutral["Signal"]->Integral(0, h_corr_rho_dt_neutral["Signal"]->GetNbinsX() + 1, 0, h_corr_rho_dt_neutral["Signal"]->GetNbinsY() + 1));

  h_corr_R_dt_charged["Signal"]->Scale(h_corr_R_dt_charged["Signal"]->GetEntries() / h_corr_R_dt_charged["Signal"]->Integral(0, h_corr_R_dt_charged["Signal"]->GetNbinsX() + 1, 0, h_corr_R_dt_charged["Signal"]->GetNbinsY() + 1));
  h_corr_R_dt_neutral["Signal"]->Scale(h_corr_R_dt_neutral["Signal"]->GetEntries() / h_corr_R_dt_neutral["Signal"]->Integral(0, h_corr_R_dt_neutral["Signal"]->GetNbinsX() + 1, 0, h_corr_R_dt_neutral["Signal"]->GetNbinsY() + 1));

  // Adding Signal to MC sum
  h_rho_charged_no_corr["MC sum"]->Add(h_rho_charged_no_corr["Signal"]);
  h_rho_neutral_no_corr["MC sum"]->Add(h_rho_neutral_no_corr["Signal"]);

  h_R_charged_no_corr["MC sum"]->Add(h_R_charged_no_corr["Signal"]);
  h_R_neutral_no_corr["MC sum"]->Add(h_R_neutral_no_corr["Signal"]);

  h_corr_rho_dt_charged["MC sum"]->Add(h_corr_rho_dt_charged["Signal"]);
  h_corr_rho_dt_neutral["MC sum"]->Add(h_corr_rho_dt_neutral["Signal"]);

  h_corr_R_dt_charged["MC sum"]->Add(h_corr_R_dt_charged["Signal"]);
  h_corr_R_dt_neutral["MC sum"]->Add(h_corr_R_dt_neutral["Signal"]);


  std::cout << "DC Wall width, no correction, charged" << std::endl;
  KloeTools::quantile_limits(h_rho_charged_no_corr["Regeneration"], 23, 27, u0, v0);
  std::cout << "Cylindrical BP Wall width, no correction, charged" << std::endl;
  KloeTools::quantile_limits(h_rho_charged_no_corr["Regeneration"], 3.0, 5.5, u0, v0);
  std::cout << "Spherical BP Wall width, no correction, charged" << std::endl;
  KloeTools::quantile_limits(h_R_charged_no_corr["Regeneration"], 7.0, 13.0, u0, v0);

  std::cout << "DC Wall width, no correction, neutral" << std::endl;
  KloeTools::quantile_limits(h_rho_neutral_no_corr["Regeneration"], 23, 27, u0, v0);
  std::cout << "Cylindrical BP Wall width, no correction, neutral" << std::endl;
  KloeTools::quantile_limits(h_rho_neutral_no_corr["Regeneration"], 3.0, 5.5, u0, v0);
  std::cout << "Spherical BP Wall width, no correction, neutral" << std::endl;
  KloeTools::quantile_limits(h_R_neutral_no_corr["Regeneration"], 7.0, 13.0, u0, v0);

  std::cout << "DC Wall width, corrected, charged" << std::endl;
  KloeTools::quantile_limits(h_rho_charged["Regeneration"], 23, 27, u0, v0);
  std::cout << "Cylindrical BP Wall width, corrected, charged" << std::endl;
  KloeTools::quantile_limits(h_rho_charged["Regeneration"], 3.0, 5.5, u0, v0);
  std::cout << "Spherical BP Wall width, corrected, charged" << std::endl;
  KloeTools::quantile_limits(h_R_charged["Regeneration"], 7.0, 13.0, u0, v0);

  std::cout << "DC Wall width, neutral" << std::endl;
  KloeTools::quantile_limits(h_rho_neutral["Regeneration"], 23, 27, u0, v0);
  std::cout << "Cylindrical BP Wall width, neutral" << std::endl;
  KloeTools::quantile_limits(h_rho_neutral["Regeneration"], 3.0, 5.5, u0, v0);
  std::cout << "Spherical BP Wall width, neutral" << std::endl;
  KloeTools::quantile_limits(h_R_neutral["Regeneration"], 7.0, 13.0, u0, v0);
  

  // Zapis wyników do pliku ROOT
  TFile output_file("regeneration_cross_section_DC.root", "RECREATE");
  for (const auto &chann : KLOE::channName)
  {
    std::string channel_name = (std::string)chann.second;
    h_rho_charged[channel_name]->Write();
    h_rho_neutral[channel_name]->Write();
    h_R_charged[channel_name]->Write();
    h_R_neutral[channel_name]->Write();
    h_rho_charged_vs_R[channel_name]->Write();
    h_rho_neutral_vs_R[channel_name]->Write();
    for (int i = 0; i < 3; ++i)
    {
      h_neu_coordinates[channel_name][i]->Write();
      h_ch_coordinates[channel_name][i]->Write();
    }

    h_xy_charged_DC[channel_name]->Write();
    h_xy_charged_BP[channel_name]->Write();

    // Automatyczny zapis nowych profili wraz ze zintegrowaną linią dopasowania sinusoidalnego!
    h_prof_polar_DC[channel_name]->Write();
    h_prof_polar_BP[channel_name]->Write();

    h_polar_DC[channel_name]->Write();
    h_polar_BP[channel_name]->Write();

    // --- TUTAJ ZACZYNA SIĘ TWOJA ISTNIEJĄCA PĘTLA RYSOWANIA ---
    TCanvas *c = new TCanvas(("c_" + channel_name).c_str(), ("Canvas " + channel_name).c_str(), 800, 800);
    c->cd();

    // 3. POZOSTAŁE HISTOGRAMY (Rysujesz tak jak wcześniej, dodaj tylko c->Clear() dla porządku)
    h_rho_charged[channel_name]->SetLineColor(kBlue);
    h_rho_charged_no_corr[channel_name]->SetLineColor(kRed);
    h_rho_charged[channel_name]->Draw();
    h_rho_charged_no_corr[channel_name]->Draw("SAME");
    c->SaveAs(("img/rho_charged_" + channel_name + ".png").c_str());

    c->Clear();
    h_rho_neutral[channel_name]->SetLineColor(kBlue);
    h_rho_neutral_no_corr[channel_name]->SetLineColor(kRed);
    h_rho_neutral[channel_name]->Draw();
    h_rho_neutral_no_corr[channel_name]->Draw("SAME");
    c->SaveAs(("img/rho_neutral_" + channel_name + ".png").c_str());

    c->Clear();
    h_R_charged[channel_name]->SetLineColor(kBlue);
    h_R_charged_no_corr[channel_name]->SetLineColor(kRed);
    h_R_charged[channel_name]->Draw();
    h_R_charged_no_corr[channel_name]->Draw("SAME");
    c->SaveAs(("img/R_charged_" + channel_name + ".png").c_str());

    c->Clear();
    h_R_neutral[channel_name]->SetLineColor(kBlue);
    h_R_neutral_no_corr[channel_name]->SetLineColor(kRed);
    h_R_neutral[channel_name]->Draw();
    h_R_neutral_no_corr[channel_name]->Draw("SAME");
    c->SaveAs(("img/R_neutral_" + channel_name + ".png").c_str());

    c->Clear();
    h_corr_rho_dt_charged[channel_name]->Draw("COLZ");
    c->SaveAs(("img/corr_rho_dt_charged_" + channel_name + ".png").c_str());

    c->Clear();
    h_corr_rho_dt_neutral[channel_name]->Draw("COLZ");
    c->SaveAs(("img/corr_rho_dt_neutral_" + channel_name + ".png").c_str());

    c->Clear();
    h_corr_R_dt_charged[channel_name]->Draw("COLZ");
    c->SaveAs(("img/corr_R_dt_charged_" + channel_name + ".png").c_str());

    c->Clear();
    h_corr_R_dt_neutral[channel_name]->Draw("COLZ");
    c->SaveAs(("img/corr_R_dt_neutral_" + channel_name + ".png").c_str());

    delete c;
    // --- KONIEC BLOKU ---
  }

  output_file.Close();

  double scaling_factor_lumi = data_entries / static_cast<double>(mc_entries);

  h_rho_charged_no_corr["Data"]->Sumw2(); // Upewniamy się, że błędy są poprawnie przeliczone dla danych bez korekcji
  h_rho_charged_no_corr["MC sum"]->Sumw2(); // Upewniamy się, że błędy są poprawnie przeliczone dla sumy
  h_rho_neutral_no_corr["Data"]->Sumw2();
  h_rho_neutral_no_corr["MC sum"]->Sumw2();
  h_R_charged_no_corr["Data"]->Sumw2();
  h_R_charged_no_corr["MC sum"]->Sumw2();
  h_R_neutral_no_corr["Data"]->Sumw2();
  h_R_neutral_no_corr["MC sum"]->Sumw2();

  h_rho_charged_no_corr["MC sum"]->Scale(scaling_factor_lumi); // Skaluje MC do danych

  KloeTools::plot_ratio_and_save(h_rho_charged_no_corr["Data"], h_rho_charged_no_corr["MC sum"], "Ratio: rho_charged", "img/ratio_rho_charged.png");
  KloeTools::plot_ratio_and_save(h_rho_neutral_no_corr["Data"], h_rho_neutral_no_corr["MC sum"], "Ratio: rho_neutral", "img/ratio_rho_neutral.png");

  KloeTools::plot_ratio_and_save(h_R_charged_no_corr["Data"], h_R_charged_no_corr["MC sum"], "Ratio: R_charged", "img/ratio_R_charged.png");
  KloeTools::plot_ratio_and_save(h_R_neutral_no_corr["Data"], h_R_neutral_no_corr["MC sum"], "Ratio: R_neutral", "img/ratio_R_neutral.png");

  // 2D ratios
  h_corr_rho_dt_charged["MC sum"]->Scale(scaling_factor_lumi); // Skaluje MC do danych
  h_corr_rho_dt_neutral["MC sum"]->Scale(scaling_factor_lumi);
  h_corr_R_dt_charged["MC sum"]->Scale(scaling_factor_lumi);
  h_corr_R_dt_neutral["MC sum"]->Scale(scaling_factor_lumi);

  h_corr_rho_dt_charged["Data"]->Sumw2(); // Upewniamy się, że błędy są poprawnie przeliczone dla danych
  h_corr_rho_dt_charged["MC sum"]->Sumw2(); // Upewniamy się, że błędy są poprawnie przeliczone dla sumy
  h_corr_rho_dt_neutral["Data"]->Sumw2();
  h_corr_rho_dt_neutral["MC sum"]->Sumw2();
  h_corr_R_dt_charged["Data"]->Sumw2();
  h_corr_R_dt_charged["MC sum"]->Sumw2();
  h_corr_R_dt_neutral["Data"]->Sumw2();
  h_corr_R_dt_neutral["MC sum"]->Sumw2();

  h_corr_rho_dt_charged["MC sum"]->Add(h_corr_rho_dt_neutral["MC sum"]); // Dodaj neutralne do naładowania sumy
  h_corr_R_dt_charged["MC sum"]->Add(h_corr_R_dt_neutral["MC sum"]);
  h_corr_rho_dt_charged["Data"]->Add(h_corr_rho_dt_neutral["Data"]); // Dodaj neutralne do naładowania sumy
  h_corr_R_dt_charged["Data"]->Add(h_corr_R_dt_neutral["Data"]);

  TH2 *h_ratio_rho_dt = KloeTools::plot_ratio_2D_and_save(h_corr_rho_dt_charged["Data"], h_corr_rho_dt_charged["MC sum"], "Ratio: rho_charged vs dt", "img/ratio_rho_charged_dt.png");
  TH2 *h_ratio_R_dt = KloeTools::plot_ratio_2D_and_save(h_corr_R_dt_charged["Data"], h_corr_R_dt_charged["MC sum"], "Ratio: R_charged vs dt", "img/ratio_R_charged_dt.png");

  TFile output_file_ratio("results/regeneration_cross_section_corr_factors.root", "RECREATE");
  h_ratio_rho_dt->Write("h_ratio_rho_dt");
  h_ratio_R_dt->Write("h_ratio_R_dt");
  output_file_ratio.Close();
}
