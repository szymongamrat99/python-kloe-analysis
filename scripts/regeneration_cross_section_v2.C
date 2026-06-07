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
#include "inc/KLOETools.h"

gROOT->SetBatch(kTRUE);

void regeneration_cross_section_v2(TString root_file_date)
{
  KLOE::setGlobalStyle();
  gStyle->SetOptStat(0);

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
  std::unordered_map<std::string, TH1 *> h_rho_charged;
  std::unordered_map<std::string, TH1 *> h_rho_neutral;

  std::unordered_map<std::string, TH1 *> h_R_charged;
  std::unordered_map<std::string, TH1 *> h_R_neutral;

  std::unordered_map<std::string, TH2 *> h_rho_charged_vs_R;
  std::unordered_map<std::string, TH2 *> h_rho_neutral_vs_R;

  std::unordered_map<std::string, TH2 *> h_xy_charged_DC;
  std::unordered_map<std::string, TH2 *> h_xy_charged_BP;

  std::unordered_map<std::string, TProfile *> h_prof_polar_DC;
  std::unordered_map<std::string, TProfile *> h_prof_polar_BP;

  std::unordered_map<std::string, TH2 *> h_polar_DC;
  std::unordered_map<std::string, TH2 *> h_polar_BP;

  std::unordered_map<std::string, std::array<TH1 *, 3>> h_neu_coordinates;
  std::unordered_map<std::string, std::array<TH1 *, 3>> h_ch_coordinates;

  for (const auto &chann : KLOE::channName)
  {
    std::string channel_name = (std::string)chann.second;

    h_rho_charged[channel_name] = new TH1F(("h_rho_charged_" + channel_name).c_str(), ("Rho Charged " + channel_name + ";#rho_{charged};Entries").c_str(), 251, -50, 50);
    h_rho_neutral[channel_name] = new TH1F(("h_rho_neutral_" + channel_name).c_str(), ("Rho Neutral " + channel_name + ";#rho_{neutral};Entries").c_str(), 251, -50, 50);

    h_R_charged[channel_name] = new TH1F(("h_R_charged_" + channel_name).c_str(), ("R Charged " + channel_name + ";R_{charged};Entries").c_str(), 251, 0, 50);
    h_R_neutral[channel_name] = new TH1F(("h_R_neutral_" + channel_name).c_str(), ("R Neutral " + channel_name + ";R_{neutral};Entries").c_str(), 251, 0, 50);

    h_rho_charged_vs_R[channel_name] = new TH2F(("h_rho_charged_vs_R_" + channel_name).c_str(), ("Rho Charged vs R " + channel_name + ";R_{charged};#rho_{charged}").c_str(), 251, 0, 50, 251, 0, 50);
    h_rho_neutral_vs_R[channel_name] = new TH2F(("h_rho_neutral_vs_R_" + channel_name).c_str(), ("Rho Neutral vs R " + channel_name + ";R_{neutral};#rho_{neutral}").c_str(), 251, 0, 50, 251, 0, 50);

    for (int i = 0; i < 3; ++i)
    {
      h_neu_coordinates[channel_name][i] = new TH1F(("h_neu_" + channel_name + "_" + std::to_string(i)).c_str(), ("Neutral Coordinate " + channel_name + " - " + std::to_string(i) + ";Coordinate;Entries").c_str(), 251, -50, 50);
      h_ch_coordinates[channel_name][i] = new TH1F(("h_ch_" + channel_name + "_" + std::to_string(i)).c_str(), ("Charged Coordinate " + channel_name + " - " + std::to_string(i) + ";Coordinate;Entries").c_str(), 251, -50, 50);
    }

    h_xy_charged_DC[channel_name] = new TH2F(("h_xy_charged_DC_" + channel_name).c_str(), ("XY Charged " + channel_name + ";X_{charged} [cm];Y_{charged} [cm]").c_str(), 251, -50, 50, 251, -50, 50);
    h_xy_charged_BP[channel_name] = new TH2F(("h_xy_charged_BP_" + channel_name).c_str(), ("XY Charged " + channel_name + " - BP;X_{charged} [cm];Y_{charged} [cm]").c_str(), 251, -50, 50, 251, -50, 50);

    h_prof_polar_DC[channel_name] = new TProfile(("h_prof_polar_DC_" + channel_name).c_str(), ("Polar Profile DC " + channel_name + ";#phi [rad];#rho [cm]").c_str(), 70, -TMath::Pi(), TMath::Pi());
    h_prof_polar_BP[channel_name] = new TProfile(("h_prof_polar_BP_" + channel_name).c_str(), ("Polar Profile BP " + channel_name + ";#phi [rad];#rho [cm]").c_str(), 70, -TMath::Pi(), TMath::Pi());

    h_prof_polar_DC[channel_name]->SetErrorOption("G");
    h_prof_polar_BP[channel_name]->SetErrorOption("G");

    h_polar_DC[channel_name] = new TH2F(("h_polar_DC_" + channel_name).c_str(), ("Polar DC " + channel_name + ";#phi [rad];#rho [cm]").c_str(), 30, -TMath::Pi(), TMath::Pi(), 100, 20, 30);
    h_polar_BP[channel_name] = new TH2F(("h_polar_BP_" + channel_name).c_str(), ("Polar BP " + channel_name + ";#phi [rad];#rho [cm]").c_str(), 30, -TMath::Pi(), TMath::Pi(), 100, 3, 6);
  }
  //

  TFitResultPtr result_mc, result_data;

  std::vector<double> displacement_MC = KloeTools::iterative_polar_fit(reader, KchrecFit, KnerecFit, Knerec, Kchrec, KchrecClosest, trk1Fit, trk2Fit, mcflag, mctruth, Bx, By, Bz, minv4gam, KnerecSix, h_prof_polar_DC["Regeneration"], "Regeneration", result_mc);
  std::vector<double> displacement_data = KloeTools::iterative_polar_fit(reader, KchrecFit, KnerecFit, Knerec, Kchrec, KchrecClosest, trk1Fit, trk2Fit, mcflag, mctruth, Bx, By, Bz, minv4gam, KnerecSix, h_prof_polar_DC["Data"], "Data", result_data);
  std::vector<double> displacement_mc_sum = KloeTools::iterative_polar_fit(reader, KchrecFit, KnerecFit, Knerec, Kchrec, KchrecClosest, trk1Fit, trk2Fit, mcflag, mctruth, Bx, By, Bz, minv4gam, KnerecSix, h_prof_polar_DC["MC sum"], "MC sum", result_data);


  h_prof_polar_DC["Regeneration"]->Reset();
  h_prof_polar_DC["Data"]->Reset();
  h_prof_polar_DC["MC sum"]->Reset();

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

  TVector3 kch_position(0, 0, 0),
      kne_position(0, 0, 0),
      displacement_vector_data(displacement_data[0], displacement_data[1], displacement_data[2]),
      displacement_vector_mc(displacement_MC[0], displacement_MC[1], displacement_MC[2]);

  std::unordered_map<std::string, Int_t> entries_count; // Map to count entries for each channel

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

    Bool_t is_mc = (*mcflag == 1);

    if (1)//t_ne < 7.0)
    {
      kch_position.SetXYZ(Kchrec[6], Kchrec[7], Kchrec[8]);

      if (is_mc)
        kch_position -= displacement_vector_mc;
      else
        kch_position -= displacement_vector_data;

      ev.R_charged = kch_position.Mag();
      ev.rho_charged = kch_position.Perp();

      // Obliczenie kąta polarnego w płaszczyźnie XY detektora
      double phi_charged = std::atan2(Kchrec[7], Kchrec[6]); // atan2(y, x) dla poprawnego zakresu kąta
      double r_geom_xy = std::sqrt(Kchrec[6] * Kchrec[6] + Kchrec[7] * Kchrec[7]);

      double rho = 0;//std::sqrt(pow(Kchrec[6], 2) + pow(Kchrec[7] + 0.8, 2));

      if (current_channel == "Data")
      {
        rho = std::sqrt(pow(Kchrec[6] - displacement_vector_data.X(), 2) + pow(Kchrec[7] - displacement_vector_data.Y(), 2));
      }
      else
      {
        rho = std::sqrt(pow(Kchrec[6] - displacement_vector_mc.X(), 2) + pow(Kchrec[7] - displacement_vector_mc.Y(), 2));
      }

      if (rho < 27 && rho > 23.0)
      {
        h_xy_charged_DC[current_channel]->Fill(kch_position.X(), kch_position.Y());
        h_prof_polar_DC[current_channel]->Fill(phi_charged, r_geom_xy); // Napełnianie profilu dla DC
        h_polar_DC[current_channel]->Fill(phi_charged, r_geom_xy); // Napełnianie histogramu 2D dla DC (do porównania)

        if (current_channel != "Data")
        {
          h_xy_charged_DC["MC sum"]->Fill(kch_position.X(), kch_position.Y());
          h_prof_polar_DC["MC sum"]->Fill(phi_charged, r_geom_xy);
          h_polar_DC["MC sum"]->Fill(phi_charged, r_geom_xy);
        }
      }

      if (r_geom_xy < 5.5 && r_geom_xy > 3.5)
      {
        h_xy_charged_BP[current_channel]->Fill(kch_position.X(), kch_position.Y());
        h_prof_polar_BP[current_channel]->Fill(phi_charged, r_geom_xy); // Napełnianie profilu dla BP
        h_polar_BP[current_channel]->Fill(phi_charged, r_geom_xy); // Napełnianie histogramu 2D dla BP (do porównania)

        if (current_channel != "Data")
        {
          h_xy_charged_BP["MC sum"]->Fill(kch_position.X(), kch_position.Y());
          h_prof_polar_BP["MC sum"]->Fill(phi_charged, r_geom_xy);
          h_polar_BP["MC sum"]->Fill(phi_charged, r_geom_xy);
        }
      }

      h_rho_charged[current_channel]->Fill(ev.rho_charged);
      h_R_charged[current_channel]->Fill(ev.R_charged);
      h_rho_charged_vs_R[current_channel]->Fill(ev.R_charged, ev.rho_charged);
    }

    if (t_ch < 7.0)
    {
      kne_position.SetXYZ(Knerec[6], Knerec[7], Knerec[8]);

      if (is_mc)
        kne_position -= displacement_vector_mc;
      else
        kne_position -= displacement_vector_data;

      ev.R_neutral = kne_position.Mag();
      ev.rho_neutral = kne_position.Perp();

      // Wypełnianie po nazwie kanału!
      h_rho_neutral[current_channel]->Fill(ev.rho_neutral);
      h_R_neutral[current_channel]->Fill(ev.R_neutral);
      h_rho_neutral_vs_R[current_channel]->Fill(ev.R_neutral, ev.rho_neutral);
    }

    // Wypełnianie histogramów współrzędnych (dla analizy jakościowej)
    for (int i = 0; i < 3; ++i)
    {
      h_neu_coordinates[current_channel][i]->Fill(Knerec[6 + i] - (is_mc ? displacement_vector_mc[i] : displacement_vector_data[i]));
      h_ch_coordinates[current_channel][i]->Fill(Kchrec[6 + i] - (is_mc ? displacement_vector_mc[i] : displacement_vector_data[i]));
    }
  }

  // Find local maxima for cylindrical radius histograms
  std::string channel = "Regeneration";
  // 1. Stare wywołanie dopasowania 2D na histogramach (pozostawione dla porównania grafiki)
  TFitResultPtr fit_result_2D_DC_Regen = KloeTools::fit_method_circle(h_xy_charged_DC[channel], -50, 50);
  TFitResultPtr fit_result_2D_BP_Regen = KloeTools::fit_method_circle(h_xy_charged_BP[channel], -50, 50);

  double x0_DC_2D = fit_result_2D_DC_Regen->IsValid() ? fit_result_2D_DC_Regen->Parameter(0) : 0.0;
  double x0_BP_2D = fit_result_2D_BP_Regen->IsValid() ? fit_result_2D_BP_Regen->Parameter(0) : 0.0;
  double y0_DC_2D = fit_result_2D_DC_Regen->IsValid() ? fit_result_2D_DC_Regen->Parameter(1) : 0.0;
  double y0_BP_2D = fit_result_2D_BP_Regen->IsValid() ? fit_result_2D_BP_Regen->Parameter(1) : 0.0;
  double r_DC_2D = fit_result_2D_DC_Regen->IsValid() ? fit_result_2D_DC_Regen->Parameter(2) : 0.0;
  double r_BP_2D = fit_result_2D_BP_Regen->IsValid() ? fit_result_2D_BP_Regen->Parameter(2) : 0.0;

  double x0_DC_2D_err = fit_result_2D_DC_Regen->IsValid() ? fit_result_2D_DC_Regen->ParError(0) : 0.0;
  double x0_BP_2D_err = fit_result_2D_BP_Regen->IsValid() ? fit_result_2D_BP_Regen->ParError(0) : 0.0;
  double y0_DC_2D_err = fit_result_2D_DC_Regen->IsValid() ? fit_result_2D_DC_Regen->ParError(1) : 0.0;
  double y0_BP_2D_err = fit_result_2D_BP_Regen->IsValid() ? fit_result_2D_BP_Regen->ParError(1) : 0.0;
  double r_DC_2D_err = fit_result_2D_DC_Regen->IsValid() ? fit_result_2D_DC_Regen->ParError(2) : 0.0;
  double r_BP_2D_err = fit_result_2D_BP_Regen->IsValid() ? fit_result_2D_BP_Regen->ParError(2) : 0.0;

  // 2. NOWE, STABILNE WYWOŁANIE: Dopasowanie sinusoidy na profilach polarnych
  // Przekazujemy promień spodziewany jako punkt startowy dla MINUITa (26.0 cm dla DC oraz 4.5 cm dla BP)
  TFitResultPtr fit_result_DC_Regen = KloeTools::fit_method_profile(h_prof_polar_DC[channel], 26.0);
  TFitResultPtr fit_result_BP_Regen = KloeTools::fit_method_profile(h_prof_polar_BP[channel], 4.5);

  double r_DC_prof = fit_result_DC_Regen->IsValid() ? fit_result_DC_Regen->Parameter(0) : 0.0;
  double r_BP_prof = fit_result_BP_Regen->IsValid() ? fit_result_BP_Regen->Parameter(0) : 0.0;
  double r_DC_prof_err = fit_result_DC_Regen->IsValid() ? fit_result_DC_Regen->ParError(0) : 0.0;
  double r_BP_prof_err = fit_result_BP_Regen->IsValid() ? fit_result_BP_Regen->ParError(0) : 0.0;

  double x0_prof_DC = fit_result_DC_Regen->IsValid() ? fit_result_DC_Regen->Parameter(1) : 0.0;
  double x0_prof_BP = fit_result_BP_Regen->IsValid() ? fit_result_BP_Regen->Parameter(1) : 0.0;
  double x0_prof_DC_err = fit_result_DC_Regen->IsValid() ? fit_result_DC_Regen->ParError(1) : 0.0;
  double x0_prof_BP_err = fit_result_BP_Regen->IsValid() ? fit_result_BP_Regen->ParError(1) : 0.0;

  double y0_prof_DC = fit_result_DC_Regen->IsValid() ? fit_result_DC_Regen->Parameter(2) : 0.0;
  double y0_prof_BP = fit_result_BP_Regen->IsValid() ? fit_result_BP_Regen->Parameter(2) : 0.0;
  double y0_prof_DC_err = fit_result_DC_Regen->IsValid() ? fit_result_DC_Regen->ParError(2) : 0.0;
  double y0_prof_BP_err = fit_result_BP_Regen->IsValid() ? fit_result_BP_Regen->ParError(2) : 0.0;

  std::cout << "\n>>> WYNIKI PORÓWNAWCZE ALIGNMENTU (Y0 Offset), Regeneracja:" << std::endl;
  std::cout << ">>> [DC] Metoda tradycyjna 2D: (" << y0_DC_2D << " +- " << y0_DC_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << y0_prof_DC << " +- " << y0_prof_DC_err << ") cm" << std::endl;
  std::cout << ">>> [BP] Metoda tradycyjna 2D: (" << y0_BP_2D << " +- " << y0_BP_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << y0_prof_BP << " +- " << y0_prof_BP_err << ") cm" << std::endl;

  std::cout << "\n>>> WYNIKI PORÓWNAWCZE ALIGNMENTU (X0 Offset), Regeneracja:" << std::endl;
  std::cout << ">>> [DC] Metoda tradycyjna 2D: (" << x0_DC_2D << " +- " << x0_DC_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << x0_prof_DC << " +- " << x0_prof_DC_err << ") cm" << std::endl;
  std::cout << ">>> [BP] Metoda tradycyjna 2D: (" << x0_BP_2D << " +- " << x0_BP_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << x0_prof_BP << " +- " << x0_prof_BP_err << ") cm" << std::endl;

  std::cout << "\n>>> WYNIKI PORÓWNAWCZE ALIGNMENTU (R), Regeneracja:" << std::endl;
  std::cout << ">>> [DC] Metoda tradycyjna 2D: (" << r_DC_2D << " +- " << r_DC_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << r_DC_prof << " +- " << r_DC_prof_err << ") cm" << std::endl;
  std::cout << ">>> [BP] Metoda tradycyjna 2D: (" << r_BP_2D << " +- " << r_BP_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << r_BP_prof << " +- " << r_BP_prof_err << ") cm" << std::endl;

  // Find local maxima for cylindrical radius histograms
  channel = "Data";
  // 1. Stare wywołanie dopasowania 2D na histogramach (pozostawione dla porównania grafiki)
  TFitResultPtr fit_result_2D_DC_Data = KloeTools::fit_method_circle(h_xy_charged_DC[channel], -50, 50);
  TFitResultPtr fit_result_2D_BP_Data = KloeTools::fit_method_circle(h_xy_charged_BP[channel], -50, 50);

  x0_DC_2D = fit_result_2D_DC_Data->IsValid() ? fit_result_2D_DC_Data->Parameter(0) : 0.0;
  x0_BP_2D = fit_result_2D_BP_Data->IsValid() ? fit_result_2D_BP_Data->Parameter(0) : 0.0;
  y0_DC_2D = fit_result_2D_DC_Data->IsValid() ? fit_result_2D_DC_Data->Parameter(1) : 0.0;
  y0_BP_2D = fit_result_2D_BP_Data->IsValid() ? fit_result_2D_BP_Data->Parameter(1) : 0.0;
  r_DC_2D = fit_result_2D_DC_Data->IsValid() ? fit_result_2D_DC_Data->Parameter(2) : 0.0;
  r_BP_2D = fit_result_2D_BP_Data->IsValid() ? fit_result_2D_BP_Data->Parameter(2) : 0.0;

  x0_DC_2D_err = fit_result_2D_DC_Data->IsValid() ? fit_result_2D_DC_Data->ParError(0) : 0.0;
  x0_BP_2D_err = fit_result_2D_BP_Data->IsValid() ? fit_result_2D_BP_Data->ParError(0) : 0.0;
  y0_DC_2D_err = fit_result_2D_DC_Data->IsValid() ? fit_result_2D_DC_Data->ParError(1) : 0.0;
  y0_BP_2D_err = fit_result_2D_BP_Data->IsValid() ? fit_result_2D_BP_Data->ParError(1) : 0.0;
  r_DC_2D_err = fit_result_2D_DC_Data->IsValid() ? fit_result_2D_DC_Data->ParError(2) : 0.0;
  r_BP_2D_err = fit_result_2D_BP_Data->IsValid() ? fit_result_2D_BP_Data->ParError(2) : 0.0;

  // 2. NOWE, STABILNE WYWOŁANIE: Dopasowanie sinusoidy na profilach polarnych
  // Przekazujemy promień spodziewany jako punkt startowy dla MINUITa (26.0 cm dla DC oraz 4.5 cm dla BP)
  TFitResultPtr fit_result_DC_Data = KloeTools::fit_method_profile(h_prof_polar_DC[channel], 26.0);
  TFitResultPtr fit_result_BP_Data = KloeTools::fit_method_profile(h_prof_polar_BP[channel], 4.5);

  r_DC_prof = fit_result_DC_Data->IsValid() ? fit_result_DC_Data->Parameter(0) : 0.0;
  r_BP_prof = fit_result_BP_Data->IsValid() ? fit_result_BP_Data->Parameter(0) : 0.0;
  r_DC_prof_err = fit_result_DC_Data->IsValid() ? fit_result_DC_Data->ParError(0) : 0.0;
  r_BP_prof_err = fit_result_BP_Data->IsValid() ? fit_result_BP_Data->ParError(0) : 0.0;

  x0_prof_DC = fit_result_DC_Data->IsValid() ? fit_result_DC_Data->Parameter(1) : 0.0;
  x0_prof_BP = fit_result_BP_Data->IsValid() ? fit_result_BP_Data->Parameter(1) : 0.0;
  x0_prof_DC_err = fit_result_DC_Data->IsValid() ? fit_result_DC_Data->ParError(1) : 0.0;
  x0_prof_BP_err = fit_result_BP_Data->IsValid() ? fit_result_BP_Data->ParError(1) : 0.0;

  y0_prof_DC = fit_result_DC_Data->IsValid() ? fit_result_DC_Data->Parameter(2) : 0.0;
  y0_prof_BP = fit_result_BP_Data->IsValid() ? fit_result_BP_Data->Parameter(2) : 0.0;
  y0_prof_DC_err = fit_result_DC_Data->IsValid() ? fit_result_DC_Data->ParError(2) : 0.0;
  y0_prof_BP_err = fit_result_BP_Data->IsValid() ? fit_result_BP_Data->ParError(2) : 0.0;

  std::cout << "\n>>> WYNIKI PORÓWNAWCZE ALIGNMENTU (Y0 Offset), Data:" << std::endl;
  std::cout << ">>> [DC] Metoda tradycyjna 2D: (" << y0_DC_2D << " +- " << y0_DC_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << y0_prof_DC << " +- " << y0_prof_DC_err << ") cm" << std::endl;
  std::cout << ">>> [BP] Metoda tradycyjna 2D: (" << y0_BP_2D << " +- " << y0_BP_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << y0_prof_BP << " +- " << y0_prof_BP_err << ") cm" << std::endl;

  std::cout << "\n>>> WYNIKI PORÓWNAWCZE ALIGNMENTU (X0 Offset), Data:" << std::endl;
  std::cout << ">>> [DC] Metoda tradycyjna 2D: (" << x0_DC_2D << " +- " << x0_DC_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << x0_prof_DC << " +- " << x0_prof_DC_err << ") cm" << std::endl;
  std::cout << ">>> [BP] Metoda tradycyjna 2D: (" << x0_BP_2D << " +- " << x0_BP_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << x0_prof_BP << " +- " << x0_prof_BP_err << ") cm" << std::endl;

  std::cout << "\n>>> WYNIKI PORÓWNAWCZE ALIGNMENTU (R), Data:" << std::endl;
  std::cout << ">>> [DC] Metoda tradycyjna 2D: (" << r_DC_2D << " +- " << r_DC_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << r_DC_prof << " +- " << r_DC_prof_err << ") cm" << std::endl;
  std::cout << ">>> [BP] Metoda tradycyjna 2D: (" << r_BP_2D << " +- " << r_BP_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << r_BP_prof << " +- " << r_BP_prof_err << ") cm" << std::endl;

  // Find local maxima for cylindrical radius histograms
  channel = "MC sum";
  // 1. Stare wywołanie dopasowania 2D na histogramach (pozostawione dla porównania grafiki)
  TFitResultPtr fit_result_2D_DC_MCSum = KloeTools::fit_method_circle(h_xy_charged_DC[channel], -50, 50);
  TFitResultPtr fit_result_2D_BP_MCSum = KloeTools::fit_method_circle(h_xy_charged_BP[channel], -50, 50);

  x0_DC_2D = fit_result_2D_DC_MCSum->IsValid() ? fit_result_2D_DC_MCSum->Parameter(0) : 0.0;
  x0_BP_2D = fit_result_2D_BP_MCSum->IsValid() ? fit_result_2D_BP_MCSum->Parameter(0) : 0.0;
  y0_DC_2D = fit_result_2D_DC_MCSum->IsValid() ? fit_result_2D_DC_MCSum->Parameter(1) : 0.0;
  y0_BP_2D = fit_result_2D_BP_MCSum->IsValid() ? fit_result_2D_BP_MCSum->Parameter(1) : 0.0;
  r_DC_2D = fit_result_2D_DC_MCSum->IsValid() ? fit_result_2D_DC_MCSum->Parameter(2) : 0.0;
  r_BP_2D = fit_result_2D_BP_MCSum->IsValid() ? fit_result_2D_BP_MCSum->Parameter(2) : 0.0;

  x0_DC_2D_err = fit_result_2D_DC_MCSum->IsValid() ? fit_result_2D_DC_MCSum->ParError(0) : 0.0;
  x0_BP_2D_err = fit_result_2D_BP_MCSum->IsValid() ? fit_result_2D_BP_MCSum->ParError(0) : 0.0;
  y0_DC_2D_err = fit_result_2D_DC_MCSum->IsValid() ? fit_result_2D_DC_MCSum->ParError(1) : 0.0;
  y0_BP_2D_err = fit_result_2D_BP_MCSum->IsValid() ? fit_result_2D_BP_MCSum->ParError(1) : 0.0;
  r_DC_2D_err = fit_result_2D_DC_MCSum->IsValid() ? fit_result_2D_DC_MCSum->ParError(2) : 0.0;
  r_BP_2D_err = fit_result_2D_BP_MCSum->IsValid() ? fit_result_2D_BP_MCSum->ParError(2) : 0.0;

  // 2. NOWE, STABILNE WYWOŁANIE: Dopasowanie sinusoidy na profilach polarnych
  // Przekazujemy promień spodziewany jako punkt startowy dla MINUITa (26.0 cm dla DC oraz 4.5 cm dla BP)
  TFitResultPtr fit_result_DC_MCSum = KloeTools::fit_method_profile(h_prof_polar_DC[channel], 26.0);
  TFitResultPtr fit_result_BP_MCSum = KloeTools::fit_method_profile(h_prof_polar_BP[channel], 4.5);

  r_DC_prof = fit_result_DC_MCSum->IsValid() ? fit_result_DC_MCSum->Parameter(0) : 0.0;
  r_BP_prof = fit_result_BP_MCSum->IsValid() ? fit_result_BP_MCSum->Parameter(0) : 0.0;
  r_DC_prof_err = fit_result_DC_MCSum->IsValid() ? fit_result_DC_MCSum->ParError(0) : 0.0;
  r_BP_prof_err = fit_result_BP_MCSum->IsValid() ? fit_result_BP_MCSum->ParError(0) : 0.0;

  x0_prof_DC = fit_result_DC_MCSum->IsValid() ? fit_result_DC_MCSum->Parameter(1) : 0.0;
  x0_prof_BP = fit_result_BP_MCSum->IsValid() ? fit_result_BP_MCSum->Parameter(1) : 0.0;
  x0_prof_DC_err = fit_result_DC_MCSum->IsValid() ? fit_result_DC_MCSum->ParError(1) : 0.0;
  x0_prof_BP_err = fit_result_BP_MCSum->IsValid() ? fit_result_BP_MCSum->ParError(1) : 0.0;

  y0_prof_DC = fit_result_DC_MCSum->IsValid() ? fit_result_DC_MCSum->Parameter(2) : 0.0;
  y0_prof_BP = fit_result_BP_MCSum->IsValid() ? fit_result_BP_MCSum->Parameter(2) : 0.0;
  y0_prof_DC_err = fit_result_DC_MCSum->IsValid() ? fit_result_DC_MCSum->ParError(2) : 0.0;
  y0_prof_BP_err = fit_result_BP_MCSum->IsValid() ? fit_result_BP_MCSum->ParError(2) : 0.0;

  std::cout << "\n>>> WYNIKI PORÓWNAWCZE ALIGNMENTU (Y0 Offset), Data:" << std::endl;
  std::cout << ">>> [DC] Metoda tradycyjna 2D: (" << y0_DC_2D << " +- " << y0_DC_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << y0_prof_DC << " +- " << y0_prof_DC_err << ") cm" << std::endl;
  std::cout << ">>> [BP] Metoda tradycyjna 2D: (" << y0_BP_2D << " +- " << y0_BP_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << y0_prof_BP << " +- " << y0_prof_BP_err << ") cm" << std::endl;

  std::cout << "\n>>> WYNIKI PORÓWNAWCZE ALIGNMENTU (X0 Offset), Data:" << std::endl;
  std::cout << ">>> [DC] Metoda tradycyjna 2D: (" << x0_DC_2D << " +- " << x0_DC_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << x0_prof_DC << " +- " << x0_prof_DC_err << ") cm" << std::endl;
  std::cout << ">>> [BP] Metoda tradycyjna 2D: (" << x0_BP_2D << " +- " << x0_BP_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << x0_prof_BP << " +- " << x0_prof_BP_err << ") cm" << std::endl;

  std::cout << "\n>>> WYNIKI PORÓWNAWCZE ALIGNMENTU (R), Data:" << std::endl;
  std::cout << ">>> [DC] Metoda tradycyjna 2D: (" << r_DC_2D << " +- " << r_DC_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << r_DC_prof << " +- " << r_DC_prof_err << ") cm" << std::endl;
  std::cout << ">>> [BP] Metoda tradycyjna 2D: (" << r_BP_2D << " +- " << r_BP_2D_err << ") cm  |  Metoda profilu polarnego 1D: (" << r_BP_prof << " +- " << r_BP_prof_err << ") cm" << std::endl;

  // Zapis wyników do pliku ROOT
  TFile output_file("regeneration_cross_section_v2.root", "RECREATE");
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
    
    // 1. RYSOWANIE HISTOGRAMU 2D DC WRAZ Z DOPASOWANĄ SINUSOIDĄ
    h_polar_DC[channel_name]->Draw("COLZ");
    
    // Wyciągamy funkcję dopasowania z powiązanego profilu 1D
    // Domyślnie ROOT nazywa funkcję dopasowania w profilu jako "prevfit_" lub nazwą użytą w fit_method_profile.
    // Najbezpieczniej pobrać ją po liście funkcji profilu:
    TF1* fit_func_DC = nullptr;
    if (h_prof_polar_DC[channel_name]->GetListOfFunctions()->GetEntries() > 0) {
        fit_func_DC = (TF1*)h_prof_polar_DC[channel_name]->GetListOfFunctions()->At(0);
    }

    if (fit_func_DC) {
        fit_func_DC->SetLineColor(kBlack);      // Czerwona linia, idealnie widoczna na COLZ
        fit_func_DC->SetLineWidth(4);         // Grubsza linia (szerokość 4 piksele)
        fit_func_DC->SetLineStyle(1);         // Linia ciągła
        fit_func_DC->Draw("SAME");            // <--- KLUCZOWE: Nakładamy funkcję na rozkład 2D
    }
    c->SaveAs(("img/polar_DC_" + channel_name + ".png").c_str());
    
    // 2. RYSOWANIE HISTOGRAMU 2D BP WRAZ Z DOPASOWANĄ SINUSOIDĄ
    c->Clear(); // Czyścimy canvas przed kolejnym rysowaniem
    h_polar_BP[channel_name]->Draw("COLZ");
    
    TF1* fit_func_BP = nullptr;
    if (h_prof_polar_BP[channel_name]->GetListOfFunctions()->GetEntries() > 0) {
        fit_func_BP = (TF1*)h_prof_polar_BP[channel_name]->GetListOfFunctions()->At(0);
    }

    if (fit_func_BP) {
        fit_func_BP->SetLineColor(kRed);
        fit_func_BP->SetLineWidth(4);
        fit_func_BP->Draw("SAME");
    }
    c->SaveAs(("img/polar_BP_" + channel_name + ".png").c_str());
    
    // 3. POZOSTAŁE HISTOGRAMY (Rysujesz tak jak wcześniej, dodaj tylko c->Clear() dla porządku)
    c->Clear();
    h_prof_polar_DC[channel_name]->Draw();
    c->SaveAs(("img/prof_polar_DC_" + channel_name + ".png").c_str());
    
    c->Clear();
    h_prof_polar_BP[channel_name]->Draw();
    c->SaveAs(("img/prof_polar_BP_" + channel_name + ".png").c_str());
    
    c->Clear();
    h_xy_charged_DC[channel_name]->Draw("COLZ");
    c->SaveAs(("img/xy_charged_DC_" + channel_name + ".png").c_str());
    
    c->Clear();
    h_xy_charged_BP[channel_name]->Draw("COLZ");
    c->SaveAs(("img/xy_charged_BP_" + channel_name + ".png").c_str());
    
    delete c;
    // --- KONIEC BLOKU ---
  }

  output_file.Close();


}
