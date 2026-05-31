#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TChain.h>

#include <string>
#include <array>

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

gROOT->SetBatch(kTRUE);

std::vector<double> SmartFit(TH1F *h, TString title, double peakPos, double lowOff, double highOff, double lowSideband, double highSideband, double rms, int padIdx, TCanvas *c, bool useExpo = false, bool useLangau = false, bool useAsym = false, bool useCrystalBall = false, bool useTwoPass = false)
{
  c->cd(padIdx);
  h->SetTitle(title);

  double xMin = peakPos - lowOff;
  double xMax = peakPos + highOff;
  double sidebandMin = peakPos - lowSideband;
  double sidebandMax = peakPos + highSideband;

  double drawMin = xMin - 1.0;
  double drawMax = xMax + 1.0;

  double max_val = h->GetMaximum();
  h->GetYaxis()->SetRangeUser(0.0, max_val * 1.5);

  h->GetXaxis()->SetRangeUser(drawMin, drawMax);

  TString sigFormula;
  if (useAsym)
  {
    sigFormula = "[0]*exp(-0.5*((x-[1])/(x<[1]?[2]:[5]))^2)";
  }
  else if (useLangau)
  {
    sigFormula = "TMath::Landau(x, [1], [2], 1) * [0]";
  }
  else if (useCrystalBall)
  {
    sigFormula = ""; // nieużywane - CB używa lambdy
  }
  else
  {
    sigFormula = "[0]*exp(-0.5*((x-[1])/[2])^2)";
  }

  // Crystal Ball ma 5 parametrow sygnalu (0-4), wiec tlo zaczyna od 5; pozostale od 3
  int bgStart = useCrystalBall ? 5 : 3;

  // Wspolna lambda CB (sygnal + tlo) -- capture by value
  auto cb_signal = [](double x, double amp, double mean, double sigma, double alpha, double n) -> double
  {
    double t = (x - mean) / sigma;
    if (t > -std::abs(alpha))
      return amp * std::exp(-0.5 * t * t);
    double absA = std::abs(alpha);
    double A = std::pow(n / absA, n) * std::exp(-0.5 * absA * absA);
    double B = n / absA - absA;
    return amp * A / std::pow(B - t, n);
  };

  TF1 *f_fit;
  if (useCrystalBall)
  {
    f_fit = new TF1("f_fit", [cb_signal, useExpo](double *x, double *p) -> double
    {
      double sig = cb_signal(x[0], p[0], p[1], p[2], p[3], p[4]);
      double bg  = useExpo ? p[5] * std::exp(p[6] * x[0]) : p[5] + p[6] * x[0];
      return sig + bg;
    }, sidebandMin, sidebandMax, 7);
  }
  else
  {
    TString bgFormula = useExpo ? TString(Form("expo(%d)", bgStart)) : TString(Form("pol1(%d)", bgStart));
    f_fit = new TF1("f_fit", sigFormula + " + " + bgFormula, sidebandMin, sidebandMax);
  }

  if (useCrystalBall)
  {
    if (useExpo)
      f_fit->SetParameters(max_val * 0.8, peakPos, 0.1 * rms, -1.5, 5.0, TMath::Log(max_val * 0.1 + 1e-9), -0.05);
    else
      f_fit->SetParameters(max_val * 0.8, peakPos, 0.1 * rms, -1.5, 5.0, max_val * 0.1, 0.0);
    f_fit->SetParLimits(3, -10.0, -0.1);  // alpha > 0
    f_fit->SetParLimits(4, 1.0, 50.0);  // n > 1
    f_fit->SetParLimits(bgStart,     0.0, 2.0 * max_val);
    f_fit->SetParLimits(bgStart + 1, -10.0, 0.0);
  }
  else if (useExpo)
  {
    f_fit->SetParameters(max_val * 0.8, peakPos, 0.1 * rms, TMath::Log(max_val * 0.1 + 1e-9), -0.05);
    f_fit->SetParLimits(4, -10.0, 0.0); // Stabilizacja: tło wykładnicze nie może rosnąć
  }
  else
  {
    f_fit->SetParameters(max_val * 0.8, peakPos, 0.1 * rms, max_val * 0.1, 0.0);
    f_fit->SetParLimits(3, 0.0, max_val); // Stabilizacja: tło nie może być ujemne
    f_fit->SetParLimits(4, -10.0, 0.0); // Stabilizacja: współczynnik kierunkowy tła nie może być zbyt duży
  }

  if (useAsym)
  {
    f_fit->SetParameter(5, 0.1 * rms);
    f_fit->SetParLimits(5, 0.01 * rms, 1.5 * rms);
  }

  f_fit->SetParLimits(0, 0.01 * max_val, 10 * max_val);
  f_fit->SetParLimits(1, peakPos - 0.5, peakPos + 0.5); // Zacieśnienie limitu pozycji piku
  f_fit->SetParLimits(2, 0.01 * rms, 0.5 * rms);        // Zacieśnienie limitu szerokości

  if (useTwoPass)
  {
    // ---- Pierwsza przepustka: tylko tlo, dopasowywane na skrzydlach (bez okolicy piku) ----
    TF1 *f_bg1 = new TF1("f_bg1", useExpo ? "expo" : "pol1", sidebandMin, sidebandMax);
    if (useExpo)
      f_bg1->SetParameters(TMath::Log(max_val * 0.1 + 1e-9), -0.05);
    else
      f_bg1->SetParameters(max_val * 0.1, 0.0);

    // Klon histogramu z wyzerowanymi binami w okolicach piku (error=0 => bin pomijany w chi2)
    TH1F *h_sb = (TH1F *)h->Clone("h_sb");
    // double exclLow  = peakPos - 2.0 * rms;
    // double exclHigh = peakPos + 2.0 * rms;
    // for (int ib = 1; ib <= h_sb->GetNbinsX(); ib++)
    // {
    //   double bc = h_sb->GetBinCenter(ib);
    //   if (bc >= exclLow && bc <= exclHigh)
    //   {
    //     h_sb->SetBinContent(ib, 0);
    //     h_sb->SetBinError(ib, 0);
    //   }
    // }
    h_sb->Fit(f_bg1, "RQ0");
    delete h_sb;

    // ---- Zafiksuj parametry tla z pierwszej przepustki ----
    f_fit->FixParameter(bgStart,     f_bg1->GetParameter(0));
    f_fit->FixParameter(bgStart + 1, f_bg1->GetParameter(1));
    delete f_bg1;
  }

  // ---- Druga przepustka (lub jedyna, gdy !useTwoPass): fituje sie sygnal (+tlo jesli niezafiksowane) ----
  // "I" uzywa calki z funkcji w binie (kluczowe przy waskich pikach i malym binowaniu)
  // "M" szuka lepszych minimow, co stabilizuje macierz kowariancji
  TFitResultPtr r = h->Fit(f_fit, "RMESLQ");

  h->Draw("E1");

  f_fit->GetXaxis()->SetRangeUser(drawMin, drawMax);

  TF1 *sigOnly;
  if (useCrystalBall)
  {
    sigOnly = new TF1("sigOnly", [cb_signal](double *x, double *p) -> double
    {
      return cb_signal(x[0], p[0], p[1], p[2], p[3], p[4]);
    }, xMin, xMax, 5);
    sigOnly->SetParameter(0, f_fit->GetParameter(0));
    sigOnly->SetParameter(1, f_fit->GetParameter(1));
    sigOnly->SetParameter(2, f_fit->GetParameter(2));
    sigOnly->SetParameter(3, f_fit->GetParameter(3));
    sigOnly->SetParameter(4, f_fit->GetParameter(4));
  }
  else
  {
    sigOnly = new TF1("sigOnly", sigFormula, xMin, xMax);
    sigOnly->SetParameter(0, f_fit->GetParameter(0));
    sigOnly->SetParameter(1, f_fit->GetParameter(1));
    sigOnly->SetParameter(2, f_fit->GetParameter(2));
    if (useAsym)
      sigOnly->SetParameter(5, f_fit->GetParameter(5));
  }
  sigOnly->SetParErrors(f_fit->GetParErrors());

  double peakMean = sigOnly->Mean(xMin, xMax);
  double variance = sigOnly->CentralMoment(2, xMin, xMax);
  double combinedSigma = (variance > 0) ? std::sqrt(variance) : f_fit->GetParameter(2);

  double peakMeanErr = f_fit->GetParError(1);
  double combinedSigmaErr = f_fit->GetParError(2);

  double low = peakMean - 3 * combinedSigma;
  double high = peakMean + 3 * combinedSigma;
  double binW = h->GetBinWidth(1);

  double Yield = (sigOnly->Integral(low, high)) / binW;
  double YieldErr = 0;

  if (r.Get())
  {
    TF1 *sigFuncForErr = (TF1 *)f_fit->Clone("sigFuncForErr");
    // ZMIANA 2: Ustawienie parametrów tła na ZERO, ale bez usuwania ich z definicji
    // Pozwala to IntegralError poprawnie odczytać korelacje z macierzy kowariancji r
    sigFuncForErr->FixParameter(bgStart,     0);
    sigFuncForErr->FixParameter(bgStart + 1, 0);

    // ZMIANA 3: Obliczanie błędu z uwzględnieniem pełnej macierzy kowariancji fita
    YieldErr = sigFuncForErr->IntegralError(low, high,
                                            r->GetParams(),
                                            r->GetCovarianceMatrix().GetMatrixArray()) /
               binW;

    delete sigFuncForErr;
  }

  TF1 *bg = new TF1("bg", useExpo ? "expo" : "pol1", xMin, xMax);
  bg->SetParameters(f_fit->GetParameter(bgStart), f_fit->GetParameter(bgStart + 1));
  bg->SetLineStyle(3);
  bg->SetLineColor(kBlue);
  bg->Draw("same");

  TLegend *l = new TLegend(0.7, 0.65, 0.85, 0.88);
  l->SetTextSize(0.035);
  l->AddEntry(f_fit, "Fit", "l");
  l->AddEntry(bg, "Background", "l");
  l->AddEntry((TObject *)0, Form("Yield: %.0f #pm %.0f", Yield, YieldErr), "");
  l->AddEntry((TObject *)0, Form("#chi^{2}/ndf: %.2f", f_fit->GetChisquare() / f_fit->GetNDF()), "");
  l->AddEntry((TObject *)0, Form("#mu: %.2f", peakMean), "");
  l->AddEntry((TObject *)0, Form("RMS: %.2f", combinedSigma), "");
  l->Draw();

  TLine *line_left = new TLine(low, 0, low, h->GetMaximum());
  line_left->SetLineStyle(2);
  line_left->SetLineColor(kBlack);
  line_left->Draw();
  TLine *line_right = new TLine(high, 0, high, h->GetMaximum());
  line_right->SetLineStyle(2);
  line_right->SetLineColor(kBlack);
  line_right->Draw();
  TLine *line_mean = new TLine(peakMean, 0, peakMean, h->GetMaximum());
  line_mean->SetLineStyle(2);
  line_mean->SetLineColor(kBlack);
  line_mean->Draw();

  return {xMin, xMax, peakMean, combinedSigma, combinedSigmaErr};
}

void DrawRegeneration(TH1F *h, double xMin, double xMax, double peakMean, double peakSigma, int padIdx, TCanvas *c)
{
  c->cd(padIdx);
  h->GetXaxis()->SetRangeUser(xMin - 1.0, xMax + 1.0);
  h->SetFillColor(kGreen - 9);
  h->SetTitle("Pure Regeneration Component (Integral Range)");
  h->Draw("HIST");

  TLine *line_left = new TLine(peakMean - 3 * peakSigma, 0, peakMean - 3 * peakSigma, h->GetMaximum());
  line_left->SetLineStyle(2);
  line_left->SetLineColor(kBlack);
  line_left->Draw();

  TLine *line_right = new TLine(peakMean + 3 * peakSigma, 0, peakMean + 3 * peakSigma, h->GetMaximum());
  line_right->SetLineStyle(2);
  line_right->SetLineColor(kBlack);
  line_right->Draw();

  TLine *line_mean = new TLine(peakMean, 0, peakMean, h->GetMaximum());
  line_mean->SetLineStyle(2);
  line_mean->SetLineColor(kBlack);
  line_mean->Draw();

  // Zliczanie
  double yield = h->Integral(h->FindBin(peakMean - 3 * peakSigma), h->FindBin(peakMean + 3 * peakSigma));

  TLegend *l = new TLegend(0.15, 0.75, 0.4, 0.88);
  l->AddEntry(h, Form("Reg. Counts: %.1f", yield), "f");
  l->Draw();
}

void regeneration_cross_section(TString root_file_date)
{
  KLOE::setGlobalStyle();
  gStyle->SetOptStat(0);

  TChain *chain = new TChain("h1");

  TString base_path = "/home/g/gamrat/root_files/" + root_file_date + "/Signal/";
  std::array<TString, 3> mc_paths = {"ALL_PHYS_SIGNAL_NoSmearing/", "ALL_PHYS2_SIGNAL_NoSmearing/", "ALL_PHYS3_SIGNAL_NoSmearing/"};
  TString data_path = base_path + "DATA_SIGNAL_NoSmearing/*.root";

  for (const auto &path : mc_paths)
  {
    chain->Add(base_path + path + "*.root");
  }
  chain->Add(data_path);

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
  std::map<TString, TH1 *> h_rho_charged;
  std::map<TString, TH1 *> h_rho_neutral;

  std::map<TString, TH1 *> h_R_charged;
  std::map<TString, TH1 *> h_R_neutral;

  std::map<TString, TH2 *> h_rho_charged_vs_R;
  std::map<TString, TH2 *> h_rho_neutral_vs_R;

  for (const auto &chann : KLOE::channName)
  {
    h_rho_charged[chann.second] = new TH1F("h_rho_charged_" + chann.second, "Rho Charged " + chann.second + ";#rho_{charged};Entries", 251, 0, 50);
    h_rho_neutral[chann.second] = new TH1F("h_rho_neutral_" + chann.second, "Rho Neutral " + chann.second + ";#rho_{neutral};Entries", 251, 0, 50);

    h_R_charged[chann.second] = new TH1F("h_R_charged_" + chann.second, "R Charged " + chann.second + ";R_{charged};Entries", 251, 0, 50);
    h_R_neutral[chann.second] = new TH1F("h_R_neutral_" + chann.second, "R Neutral " + chann.second + ";R_{neutral};Entries", 251, 0, 50);

    h_rho_charged_vs_R[chann.second] = new TH2F("h_rho_charged_vs_R_" + chann.second, "Rho Charged vs R " + chann.second + ";R_{charged};#rho_{charged}", 251, 0, 50, 251, 0, 50);
    h_rho_neutral_vs_R[chann.second] = new TH2F("h_rho_neutral_vs_R_" + chann.second, "Rho Neutral vs R " + chann.second + ";R_{neutral};#rho_{neutral}", 251, 0, 50, 251, 0, 50);
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

  TVector3 kch_position(0, 0, 0),
      kne_position(0, 0, 0),
      ip_position(0, 0, 0);

  std::map<TString, Int_t> entries_count; // Map to count entries for each channel

  while (reader.Next())
  {
    // Cut setting
    std::array<Double_t, 3> distNeutralCharged = {KchrecFit[6] - KnerecFit[6],
                                                  KchrecFit[7] - KnerecFit[7],
                                                  KchrecFit[8] - KnerecFit[8]},
                            distNeutralIP = {Knerec[6] - *Bx,
                                             Knerec[7] - *By,
                                             Knerec[8] - KchrecClosest[8]},
                            distChargedIP = {KchrecClosest[6] - *Bx,
                                             KchrecClosest[7] - *By,
                                             KchrecClosest[8] - *Bz};

    Float_t
        rho_pm = std::sqrt(distChargedIP[0] * distChargedIP[0] + distChargedIP[1] * distChargedIP[1]),
        rho_00 = std::sqrt(distNeutralIP[0] * distNeutralIP[0] + distNeutralIP[1] * distNeutralIP[1]),
        rho = std::sqrt(std::pow(rho_pm, 2) + std::pow(rho_00, 2));

    Double_t radius00 = std::sqrt(std::pow(Knerec[6] - *Bx, 2) +
                                  std::pow(Knerec[7] - *By, 2)),
             radiuspm = std::sqrt(std::pow(KchrecClosest[6] - *Bx, 2) +
                                  std::pow(KchrecClosest[7] - *By, 2)),
             zdist00 = std::abs(Knerec[8] - KchrecClosest[8]),
             zdistpm = std::abs(KchrecClosest[8] - *Bz);

    Double_t fiducialVolume = std::sqrt(std::pow(distNeutralCharged[0], 2) + std::pow(distNeutralCharged[1], 2)) < 2.05 && std::abs(distNeutralCharged[2]) < 2.45,
             fiducialVolumeClose = radius00 < 1.5 && radiuspm < 2.0 && zdist00 < 1.5 && zdistpm < 1.5;

    TVector3 KchrecVec = {KchrecFit[0], KchrecFit[1], KchrecFit[2]};
    TVector3 trk1VecFit = {trk1Fit[0], trk1Fit[1], trk1Fit[2]};
    TVector3 trk2VecFit = {trk2Fit[0], trk2Fit[1], trk2Fit[2]};

    Double_t phiTrk1Angle = cos(trk1VecFit.Angle(KchrecVec)),
             phiTrk2Angle = cos(trk2VecFit.Angle(KchrecVec));

    Bool_t global_cut = ((fiducialVolume && (abs(phiTrk1Angle) < 0.8 || abs(phiTrk2Angle) < 0.8)) || !fiducialVolume) && ((fiducialVolumeClose && rho > 1.5) || !fiducialVolumeClose) && ellipse_cut(*minv4gam - PhysicsConstants::mK0, KnerecSix[5] - PhysicsConstants::mK0);

    // ---

    kch_position.SetXYZ(Kchrec[6], Kchrec[7], Kchrec[8]);
    kne_position.SetXYZ(Knerec[6], Knerec[7], Knerec[8]);
    ip_position.SetXYZ(ip[0], ip[1], ip[2]);

    kch_position -= ip_position;
    kne_position -= ip_position;

    Double_t R_charged = kch_position.Mag(),
             R_neutral = kne_position.Mag(),
             rho_charged = kch_position.Perp(),
             rho_neutral = kne_position.Perp();

    Double_t mean_BP_Ch_MC = 10.44, sigma_BP_Ch_MC = 0.69,
             mean_BP_Ch_Data = 10.50, sigma_BP_Ch_Data = 0.76,
             mean_BP_Ne_MC = 10.30, sigma_BP_Ne_MC = 1.11,
             mean_BP_Ne_Data = 10.44, sigma_BP_Ne_Data = 1.05;

    if (*mcflag == 1 && global_cut) // MC
    {
      if (*mctruth == -1 || *mctruth == 0)
        continue; // Skip events with mctruth -1 or 0

      if (*time_ne_CM_signal_fit < 7.0)
      {
        h_rho_charged[KLOE::channName.at(*mctruth)]->Fill(rho_charged);
        h_R_charged[KLOE::channName.at(*mctruth)]->Fill(R_charged);

        h_rho_charged_vs_R[KLOE::channName.at(*mctruth)]->Fill(R_charged, rho_charged);
      }

      if (*time_ch_CM_signal_fit < 7.0)
      {
        h_rho_neutral[KLOE::channName.at(*mctruth)]->Fill(rho_neutral);
        h_R_neutral[KLOE::channName.at(*mctruth)]->Fill(R_neutral);

        h_rho_neutral_vs_R[KLOE::channName.at(*mctruth)]->Fill(R_neutral, rho_neutral);
      }

      entries_count[KLOE::channName.at(*mctruth)]++; // Increment count for this channel
    }

    if (*mcflag == 0 && global_cut) // Data
    {
      if (*time_ne_CM_signal_fit < 7.0)
      {
        h_rho_charged["Data"]->Fill(rho_charged);
        h_R_charged["Data"]->Fill(R_charged);
        h_rho_charged_vs_R["Data"]->Fill(R_charged, rho_charged);
      }

      if (*time_ch_CM_signal_fit < 7.0)
      {
        h_rho_neutral["Data"]->Fill(rho_neutral);
        h_R_neutral["Data"]->Fill(R_neutral);
        h_rho_neutral_vs_R["Data"]->Fill(R_neutral, rho_neutral);
      }

      entries_count["Data"]++; // Increment count for data
    }
  }

  Long_t totalMC_entries = 0;

  for (const auto &chann : KLOE::channName)
  {
    if (chann.second == "Data" || chann.second == "MC sum")
      continue; // Skip data and sum itself

    std::cout << "Entries for " << chann.second << ": " << entries_count[chann.second] << std::endl;
    totalMC_entries += entries_count[chann.second];

    std::cout << "Total MC entries: " << totalMC_entries << std::endl;
    std::cout << "Entries for Data: " << entries_count["Data"] << std::endl;
  }

  Double_t scale_factor = entries_count["Data"] / static_cast<Double_t>(totalMC_entries);

  // --- MC sum ---
  for (const auto &chann : KLOE::channName)
  {
    if (chann.second == "Data" || chann.second == "MC sum")
      continue; // Skip data and sum itself

    h_rho_charged["MC sum"]->Add(h_rho_charged[chann.second], scale_factor);
    h_rho_neutral["MC sum"]->Add(h_rho_neutral[chann.second], scale_factor);

    h_R_charged["MC sum"]->Add(h_R_charged[chann.second], scale_factor);
    h_R_neutral["MC sum"]->Add(h_R_neutral[chann.second], scale_factor);
    h_rho_charged_vs_R["MC sum"]->Add(h_rho_charged_vs_R[chann.second], scale_factor);
    h_rho_neutral_vs_R["MC sum"]->Add(h_rho_neutral_vs_R[chann.second], scale_factor);
  };

  // --- SEKCJA RYSOWANIA ---
  TCanvas *c1 = new TCanvas("c1", "Analiza", 800, 1000);
  c1->Divide(1, 3);
  TString outName = "regeneration_results_" + root_file_date + ".pdf";
  c1->Print(outName + "[");

  struct PlotGroup
  {
    TString name;
    TH1F *data;
    TH1F *mc;
    TH1F *reg;
    double pos;
    double lowOff;
    double highOff;
    double lowSideband;
    double highSideband;
    TF1 *weight_func = nullptr;
    double rangeLeft  = 0.0;
    double rangeRight = 0.0;
  };
  std::vector<PlotGroup> groups = {
      {"R_Charged", (TH1F *)h_R_charged["Data"], (TH1F *)h_R_charged["MC sum"], (TH1F *)h_R_charged["Regeneration"], 10.0, 3.0, 3.0, 9.0, 12.0},
      {"Rho_Charged", (TH1F *)h_rho_charged["Data"], (TH1F *)h_rho_charged["MC sum"], (TH1F *)h_rho_charged["Regeneration"], 25.0, 4.0, 4.0, 15.0, 15.0},
      {"R_Neutral", (TH1F *)h_R_neutral["Data"], (TH1F *)h_R_neutral["MC sum"], (TH1F *)h_R_neutral["Regeneration"], 10.0, 3.0, 3.0, 9.0, 12.0},
      {"Rho_Neutral", (TH1F *)h_rho_neutral["Data"], (TH1F *)h_rho_neutral["MC sum"], (TH1F *)h_rho_neutral["Regeneration"], 25.0, 4.0, 4.0, 15.0, 15.0}};

  TFile *file = new TFile("regeneration_correction_functions.root", "RECREATE");
  std::ofstream resultFile("regeneration_fit_results.txt");


  for (auto &g : groups)
  {
    double rms = g.reg->GetRMS();

    // 1. Fit Danych
    auto resData = SmartFit(g.data, g.name + " (DATA)", g.pos, g.lowOff, g.highOff, g.lowSideband, g.highSideband, rms, 1, c1, true, false, false, true, false);
    TF1 *fData = (TF1 *)g.data->GetFunction("f_fit");

    // 2. Fit MC
    auto resMC = SmartFit(g.mc, g.name + " (MC SUM)", g.pos, g.lowOff, g.highOff, g.lowSideband, g.highSideband, rms, 2, c1, true, false, false, true, false);
    TF1 *fMC = (TF1 *)g.mc->GetFunction("f_fit");

    // 3. Wyznaczanie analitycznej funkcji poprawkowej
    c1->cd(3);
    gPad->Clear();

    // Zakres taki sam jak w histogramach powyżej (resMC[0] do resMC[1])
    Double_t xLow = resMC[0] - 1.0;
    Double_t xHigh = resMC[1] + 1.0;

    TF1 *f_weight = new TF1("f_weight", "([0]*exp(-0.5*((x-[1])/[2])^2)) / ([3]*exp(-0.5*((x-[4])/[5])^2))", xLow, xHigh);
    f_weight->SetParameters(fData->GetParameter(0), fData->GetParameter(1), fData->GetParameter(2),
                            fMC->GetParameter(0), fMC->GetParameter(1), fMC->GetParameter(2));
    f_weight->SetParError(0, fData->GetParError(0));
    f_weight->SetParError(1, fData->GetParError(1));
    f_weight->SetParError(2, fData->GetParError(2));
    f_weight->SetParError(3, fMC->GetParError(0));
    f_weight->SetParError(4, fMC->GetParError(1));
    f_weight->SetParError(5, fMC->GetParError(2));

    fMC->Write(g.name + "_Fit_MC");
    fData->Write(g.name + "_Fit_Data");
    f_weight->Write(g.name + "_WeightFunction");

    // Dynamiczne granice 3-sigma (do narysowania linii pionowych)
    Double_t leftMax = std::min(fData->GetParameter(1) - 2 * fData->GetParameter(2), fMC->GetParameter(1) - 2 * fMC->GetParameter(2));
    Double_t rightMax = std::max(fData->GetParameter(1) + 2 * fData->GetParameter(2), fMC->GetParameter(1) + 2 * fMC->GetParameter(2));

    resultFile << g.name << " Fit Results:" << std::endl;
    resultFile << "Weight Function Parameters: ";
    for (int i = 0; i < 6; i++)
    {
      resultFile << f_weight->GetParameter(i) << " +/- " << f_weight->GetParError(i) << "; ";
    }
    resultFile << std::endl;
    resultFile << "Range: " << leftMax << " to " << rightMax << std::endl;
    resultFile << "----------------------------------------" << std::endl;

    const int nSteps = 200;
    TGraphErrors *gr_band = new TGraphErrors();
    double step = (xHigh - xLow) / (nSteps - 1);

    double yMaxFound = 0;
    const double yAbsoluteLimit = 12.0;
    int nRealPoints = 0;

    for (int i = 0; i < nSteps; i++)
    {
      double x = xLow + i * step;
      double valData = fData->GetParameter(0) * exp(-0.5 * pow((x - fData->GetParameter(1)) / fData->GetParameter(2), 2));
      double valMC = fMC->GetParameter(0) * exp(-0.5 * pow((x - fMC->GetParameter(1)) / fMC->GetParameter(2), 2));

      // Zabezpieczenie przed dzieleniem przez zero i "wystrzałami" w nieskończoność
      if (valMC <= 1e-5)
        continue;

      double ratio = valData / valMC;
      if (ratio > yAbsoluteLimit * 2)
        continue; // Pomijaj punkty, które są daleko poza skalą

      double relErrData = fData->GetParError(0) / fData->GetParameter(0);
      double relErrMC = fMC->GetParError(0) / fMC->GetParameter(0);
      double totalErr = ratio * std::sqrt(relErrData * relErrData + relErrMC * relErrMC);

      gr_band->SetPoint(nRealPoints, x, ratio);
      gr_band->SetPointError(nRealPoints, 0, totalErr);
      nRealPoints++;

      if (ratio + totalErr > yMaxFound && ratio + totalErr < yAbsoluteLimit)
      {
        yMaxFound = ratio + totalErr;
      }
    }

    f_weight->SetTitle("Correction Factor W(" + g.name + ") with 1#sigma Band;" + g.name + ";W(x) = Data/MC");
    f_weight->SetLineColor(kRed);
    f_weight->SetLineWidth(2);
    f_weight->SetMinimum(0);

    if (yMaxFound > 0)
      f_weight->SetMaximum(yMaxFound * 1.3);
    else
      f_weight->SetMaximum(yAbsoluteLimit);

    // Najpierw rysujemy funkcję, która ustawia osie dokładnie w zakresie xLow - xHigh
    f_weight->Draw("L");

    if (nRealPoints > 0)
    {
      gr_band->SetFillColorAlpha(kRed - 7, 0.35);
      gr_band->SetMarkerSize(0);
      gr_band->Draw("E3 SAME");
    }

    // Linie pionowe wskazujące zakres 3-sigma (istotny fizycznie)
    TLine *l_low = new TLine(leftMax, 0, leftMax, f_weight->GetMaximum());
    TLine *l_high = new TLine(rightMax, 0, rightMax, f_weight->GetMaximum());
    l_low->SetLineStyle(2);
    l_low->SetLineColor(kBlack);
    l_low->Draw();
    l_high->SetLineStyle(2);
    l_high->SetLineColor(kBlack);
    l_high->Draw();

    TLegend *lw = new TLegend(0.15, 0.7, 0.45, 0.88);
    lw->SetBorderSize(0);
    lw->SetFillStyle(0);
    lw->AddEntry(f_weight, "Analytical Weight W(x)", "l");
    lw->AddEntry(gr_band, "1#sigma Uncertainty", "f");
    lw->AddEntry((TObject *)0, Form("Peak Ratio: %.3f", f_weight->Eval(fData->GetParameter(1))), "");
    lw->Draw();

    c1->Print(outName);

    g.weight_func = f_weight; // zachowaj do sekcji korekcji
    g.rangeLeft   = leftMax;
    g.rangeRight  = rightMax;
    delete gr_band;
  }

  c1->Print(outName + "]");
  delete c1;

  file->Close();
  resultFile.close();

  std::cout << "Gotowe! Wynik w: " << outName << std::endl;

  // --- Drugi przebieg: korekcja zdarzenie po zdarzeniu dla kanalu Regeneracji ---
  // Indeksy grup: 0=R_Charged, 1=Rho_Charged, 2=R_Neutral, 3=Rho_Neutral
  auto applyWeight = [&](int idx, double val) -> double
  {
    if (!groups[idx].weight_func) return 1.0;
    if (val >= groups[idx].rangeLeft && val <= groups[idx].rangeRight)
      return groups[idx].weight_func->Eval(val);
    return 1.0;
  };

  // Returns sigma_w(val) propagated from fit parameter errors (Data i MC fity niezalezne)
  auto weightErr = [&](int idx, double val) -> double
  {
    if (!groups[idx].weight_func) return 0.0;
    if (val < groups[idx].rangeLeft || val > groups[idx].rangeRight) return 0.0;
    TF1   *fw   = groups[idx].weight_func;
    double w    = fw->Eval(val);
    double A_D   = fw->GetParameter(0), mu_D   = fw->GetParameter(1), sig_D   = fw->GetParameter(2);
    double eA_D  = fw->GetParError(0),  emu_D  = fw->GetParError(1),  esig_D  = fw->GetParError(2);
    double A_MC  = fw->GetParameter(3), mu_MC  = fw->GetParameter(4), sig_MC  = fw->GetParameter(5);
    double eA_MC = fw->GetParError(3),  emu_MC = fw->GetParError(4),  esig_MC = fw->GetParError(5);
    double t_D   = (sig_D  > 0) ? (val - mu_D)  / sig_D  : 0.0;
    double t_MC  = (sig_MC > 0) ? (val - mu_MC) / sig_MC : 0.0;
    double relErr2 = std::pow(eA_D  / A_D,  2)
                   + std::pow(emu_D  * t_D  / sig_D,  2)
                   + std::pow(esig_D * t_D  * t_D / sig_D,  2)
                   + std::pow(eA_MC / A_MC, 2)
                   + std::pow(emu_MC * t_MC / sig_MC, 2)
                   + std::pow(esig_MC * t_MC * t_MC / sig_MC, 2);
    return w * std::sqrt(relErr2);
  };

  std::map<TString, TH1 *> h_rho_charged_corr, h_rho_neutral_corr, h_R_charged_corr, h_R_neutral_corr;
  for (const auto &chann : KLOE::channName)
  {
    h_rho_charged_corr[chann.second] = new TH1F("h_rho_ch_corr_" + chann.second, ";#rho_{charged} [cm];Entries", 251, 0, 50);
    h_rho_neutral_corr[chann.second] = new TH1F("h_rho_ne_corr_" + chann.second, ";#rho_{neutral} [cm];Entries", 251, 0, 50);
    h_R_charged_corr[chann.second]   = new TH1F("h_R_ch_corr_"   + chann.second, ";R_{charged} [cm];Entries",   251, 0, 50);
    h_R_neutral_corr[chann.second]   = new TH1F("h_R_ne_corr_"   + chann.second, ";R_{neutral} [cm];Entries",   251, 0, 50);
  }

  // Histogramy akumulujace sum(sigma_w^2) per bin dla kanacu Regeneracji
  TH1F *h_w2err_rho_ch = new TH1F("h_w2err_rho_ch", "", 251, 0, 50);
  TH1F *h_w2err_R_ch   = new TH1F("h_w2err_R_ch",   "", 251, 0, 50);
  TH1F *h_w2err_rho_ne = new TH1F("h_w2err_rho_ne", "", 251, 0, 50);
  TH1F *h_w2err_R_ne   = new TH1F("h_w2err_R_ne",   "", 251, 0, 50);

  reader.Restart();
  while (reader.Next())
  {
    std::array<Double_t, 3> distNeutralCharged = {KchrecFit[6] - KnerecFit[6],
                                                  KchrecFit[7] - KnerecFit[7],
                                                  KchrecFit[8] - KnerecFit[8]},
                            distNeutralIP = {Knerec[6] - *Bx,
                                             Knerec[7] - *By,
                                             Knerec[8] - KchrecClosest[8]},
                            distChargedIP = {KchrecClosest[6] - *Bx,
                                             KchrecClosest[7] - *By,
                                             KchrecClosest[8] - *Bz};

    Float_t
        rho_pm = std::sqrt(distChargedIP[0] * distChargedIP[0] + distChargedIP[1] * distChargedIP[1]),
        rho_00 = std::sqrt(distNeutralIP[0] * distNeutralIP[0] + distNeutralIP[1] * distNeutralIP[1]),
        rho    = std::sqrt(std::pow(rho_pm, 2) + std::pow(rho_00, 2));

    Double_t radius00 = std::sqrt(std::pow(Knerec[6] - *Bx, 2) + std::pow(Knerec[7] - *By, 2)),
             radiuspm = std::sqrt(std::pow(KchrecClosest[6] - *Bx, 2) + std::pow(KchrecClosest[7] - *By, 2)),
             zdist00  = std::abs(Knerec[8] - KchrecClosest[8]),
             zdistpm  = std::abs(KchrecClosest[8] - *Bz);

    Double_t fiducialVolume      = std::sqrt(std::pow(distNeutralCharged[0], 2) + std::pow(distNeutralCharged[1], 2)) < 2.05 && std::abs(distNeutralCharged[2]) < 2.45,
             fiducialVolumeClose = radius00 < 1.5 && radiuspm < 2.0 && zdist00 < 1.5 && zdistpm < 1.5;

    TVector3 KchrecVec  = {KchrecFit[0], KchrecFit[1], KchrecFit[2]};
    TVector3 trk1VecFit = {trk1Fit[0], trk1Fit[1], trk1Fit[2]};
    TVector3 trk2VecFit = {trk2Fit[0], trk2Fit[1], trk2Fit[2]};

    Double_t phiTrk1Angle = cos(trk1VecFit.Angle(KchrecVec)),
             phiTrk2Angle = cos(trk2VecFit.Angle(KchrecVec));

    Bool_t global_cut = ((fiducialVolume && (abs(phiTrk1Angle) < 0.8 || abs(phiTrk2Angle) < 0.8)) || !fiducialVolume) && ((fiducialVolumeClose && rho > 1.5) || !fiducialVolumeClose) && ellipse_cut(*minv4gam - PhysicsConstants::mK0, KnerecSix[5] - PhysicsConstants::mK0);

    if (!(*mcflag == 1 && global_cut)) continue;
    if (*mctruth == -1 || *mctruth == 0) continue;

    kch_position.SetXYZ(Kchrec[6], Kchrec[7], Kchrec[8]);
    kne_position.SetXYZ(Knerec[6], Knerec[7], Knerec[8]);
    ip_position.SetXYZ(0,0,0);//ip[0], ip[1], ip[2]);
    kch_position -= ip_position;
    kne_position -= ip_position;

    Double_t R_charged   = kch_position.Mag(),
             R_neutral   = kne_position.Mag(),
             rho_charged = kch_position.Perp(),
             rho_neutral = kne_position.Perp();

    TString ch   = KLOE::channName.at(*mctruth);
    bool isRegen = (ch == "Regeneration");

    if (*time_ne_CM_signal_fit < 7.0)
    {
      double w_ch = 1.0;
      if (isRegen)
      {
        double w1 = applyWeight(0, R_charged),   ew1 = weightErr(0, R_charged);
        double w2 = applyWeight(1, rho_charged),  ew2 = weightErr(1, rho_charged);
        w_ch = w1 * w2;
        double sigma2 = std::pow(w2 * ew1, 2) + std::pow(w1 * ew2, 2);
        h_w2err_rho_ch->AddBinContent(h_w2err_rho_ch->FindBin(rho_charged), sigma2);
        h_w2err_R_ch->AddBinContent(h_w2err_R_ch->FindBin(R_charged),       sigma2);
      }
      h_rho_charged_corr[ch]->Fill(rho_charged, w_ch);
      h_R_charged_corr[ch]->Fill(R_charged, w_ch);
    }

    if (*time_ch_CM_signal_fit < 7.0)
    {
      double w_ne = 1.0;
      if (isRegen)
      {
        double w1 = applyWeight(2, R_neutral),   ew1 = weightErr(2, R_neutral);
        double w2 = applyWeight(3, rho_neutral),  ew2 = weightErr(3, rho_neutral);
        w_ne = w1 * w2;
        double sigma2 = std::pow(w2 * ew1, 2) + std::pow(w1 * ew2, 2);
        h_w2err_rho_ne->AddBinContent(h_w2err_rho_ne->FindBin(rho_neutral), sigma2);
        h_w2err_R_ne->AddBinContent(h_w2err_R_ne->FindBin(R_neutral),       sigma2);
      }
      h_rho_neutral_corr[ch]->Fill(rho_neutral, w_ne);
      h_R_neutral_corr[ch]->Fill(R_neutral, w_ne);
    }
  }

  for (const auto &chann : KLOE::channName)
  {
    if (chann.second == "Data" || chann.second == "MC sum") continue;
    h_rho_charged_corr["MC sum"]->Add(h_rho_charged_corr[chann.second], scale_factor);
    h_rho_neutral_corr["MC sum"]->Add(h_rho_neutral_corr[chann.second], scale_factor);
    h_R_charged_corr["MC sum"]->Add(h_R_charged_corr[chann.second], scale_factor);
    h_R_neutral_corr["MC sum"]->Add(h_R_neutral_corr[chann.second], scale_factor);
  }

  // Propagacja niepewnosci funkcji wagowej do bledow MC sum (dodajemy w kwadraturze)
  double sf2 = scale_factor * scale_factor;
  auto addWeightErrToBin = [&](TH1 *h, TH1F *h_w2err)
  {
    for (int ib = 1; ib <= h->GetNbinsX(); ib++)
    {
      double e_old = h->GetBinError(ib);
      h->SetBinError(ib, std::sqrt(e_old * e_old + sf2 * h_w2err->GetBinContent(ib)));
    }
  };
  addWeightErrToBin(h_rho_charged_corr["MC sum"], h_w2err_rho_ch);
  addWeightErrToBin(h_R_charged_corr["MC sum"],   h_w2err_R_ch);
  addWeightErrToBin(h_rho_neutral_corr["MC sum"], h_w2err_rho_ne);
  addWeightErrToBin(h_R_neutral_corr["MC sum"],   h_w2err_R_ne);
  // ---

  TString max_chann_charged_rho = "", max_chann_neutral_rho = "", max_chann_charged_R = "", max_chann_neutral_R = "";
  double max_charged_rho = 0, max_neutral_rho = 0, max_charged_R = 0, max_neutral_R = 0;

  // --- Checking which histogram is highest ---
  for (const auto &chann : KLOE::channName)
  {
    double max_charged_rho_tmp = h_rho_charged[chann.second]->GetMaximum();
    double max_neutral_rho_tmp = h_rho_neutral[chann.second]->GetMaximum();
    double max_charged_R_tmp = h_R_charged[chann.second]->GetMaximum();
    double max_neutral_R_tmp = h_R_neutral[chann.second]->GetMaximum();

    if (max_charged_rho_tmp > max_charged_rho)
    {
      max_charged_rho = max_charged_rho_tmp;
      max_chann_charged_rho = chann.second;
    }
    if (max_neutral_rho_tmp > max_neutral_rho)
    {
      max_neutral_rho = max_neutral_rho_tmp;
      max_chann_neutral_rho = chann.second;
    }
    if (max_charged_R_tmp > max_charged_R)
    {
      max_charged_R = max_charged_R_tmp;
      max_chann_charged_R = chann.second;
    }
    if (max_neutral_R_tmp > max_neutral_R)
    {
      max_neutral_R = max_neutral_R_tmp;
      max_chann_neutral_R = chann.second;
    }
  }

  h_rho_charged["Regeneration"]->GetXaxis()->SetRangeUser(0, 50);
  h_R_charged["Regeneration"]->GetXaxis()->SetRangeUser(0, 50);

  h_rho_neutral["Regeneration"]->GetXaxis()->SetRangeUser(0, 50);
  h_R_neutral["Regeneration"]->GetXaxis()->SetRangeUser(0, 50);

  h_rho_charged["MC sum"]->GetXaxis()->SetRangeUser(0, 50);
  h_R_charged["MC sum"]->GetXaxis()->SetRangeUser(0, 50);

  h_rho_neutral["MC sum"]->GetXaxis()->SetRangeUser(0, 50);
  h_R_neutral["MC sum"]->GetXaxis()->SetRangeUser(0, 50);

  h_rho_charged["Data"]->GetXaxis()->SetRangeUser(0, 50);
  h_R_charged["Data"]->GetXaxis()->SetRangeUser(0, 50);

  h_rho_neutral["Data"]->GetXaxis()->SetRangeUser(0, 50);
  h_R_neutral["Data"]->GetXaxis()->SetRangeUser(0, 50);

  TCanvas *c2 = new TCanvas("c_charged", "Charged", 800, 600);

  h_rho_charged[max_chann_charged_rho]->SetTitle("Rho Charged (Full Scale);#rho_{charged};Entries");
  h_rho_charged[max_chann_charged_rho]->GetYaxis()->SetRangeUser(0.0, 1.2 * max_charged_rho);

  if (max_chann_charged_rho != "Data")
    h_rho_charged[max_chann_charged_rho]->Draw("HIST");
  else
    h_rho_charged[max_chann_charged_rho]->Draw("PE1");

  // --- Draw both distributions full scale ---
  TLegend *l = new TLegend(0.6, 0.7, 0.85, 0.88);
  for (const auto &chann : KLOE::channName)
  {
    h_rho_charged[chann.second]->SetLineColor(KLOE::channColor.at(chann.second));

    if (chann.second == max_chann_charged_rho)
      continue;

    if (chann.second != "Data")
      h_rho_charged[chann.second]->Draw("HISTSAME");
    else
      h_rho_charged[chann.second]->Draw("PE1SAME");

    if (chann.second != "Data")
      l->AddEntry(h_rho_charged[chann.second], KLOE::channTitle.at(chann.second), "l");
    else
      l->AddEntry(h_rho_charged[chann.second], KLOE::channTitle.at(chann.second), "lep1");
  }
  l->Draw();

  c2->Print("full_scale_charged_rho.pdf");
  c2->Clear();
  l->Clear();

  // --- Draw both distributions full scale ---

  h_R_charged[max_chann_charged_R]->SetTitle("R Charged (Full Scale);R_{charged};Entries");
  h_R_charged[max_chann_charged_R]->GetYaxis()->SetRangeUser(0.0, 1.2 * max_charged_R);

  if (max_chann_charged_R != "Data")
    h_R_charged[max_chann_charged_R]->Draw("HIST");
  else
    h_R_charged[max_chann_charged_R]->Draw("PE1");

  for (const auto &chann : KLOE::channName)
  {
    h_R_charged[chann.second]->SetLineColor(KLOE::channColor.at(chann.second));

    if (chann.second == max_chann_charged_R)
      continue;

    if (chann.second != "Data")
      h_R_charged[chann.second]->Draw("HISTSAME");
    else
      h_R_charged[chann.second]->Draw("PE1SAME");

    if (chann.second != "Data")
      l->AddEntry(h_R_charged[chann.second], KLOE::channTitle.at(chann.second), "l");
    else
      l->AddEntry(h_R_charged[chann.second], KLOE::channTitle.at(chann.second), "lep1");
  }
  l->Draw();

  c2->Print("full_scale_charged_R.pdf");
  c2->Clear();
  l->Clear();

  // --- Repeat for neutral ---

  h_rho_neutral[max_chann_neutral_rho]->SetTitle("Rho Neutral (Full Scale);#rho_{neutral};Entries");
  h_rho_neutral[max_chann_neutral_rho]->GetYaxis()->SetRangeUser(0.0, 1.2 * max_neutral_rho);

  if (max_chann_neutral_rho != "Data")
    h_rho_neutral[max_chann_neutral_rho]->Draw("HIST");
  else
    h_rho_neutral[max_chann_neutral_rho]->Draw("PE1");

  // --- Draw both distributions full scale ---
  for (const auto &chann : KLOE::channName)
  {
    h_rho_neutral[chann.second]->SetLineColor(KLOE::channColor.at(chann.second));

    if (chann.second == max_chann_neutral_rho)
      continue;

    if (chann.second != "Data")
      h_rho_neutral[chann.second]->Draw("HISTSAME");
    else
      h_rho_neutral[chann.second]->Draw("PE1SAME");

    if (chann.second != "Data")
      l->AddEntry(h_rho_neutral[chann.second], KLOE::channTitle.at(chann.second), "l");
    else
      l->AddEntry(h_rho_neutral[chann.second], KLOE::channTitle.at(chann.second), "lep1");
  }
  l->Draw();

  c2->Print("full_scale_neutral_rho.pdf");
  c2->Clear();
  l->Clear();

  // --- Draw both distributions full scale ---
  h_R_neutral[max_chann_neutral_R]->SetTitle("R Neutral (Full Scale);R_{neutral};Entries");
  h_R_neutral[max_chann_neutral_R]->GetYaxis()->SetRangeUser(0.0, 1.2 * max_neutral_R);

  if (max_chann_neutral_R != "Data")
    h_R_neutral[max_chann_neutral_R]->Draw("HIST");
  else
    h_R_neutral[max_chann_neutral_R]->Draw("PE1");

  for (const auto &chann : KLOE::channName)
  {
    h_R_neutral[chann.second]->SetLineColor(KLOE::channColor.at(chann.second));

    if (chann.second == max_chann_neutral_R)
      continue;

    if (chann.second != "Data")
      h_R_neutral[chann.second]->Draw("HISTSAME");
    else
      h_R_neutral[chann.second]->Draw("PE1SAME");

    if (chann.second != "Data")
      l->AddEntry(h_R_neutral[chann.second], KLOE::channTitle.at(chann.second), "l");
    else
      l->AddEntry(h_R_neutral[chann.second], KLOE::channTitle.at(chann.second), "lep1");
  }
  l->Draw();

  c2->Print("full_scale_neutral_R.pdf");
  c2->Clear();
  l->Clear();

  for (const auto &chann : KLOE::channName)
  {
    h_rho_neutral_vs_R[chann.second]->Draw("COLZ");

    c2->Print("rho_neutral_vs_R_" + chann.second + ".pdf");
    c2->Clear();
  }

  for (const auto &chann : KLOE::channName)
  {
    h_rho_charged_vs_R[chann.second]->Draw("COLZ");

    c2->Print("rho_charged_vs_R_" + chann.second + ".pdf");
    c2->Clear();
  }

  // --- SEKCJA RYSOWANIA: Data vs MC sum po korekcji event-by-event + residua ---
  struct CorrPlotGroup
  {
    TString name;
    TString fileName;
    TString xLabel;
    TH1 *data;
    TH1 *mc_orig;
    TH1 *mc_corr;
  };

  std::vector<CorrPlotGroup> corrPlots = {
    {"R_{charged}",    "R_charged",   "R_{charged} [cm]",    h_R_charged["Data"],   h_R_charged["MC sum"],   h_R_charged_corr["MC sum"]},
    {"#rho_{charged}", "Rho_charged", "#rho_{charged} [cm]", h_rho_charged["Data"], h_rho_charged["MC sum"], h_rho_charged_corr["MC sum"]},
    {"R_{neutral}",    "R_neutral",   "R_{neutral} [cm]",    h_R_neutral["Data"],   h_R_neutral["MC sum"],   h_R_neutral_corr["MC sum"]},
    {"#rho_{neutral}", "Rho_neutral", "#rho_{neutral} [cm]", h_rho_neutral["Data"], h_rho_neutral["MC sum"], h_rho_neutral_corr["MC sum"]},
  };

  for (auto &g : groups)
  {
    delete g.weight_func;
    g.weight_func = nullptr;
  }

  TCanvas *c_corr = new TCanvas("c_corr", "Corrected Comparison", 800, 800);

  for (auto &cg : corrPlots)
  {
    TH1F *h_resid = (TH1F *)cg.data->Clone("h_resid_" + cg.fileName);
    h_resid->SetDirectory(nullptr);
    h_resid->Reset();
    for (int ib = 1; ib <= cg.data->GetNbinsX(); ib++)
    {
      double data_val = cg.data->GetBinContent(ib);
      double data_err = sqrt(pow(cg.data->GetBinError(ib), 2) + pow(cg.mc_corr->GetBinError(ib), 2));
      double mc_val   = cg.mc_corr->GetBinContent(ib);
      if (data_err > 0)
      {
        h_resid->SetBinContent(ib, (data_val - mc_val) / data_err);
        h_resid->SetBinError(ib, 1.0);
      }
    }

    c_corr->Clear();

    TPad *pad_top = new TPad("top_" + cg.fileName, "", 0, 0.30, 1, 1);
    pad_top->SetBottomMargin(0.015);
    pad_top->SetLeftMargin(0.12);
    pad_top->Draw();

    TPad *pad_bot = new TPad("bot_" + cg.fileName, "", 0, 0, 1, 0.30);
    pad_bot->SetTopMargin(0.015);
    pad_bot->SetBottomMargin(0.35);
    pad_bot->SetLeftMargin(0.12);
    pad_bot->Draw();

    pad_top->cd();
    double yMax = std::max({cg.data->GetMaximum(), cg.mc_orig->GetMaximum(), cg.mc_corr->GetMaximum()});
    cg.data->GetXaxis()->SetRangeUser(0, 50);
    cg.data->GetXaxis()->SetLabelSize(0);
    cg.data->GetYaxis()->SetRangeUser(0, 1.35 * yMax);
    cg.data->SetTitle("MC_{corr} vs Data: " + cg.name + ";;Entries");
    ((TH1F *)cg.data)->SetMarkerStyle(20);
    ((TH1F *)cg.data)->SetMarkerSize(0.7);
    cg.data->Draw("PE1");

    cg.mc_orig->SetLineColor(kGray + 1);
    cg.mc_orig->SetLineStyle(2);
    cg.mc_orig->GetXaxis()->SetRangeUser(0, 50);
    cg.mc_orig->Draw("HIST SAME");

    cg.mc_corr->SetLineColor(kRed + 1);
    cg.mc_corr->SetLineWidth(2);
    cg.mc_corr->GetXaxis()->SetRangeUser(0, 50);
    cg.mc_corr->Draw("HIST SAME");

    TLegend *l_corr = new TLegend(0.65, 0.62, 0.90, 0.88);
    l_corr->SetBorderSize(0);
    l_corr->SetFillStyle(0);
    l_corr->SetTextSize(0.040);
    l_corr->AddEntry(cg.data,    "Data",           "lep");
    l_corr->AddEntry(cg.mc_orig, "MC sum (oryg.)", "l");
    l_corr->AddEntry(cg.mc_corr, "MC sum (skor.)", "l");
    l_corr->Draw();

    pad_bot->cd();
    h_resid->SetTitle(";" + cg.xLabel + ";(Data-MC_{corr})/#sigma");
    h_resid->GetXaxis()->SetRangeUser(0, 50);
    h_resid->GetXaxis()->SetTitleSize(0.13);
    h_resid->GetXaxis()->SetLabelSize(0.10);
    h_resid->GetYaxis()->SetTitleSize(0.11);
    h_resid->GetYaxis()->SetTitleOffset(0.55);
    h_resid->GetYaxis()->SetLabelSize(0.09);
    h_resid->GetYaxis()->SetNdivisions(505);
    h_resid->GetYaxis()->SetRangeUser(-5, 5);
    h_resid->SetLineColor(kBlack);
    h_resid->SetMarkerColor(kBlack);
    h_resid->SetMarkerStyle(20);
    h_resid->Draw("PE1");

    TLine *zero_line = new TLine(0, 0, 50, 0);
    zero_line->SetLineStyle(2);
    zero_line->SetLineColor(kRed);
    zero_line->Draw();

    c_corr->Print("full_scale_corrected_" + cg.fileName + ".pdf");

    delete h_resid;
  }

  delete c_corr;
}