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

std::vector<double> SmartFit(TH1F *h, TString title, double peakPos, double lowOff, double highOff, double rms, int padIdx, TCanvas *c, bool useExpo = false, bool useLangau = false, bool useAsym = false, TFitResultPtr &resOut = *(TFitResultPtr *)nullptr)
{
  c->cd(padIdx);
  h->SetTitle(title);

  double xMin = peakPos - lowOff;
  double xMax = peakPos + highOff;

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
  else
  {
    sigFormula = "[0]*exp(-0.5*((x-[1])/[2])^2)";
  }

  TString bgFormula = useExpo ? "expo(3)" : "pol1(3)";
  TF1 *f_fit = new TF1("f_fit", sigFormula + " + " + bgFormula, xMin, xMax);

  if (useExpo)
  {
    f_fit->SetParameters(max_val * 0.8, peakPos, 0.1 * rms, TMath::Log(max_val * 0.1 + 1e-9), -0.05);
    f_fit->SetParLimits(4, -10.0, 0.0); // Stabilizacja: tło wykładnicze nie może rosnąć
  }
  else
  {
    f_fit->SetParameters(max_val * 0.8, peakPos, 0.1 * rms, max_val * 0.1, 0.0);
    f_fit->SetParLimits(3, 0.0, max_val); // Stabilizacja: tło nie może być ujemne
    f_fit->SetParLimits(4, -10.0, 0.0); // Stabilizacja: współczynnik liniowy tła w rozsądnym zakresie
  }

  if (useAsym)
  {
    f_fit->SetParameter(5, 0.1 * rms);
    f_fit->SetParLimits(5, 0.01 * rms, 2.0 * rms);
  }

  f_fit->SetParLimits(0, 0.01 * max_val, 10 * max_val);
  f_fit->SetParLimits(1, peakPos - 0.5, peakPos + 0.5); // Zacieśnienie limitu pozycji piku
  f_fit->SetParLimits(2, 0.01 * rms, 1.5 * rms);        // Zacieśnienie limitu szerokości

  // ZMIANA 1: Opcja "I" (Integral) oraz "M" (Improve)
  // "I" używa całki z funkcji w binie (kluczowe przy wąskich pikach i małym binowaniu)
  // "M" szuka lepszych minimów, co stabilizuje macierz kowariancji
  TFitResultPtr r = h->Fit(f_fit, "RMSLQI");
  resOut = r;

  h->Draw("E1");

  f_fit->GetXaxis()->SetRangeUser(drawMin, drawMax);

  TF1 *sigOnly = new TF1("sigOnly", sigFormula, xMin, xMax);
  sigOnly->SetParameter(0, f_fit->GetParameter(0));
  sigOnly->SetParameter(1, f_fit->GetParameter(1));
  sigOnly->SetParameter(2, f_fit->GetParameter(2));
  if (useAsym)
    sigOnly->SetParameter(5, f_fit->GetParameter(5));
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
    sigFuncForErr->FixParameter(3, 0);
    sigFuncForErr->FixParameter(4, 0);

    // ZMIANA 3: Obliczanie błędu z uwzględnieniem pełnej macierzy kowariancji fita
    YieldErr = sigFuncForErr->IntegralError(low, high,
                                            r->GetParams(),
                                            r->GetCovarianceMatrix().GetMatrixArray()) /
               binW;

    delete sigFuncForErr;
  }

  TF1 *bg = new TF1("bg", useExpo ? "expo" : "pol1", xMin, xMax);
  bg->SetParameters(f_fit->GetParameter(3), f_fit->GetParameter(4));
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
  };
  std::vector<PlotGroup> groups = {
      {"R_Charged", (TH1F *)h_R_charged["Data"], (TH1F *)h_R_charged["MC sum"], (TH1F *)h_R_charged["Regeneration"], 10.0, 4.0, 6.0},
      {"Rho_Charged", (TH1F *)h_rho_charged["Data"], (TH1F *)h_rho_charged["MC sum"], (TH1F *)h_rho_charged["Regeneration"], 25.0, 6.0, 6.0},
      {"R_Neutral", (TH1F *)h_R_neutral["Data"], (TH1F *)h_R_neutral["MC sum"], (TH1F *)h_R_neutral["Regeneration"], 10.0, 4.0, 6.0},
      {"Rho_Neutral", (TH1F *)h_rho_neutral["Data"], (TH1F *)h_rho_neutral["MC sum"], (TH1F *)h_rho_neutral["Regeneration"], 25.0, 6.0, 6.0}};

  TFile *file = new TFile("regeneration_correction_functions.root", "RECREATE");
  std::ofstream resultFile("regeneration_fit_results.txt");

  for (auto &g : groups)
  {
    double rms = g.reg->GetRMS();

    // 1. Fit Danych
    TFitResultPtr fitResultData;
    auto resData = SmartFit(g.data, g.name + " (DATA)", g.pos, g.lowOff, g.highOff, rms, 1, c1, false, false, false, fitResultData);
    TF1 *fData = (TF1 *)g.data->GetFunction("f_fit");

    // 2. Fit MC
    TFitResultPtr fitResultMC;
    auto resMC = SmartFit(g.mc, g.name + " (MC SUM)", g.pos, g.lowOff, g.highOff, rms, 2, c1, false, false, false, fitResultMC);
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
    Double_t leftMax = std::min(fData->GetParameter(1) - 3 * fData->GetParameter(2), fMC->GetParameter(1) - 3 * fMC->GetParameter(2));
    Double_t rightMax = std::max(fData->GetParameter(1) + 3 * fData->GetParameter(2), fMC->GetParameter(1) + 3 * fMC->GetParameter(2));

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

    // --- ZMIENNE NA GRANICE BEZPIECZNEGO PRZEDZIAŁU (Wykonaj przed pętlą) ---
    double safeXMin = -999.0;
    double safeXMax = -999.0;
    bool foundFirstSafe = false;

    // --- WYCIĄGANIE KOWARIANCJI Z WYNIKÓW FITU (Wykonaj koniecznie PRZED pętlą) ---
    double cov01_Data = 0, cov02_Data = 0, cov12_Data = 0;
    double cov01_MC = 0, cov02_MC = 0, cov12_MC = 0;

    if (fitResultData.Get() != nullptr && fitResultData->IsValid())
    {
      cov01_Data = fitResultData->CovMatrix(0, 1); // cov(A, mu)
      cov02_Data = fitResultData->CovMatrix(0, 2); // cov(A, sigma)
      cov12_Data = fitResultData->CovMatrix(1, 2); // cov(mu, sigma)
    }

    if (fitResultMC.Get() != nullptr && fitResultMC->IsValid())
    {
      cov01_MC = fitResultMC->CovMatrix(0, 1); // cov(A, mu)
      cov02_MC = fitResultMC->CovMatrix(0, 2); // cov(A, sigma)
      cov12_MC = fitResultMC->CovMatrix(1, 2); // cov(mu, sigma)
    }

    for (int i = 0; i < nSteps; i++)
    {
      double x = xLow + i * step;

      // --- DANE DLA DATA ---
      double A_Data = fData->GetParameter(0);
      double mu_Data = fData->GetParameter(1);
      double sigma_Data = fData->GetParameter(2);

      double uA_Data = fData->GetParError(0);
      double uMu_Data = fData->GetParError(1);
      double uSigma_Data = fData->GetParError(2);

      double valData = A_Data * std::exp(-0.5 * std::pow((x - mu_Data) / sigma_Data, 2));

      // --- DANE DLA MC ---
      double A_MC = fMC->GetParameter(0);
      double mu_MC = fMC->GetParameter(1);
      double sigma_MC = fMC->GetParameter(2);

      double uA_MC = fMC->GetParError(0);
      double uMu_MC = fMC->GetParError(1);
      double uSigma_MC = fMC->GetParError(2);

      double valMC = A_MC * std::exp(-0.5 * std::pow((x - mu_MC) / sigma_MC, 2));

      // Zabezpieczenie przed dzieleniem przez zero
      if (valMC <= 1e-5)
        continue;

      double ratio = valData / valMC;
      if (ratio > yAbsoluteLimit * 2)
        continue; // Pomijaj punkty daleko poza skalą

      // --- PROPAGACJA BŁĘDU DLA DATA (Z PEŁNĄ KOWARIANCJĄ) ---
      // Pochodne cząstkowe (jakobian)
      double df_dA_Data = std::exp(-0.5 * std::pow((x - mu_Data) / sigma_Data, 2));
      double df_dmu_Data = valData * (x - mu_Data) / (sigma_Data * sigma_Data);
      double df_dsigma_Data = valData * std::pow(x - mu_Data, 2) / std::pow(sigma_Data, 3);

      // Pełna wariancja (człony kwadratowe + człony kowariancji)
      double varData = std::pow(df_dA_Data * uA_Data, 2) +
                       std::pow(df_dmu_Data * uMu_Data, 2) +
                       std::pow(df_dsigma_Data * uSigma_Data, 2) +
                       2.0 * df_dA_Data * df_dmu_Data * cov01_Data +
                       2.0 * df_dA_Data * df_dsigma_Data * cov02_Data +
                       2.0 * df_dmu_Data * df_dsigma_Data * cov12_Data;

      // Wyznaczenie błędu względnego z zabezpieczeniem przed ujemną wariancją numeryczną
      double relErrData = (valData > 1e-5 && varData > 0) ? std::sqrt(varData) / valData : 0;

      // --- PROPAGACJA BŁĘDU DLA MC (Z PEŁNĄ KOWARIANCJĄ) ---
      // Pochodne cząstkowe (jakobian)
      double df_dA_MC = std::exp(-0.5 * std::pow((x - mu_MC) / sigma_MC, 2));
      double df_dmu_MC = valMC * (x - mu_MC) / (sigma_MC * sigma_MC);
      double df_dsigma_MC = valMC * std::pow(x - mu_MC, 2) / std::pow(sigma_MC, 3);

      // Pełna wariancja
      double varMC = std::pow(df_dA_MC * uA_MC, 2) +
                     std::pow(df_dmu_MC * uMu_MC, 2) +
                     std::pow(df_dsigma_MC * uSigma_MC, 2) +
                     2.0 * df_dA_MC * df_dmu_MC * cov01_MC +
                     2.0 * df_dA_MC * df_dsigma_MC * cov02_MC +
                     2.0 * df_dmu_MC * df_dsigma_MC * cov12_MC;

      double relErrMC = (valMC > 1e-5 && varMC > 0) ? std::sqrt(varMC) / valMC : 0;

      // --- CAŁKOWITY BŁĄD DLA RATIO ---
      double relErrRatio = std::sqrt(relErrData * relErrData + relErrMC * relErrMC);
      double totalErr = ratio * relErrRatio; // Błąd bezwzględny dla ratio

      // Zapisujemy granice przedziału, w którym błąd jest bezpieczny
      if (!foundFirstSafe && relErrRatio <= 1.0)
      {
        safeXMin = x; // Pierwszy punkt pętli spełniający kryterium < 50%
        foundFirstSafe = true;
      }

      if (relErrRatio <= 1.00)
      {
        safeXMax = x; // Nadpisuje się aż do ostatniego poprawnego punktu
      }

      // --- ZAPIS DO GRAPHU ---
      gr_band->SetPoint(nRealPoints, x, ratio);
      gr_band->SetPointError(nRealPoints, 0, totalErr);
      nRealPoints++;

      if (ratio + totalErr > yMaxFound && ratio + totalErr < yAbsoluteLimit)
      {
        yMaxFound = ratio + totalErr;
      }
    }

    // --- WYPISANIE WYNIKÓW W KONSOLI ---
    std::cout << "\n=======================================================" << std::endl;
    std::cout << " ZAKRES STOSOWALNOŚCI POPRAWKI (Błąd względny < 50%):" << std::endl;
    std::cout << " NAZWA GRUPY: " << g.name << std::endl;
    if (foundFirstSafe)
    {
      std::cout << " X_min = " << safeXMin << std::endl;
      std::cout << " X_max = " << safeXMax << std::endl;
      std::cout << " Szerokość okna bezpiecznego: " << (safeXMax - safeXMin) << std::endl;
    }
    else
    {
      std::cout << " ALERT: Ani jeden punkt nie spełnił kryterium błędu < 50%!" << std::endl;
    }
    std::cout << "=======================================================\n"
              << std::endl;

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
    TLine *l_low = new TLine(safeXMin, 0, safeXMin, f_weight->GetMaximum());
    TLine *l_high = new TLine(safeXMax, 0, safeXMax, f_weight->GetMaximum());
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

    delete f_weight;
    delete gr_band;
  }

  c1->Print(outName + "]");
  delete c1;

  file->Close();
  resultFile.close();

  std::cout << "Gotowe! Wynik w: " << outName << std::endl;

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
}