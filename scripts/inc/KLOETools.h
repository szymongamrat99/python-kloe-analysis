#ifndef KLOETOOLS_H
#define KLOETOOLS_H

#include "TChain.h"
#include "TString.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TVector3.h"
#include <array>
#include <iostream>
#include <functional>
#include <cmath>

#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TF1.h>
#include <TF2.h>
#include <TEllipse.h>
#include <TMath.h>
#include <TFitResult.h>

#include "const.h"

struct EventKinematics
{
  double R_charged = 0.0;
  double R_neutral = 0.0;
  double rho_charged = 0.0;
  double rho_neutral = 0.0;
};

namespace KloeTools
{

  // Funkcja ładująca pliki MC i Data do podanego TChain
  void LoadSignalChains(TChain *chain, const TString &root_file_date)
  {
    // Sprawdzenie na wypadek, gdyby przekazano pusty wskaźnik
    if (!chain)
    {
      std::cerr << "[-] BŁĄD: Przekazano pusty wskaźnik TChain do LoadSignalChains!" << std::endl;
      return;
    }

    TString base_path = "/data/ssd/gamrat/CNAF_Produced_Files/root_files/" + root_file_date + "/Signal/";
    std::array<TString, 3> mc_paths = {
        "ALL_PHYS_SIGNAL_NoSmearing/",
        "ALL_PHYS2_SIGNAL_NoSmearing/",
        "ALL_PHYS3_SIGNAL_NoSmearing/"};
    TString data_path = base_path + "DATA_SIGNAL_NoSmearing/*.root";

    // Ładowanie plików Monte Carlo
    for (const auto &path : mc_paths)
    {
      TString full_mc_path = base_path + path + "*.root";
      chain->Add(full_mc_path);
      std::cout << "[+] Dodano do chain: " << full_mc_path << std::endl;
    }

    // Ładowanie plików danych eksperymentalnych
    chain->Add(data_path);
    std::cout << "[+] Dodano do chain: " << data_path << std::endl;
  }

  bool PassesGlobalCuts(
      const TTreeReaderArray<Double_t> &KchrecFit, const TTreeReaderArray<Double_t> &KnerecFit,
      const TTreeReaderArray<Double_t> &Knerec, const TTreeReaderArray<Double_t> &KchrecClosest,
      const TTreeReaderArray<Double_t> &trk1Fit, const TTreeReaderArray<Double_t> &trk2Fit,
      double Bx, double By, double Bz, double minv4gam, double KnerecSix_5,
      std::function<bool(double, double)> ellipse_cut)
  {
    // 1. Definicja odległości
    double distNeutralCharged[3] = {KchrecFit[6] - KnerecFit[6], KchrecFit[7] - KnerecFit[7], KchrecFit[8] - KnerecFit[8]};
    double distNeutralIP[3] = {Knerec[6] - Bx, Knerec[7] - By, Knerec[8] - KchrecClosest[8]};
    double distChargedIP[3] = {KchrecClosest[6] - Bx, KchrecClosest[7] - By, KchrecClosest[8] - Bz};

    // 2. Promienie cylindryczne i odległości Z
    double rho_pm2 = distChargedIP[0] * distChargedIP[0] + distChargedIP[1] * distChargedIP[1];
    double rho_002 = distNeutralIP[0] * distNeutralIP[0] + distNeutralIP[1] * distNeutralIP[1];
    double rho_total = std::sqrt(rho_pm2 + rho_002);

    double radius00 = std::sqrt(distNeutralIP[0] * distNeutralIP[0] + distNeutralIP[1] * distNeutralIP[1]);
    double radiuspm = std::sqrt(distChargedIP[0] * distChargedIP[0] + distChargedIP[1] * distChargedIP[1]);
    double zdist00 = std::abs(Knerec[8] - KchrecClosest[8]);
    double zdistpm = std::abs(KchrecClosest[8] - Bz);

    // 3. Sprawdzanie Fiducial Volume
    bool fiducialVolume = (distNeutralCharged[0] * distNeutralCharged[0] + distNeutralCharged[1] * distNeutralCharged[1] < 4.2025) && (std::abs(distNeutralCharged[2]) < 2.45);
    bool fiducialVolumeClose = radius00 < 1.5 && radiuspm < 2.0 && zdist00 < 1.5 && zdistpm < 1.5;

    // 4. Szybki test elipsy i bliskiego FV
    if (fiducialVolumeClose && rho_total <= 1.5)
      return false;
    if (!ellipse_cut(minv4gam - 497.611, KnerecSix_5 - 497.611))
      return false; // Wykorzystanie stałej masy K0

    // 5. Test kątowy (tylko jeśli wpadliśmy w fiducialVolume)
    if (fiducialVolume)
    {
      TVector3 KchrecVec(KchrecFit[0], KchrecFit[1], KchrecFit[2]);
      TVector3 trk1VecFit(trk1Fit[0], trk1Fit[1], trk1Fit[2]);
      TVector3 trk2VecFit(trk2Fit[0], trk2Fit[1], trk2Fit[2]);

      if (std::abs(std::cos(trk1VecFit.Angle(KchrecVec))) >= 0.8 &&
          std::abs(std::cos(trk2VecFit.Angle(KchrecVec))) >= 0.8)
        return false;
    }

    return true; // Zdarzenie przeszło pomyślnie!
  }

  double get_local_poly_max(TH1 *h, double x_min, double x_max)
  {
    int max_bin = h->GetMaximumBin();
    double max_x = h->GetBinCenter(max_bin);
    double bin_width = h->GetBinWidth(max_bin);

    // Definiujemy bardzo wąski zakres: tylko +- 3 biny wokół czubka
    TF1 *f_clover = new TF1("f_clover", "pol2", x_min, x_max);

    // Fitujemy bez rysowania ("Q" - quiet, "R" - restricted range)
    h->Fit(f_clover, "QR");

    double a = f_clover->GetParameter(2);
    double b = f_clover->GetParameter(1);

    // Wierzchołek paraboli: -b / 2a
    return -b / (2.0 * a);
  }

  double get_local_gauss_max(TH1 *h, double x_min, double x_max)
  {
    int max_bin = h->GetMaximumBin();
    double max_x = h->GetBinCenter(max_bin);
    double bin_width = h->GetBinWidth(max_bin);

    // Definiujemy bardzo wąski zakres: tylko +- 3 biny wokół czubka
    TF1 *f_gauss = new TF1("f_gauss", "gaus", x_min, x_max);

    // Fitujemy bez rysowania ("Q" - quiet, "R" - restricted range)
    h->Fit(f_gauss, "QR");

    return f_gauss->GetParameter(1); // Środek Gaussa
  }

  double f_circle_dist(double *x, double *par)
  {
    double dx = x[0] - par[0];
    double dy = x[1] - par[1];
    double r = std::sqrt(dx * dx + dy * dy);

    // Zwracamy odległość od brzegu okręgu.
    return r - par[2];
  }

  // Zmieniamy typ zwracany na double, ale wewnątrz wyciągamy też y0, bo to ono nas interesuje!
  TFitResultPtr fit_method_circle(TH2 *h, double x_min, double x_max)
  {
    TF2 *f_circle = new TF2("f_circle", f_circle_dist, x_min, x_max, x_min, x_max, 3);

    // Nazywamy parametry, żeby Minuit ładnie pisał w logach
    f_circle->SetParName(0, "x0");
    f_circle->SetParName(1, "y0");
    f_circle->SetParName(2, "R");

    // Inicjalizacja parametrów na wartości BLISKIE RZECZYWISTOŚCI w KLOE:
    f_circle->SetParameter(0, 0.0);  // Spodziewany środek X ok. 0
    f_circle->SetParameter(1, -0.5); // Spodziewany środek Y ok. 1.5 cm (Twój offset)
    f_circle->SetParameter(2, 25.0); // Rura wiązki KLOE ma promień ok. 10 cm (lub zmień na ~4.5 dla KLOE-2 IT)

    // Bardzo ważne: Nakładamy ograniczenia, żeby fit nie uciekł w kosmos
    f_circle->SetParLimits(0, -5.0, 5.0);
    f_circle->SetParLimits(1, -5.0, 5.0);
    f_circle->SetParLimits(2, 0, 27); // Promień musi być fizyczny

    // FITOWANIE:
    // "Q" - quiet (wyłącz potok tekstu)
    // "R" - użyj zakresu zdefiniowanego w TF2 (czyli x_min, x_max)
    // "W" - ustawia równe wagi dla binów, dzięki czemu algorytm szuka geometrycznego kształtu chmury punktów, a nie dopasowuje wysokości słupków!
    TFitResultPtr fit_result = h->Fit(f_circle, "QRNLWS");

    if (!fit_result->IsValid())
    {
      std::cerr << "[-] BŁĄD: Dopasowanie okręgu nie powiodło się! Sprawdź histogram i zakresy." << std::endl;
      return fit_result; // Zwracamy wynik, nawet jeśli jest niepoprawny, żeby można było debugować
    }

    double x0 = f_circle->GetParameter(0);
    double y0 = f_circle->GetParameter(1);
    double r = f_circle->GetParameter(2);

// 4. SAMODZIELNIE nakładamy dokładnie jeden okrąg o dopasowanych wymiarach
#include "TEllipse.h"
    TEllipse *visual_circle = new TEllipse(x0, y0, r, r);
    visual_circle->SetFillStyle(0);    // Środek przezroczysty, żeby nie zasłonić danych
    visual_circle->SetLineColor(kRed); // Wyrazisty czerwony kolor
    visual_circle->SetLineWidth(3);    // Grubsza linia dla widoczności

    h->GetListOfFunctions()->Add(visual_circle); // Dodajemy okrąg do wykresu
    h->GetListOfFunctions()->Remove(f_circle);   // Usuwamy funkcję fitującą, bo już mamy wizualizację

    // Zwracamy wynik dopasowania, aby można było go wykorzystać do dalszej analizy, np. porównania z metodą profilu polarnego
    return fit_result;
  }

  inline double f_polar_circle(double *x, double *par)
  {
    double phi = x[0];
    return par[0] + par[1] * std::cos(phi) + par[2] * std::sin(phi);
  }

  inline TFitResultPtr fit_method_profile(TProfile *prof, double expected_r)
  {
    if (!prof || prof->GetEntries() == 0)
    {
      std::cerr << "[-] BŁĄD: Profil jest pusty lub nie istnieje w fit_method_profile!" << std::endl;
      return TFitResultPtr();
    }

    // Czyszczenie binów o niskiej statystyce (usuwanie rozproszonego tła), aby ustabilizować wyznaczanie R
    for (int b = 1; b <= prof->GetNbinsX(); ++b)
    {
      if (prof->GetBinEntries(b) < 5)
      { // Próg statystyczny: minimum 5 zliczeń w pasku kątowym
        prof->SetBinContent(b, 0);
        prof->SetBinError(b, 0);
      }
    }

    // Definiujemy jednowymiarową funkcję TF1 w pełnym zakresie kątowym -PI do PI
    TF1 *f_polar = new TF1("f_polar_fit", f_polar_circle, -TMath::Pi(), TMath::Pi(), 3);
    f_polar->SetParName(0, "R_true");
    f_polar->SetParName(1, "x0");
    f_polar->SetParName(2, "y0");

    // Parametry startowe dopasowane do przekazanego promienia spodziewanego (np. 4.5 dla BP lub 26.0 dla DC)
    f_polar->SetParameter(0, expected_r);
    f_polar->SetParameter(1, 0.0);
    f_polar->SetParameter(2, 1.5); // Sugerowany punkt startowy dla offsetu Y

    // Limity stabilizujące
    f_polar->SetParLimits(0, expected_r - 4.0, expected_r + 4.0);
    f_polar->SetParLimits(1, -5.0, 5.0);
    f_polar->SetParLimits(2, -5.0, 5.0);

    // Wywołanie stabilnego dopasowania 1D ("Q" - quiet, "R" - zakres z funkcji)
    // Wynikowa funkcja TF1 automatycznie doda się do listy funkcji profilu i zapisze w pliku ROOT!
    TFitResultPtr fit_result = prof->Fit(f_polar, "QRGS");

    if (fit_result->IsValid())
    {
      double r_true = f_polar->GetParameter(0);
      double x0 = f_polar->GetParameter(1);
      double y0 = f_polar->GetParameter(2);
    }

    return fit_result;
  }

  inline std::vector<double> iterative_polar_fit(TTreeReader &reader, const TTreeReaderArray<double> &KchrecFit, const TTreeReaderArray<double> &KnerecFit, const TTreeReaderArray<double> &Knerec, TTreeReaderArray<double> &Kchrec, TTreeReaderArray<double> &KchrecClosest, TTreeReaderArray<double> &trk1Fit, TTreeReaderArray<double> &trk2Fit, TTreeReaderValue<int> &mcflag, TTreeReaderValue<int> &mctruth, TTreeReaderValue<double> &Bx, TTreeReaderValue<double> &By, TTreeReaderValue<double> &Bz, TTreeReaderValue<double> &minv4gam, TTreeReaderArray<double> &KnerecSix, TProfile *profile, TString channel, TFitResultPtr &result, int iterations = 20)
  {
    TVector3 displacement_vector(0, 0, 0),
        kch_position(0, 0, 0),
        displacement_vector_err(0, 0, 0);

    TFitResultPtr tmp_result = nullptr;

    // Włączamy opcję "G" dla profilu na starcie, jeśli jeszcze nie jest włączona
    profile->SetErrorOption("G");

    for (int i = 0; i < iterations; i++)
    {
      // 1. KRYTYCZNE: Resetujemy reader na początek pliku przed każdą iteracją danych
      reader.SetEntry(0);

      // 2. KRYTYCZNE: Czyszczenie profilu z danych z poprzedniej iteracji
      profile->Reset();

      // Definicja lambdy dla cięcia
      double u0 = -92.9524, v0 = 7.1644, su = 33.7633, sv = 42.0163, rho = 0.556;
      double n_sigma_cut = 2.5;

      auto ellipse_cut = [=](double u, double v)
      {
        double du = u - u0;
        double dv = v - v0;
        double inv_rho2 = 1.0 / (1.0 - rho * rho);
        double dist2 = inv_rho2 * ((du * du) / (su * su) + (dv * dv) / (sv * sv) - 2.0 * rho * (du * dv) / (su * sv));
        return std::sqrt(std::max(0.0, dist2)) > n_sigma_cut;
      };

      // Pętla po zdarzeniach
      while (reader.Next())
      {
        bool global_cut = KloeTools::PassesGlobalCuts(KchrecFit, KnerecFit, Knerec, KchrecClosest, trk1Fit, trk2Fit,
                                                      *Bx, *By, *Bz, *minv4gam, KnerecSix[5], ellipse_cut);
        if (!global_cut)
          continue;

        kch_position.SetXYZ(Kchrec[6], Kchrec[7], Kchrec[8]);

        // Korekta pozycji DC o dotychczasowe, skumulowane przesunięcie
        kch_position -= displacement_vector;

        double rho_charged = kch_position.Perp();
        double phi_charged = std::atan2(kch_position.Y(), kch_position.X());
        double r_geom_xy = kch_position.Perp(); // To samo co kch_position.Perp()

        // Fiducial Volume zdefiniowane w układzie własnym komory (symetrycznie wokół 25.4 cm)
        if (rho_charged < 27 && rho_charged > 23)
        {
          if (KLOE::channName.at(*mctruth) == channel)
            profile->Fill(phi_charged, r_geom_xy);

          if (channel == "MC sum" && KLOE::channName.at(*mctruth) != "Data" && *mctruth != 0 && *mctruth != -1)
          {
            profile->Fill(phi_charged, r_geom_xy);
          }
        }
      }

      // Wykonujemy dopasowanie na oczyszczonym profilu
      tmp_result = KloeTools::fit_method_profile(profile, 25.4);

      // 3. BEZPIECZEŃSTWO: Sprawdzamy wynik fitu natychmiast po jego wykonaniu
      if (tmp_result->IsValid())
      {
        std::cout << "[+] Iteracja " << i << ": Dopasowanie profilu zakończone sukcesem." << std::endl;

        // Kumulujemy przesunięcie (poprawka z tej iteracji dąży do zera)
        displacement_vector += TVector3(tmp_result->Parameter(1), tmp_result->Parameter(2), 0);
        displacement_vector_err.SetXYZ(tmp_result->ParError(1), tmp_result->ParError(2), 0);

        std::cout << "    Aktualne skumulowane przesunięcie: (" << displacement_vector.X() << " ± " << displacement_vector_err.X() << ", "
                  << displacement_vector.Y() << " ± " << displacement_vector_err.Y() << ", 0)" << std::endl;
        std::cout << "    Aktualny dofitowany promień R: " << tmp_result->Parameter(0) << " ± " << tmp_result->ParError(0) << std::endl;

        result = tmp_result; // Nadpisujemy ostateczny wynik referencyjny
      }
      else
      {
        std::cerr << "[-] BŁĄD: Dopasowanie profilu w iteracji " << i << " nie powiodło się! Przerywam." << std::endl;
        break;
      }

      // Kryterium stopu: jeśli poprawki w danej iteracji są mniejsze niż błąd statystyczny, można wyjść wcześniej
      if (std::abs(tmp_result->Parameter(1)) < tmp_result->ParError(1) &&
          std::abs(tmp_result->Parameter(2)) < tmp_result->ParError(2))
      {
        std::cout << "[+] Osiągnięto zbieżność iteracyjną w kroku " << i << ". Przerywam dalsze pętle." << std::endl;
        break;
      }
    }

    reader.SetEntry(0); // Zostawiamy reader w stanie czystym na koniec
    return {displacement_vector.X(), displacement_vector.Y(), displacement_vector.Z(), displacement_vector_err.X(), displacement_vector_err.Y(), displacement_vector_err.Z()};
  }

  void quantile_limits(TH1 *h_data_clean, double range_min, double range_max, double &r_min, double &r_max)
  {
    // 1. Znajdujemy numery binów dla zdefiniowanego przez nas zakresu
    int bin_start = h_data_clean->GetXaxis()->FindBin(range_min);
    int bin_end = h_data_clean->GetXaxis()->FindBin(range_max);
    int n_bins_sub = bin_end - bin_start + 1;

    // 2. Tworzymy tymczasowy, lokalny histogram o dokładnych granicach, które nas interesują
    TH1D *h_sub = new TH1D("h_sub_temp", "Tymczasowy wycinek", n_bins_sub, range_min, range_max);

    // 3. Przepisujemy zawartość bin po binie TYLKO z wybranego obszaru
    for (int i = 0; i < n_bins_sub; ++i)
    {
      int source_bin = bin_start + i;
      h_sub->SetBinContent(i + 1, h_data_clean->GetBinContent(source_bin));
      h_sub->SetBinError(i + 1, h_data_clean->GetBinError(source_bin));
    }

    // 4. Definiujemy poziomy kwantyli (95% sygnału wewnątrz tego konkretnego okna)
    const int nq = 2;
    Double_t xq[nq] = {0.025, 0.975};
    Double_t yq[nq];

    // 5. Teraz GetQuantiles jest zmuszony liczyć kwantyle TYLKO z wyciętego fragmentu
    h_sub->GetQuantiles(nq, yq, xq);

    r_min = yq[0];
    r_max = yq[1];
    double efektywna_szerokosc = r_max - r_min;

    std::cout << "\n[+] POPRAWIONA METODA KWANTYLI (Wycinek [" << range_min << ", " << range_max << "] cm):" << std::endl;
    std::cout << "    Granica lewa (2.5%):  " << r_min << " cm" << std::endl;
    std::cout << "    Granica prawa (97.5%): " << r_max << " cm" << std::endl;
    std::cout << "    Efektywna szerokość ścianki w tym oknie: " << efektywna_szerokosc << " cm" << std::endl;

    // 6. Sprzątamy pamięć
    delete h_sub;
  }

  void plot_ratio_and_save(TH1 *h_numerator, TH1 *h_mc_sum, TH1 *h_regeneration, TString title, TString filename)
  {
    // 1. Sprawdzenie czy histogramy istnieją
    if (!h_numerator || !h_mc_sum || !h_regeneration)
    {
      std::cout << "[-] Błąd: Któryś z histogramów nie istnieje!" << std::endl;
      return;
    }

    // 2. Tworzenie histogramu ilorazu jako kopii licznika
    TH1D *h_ratio = (TH1D *)h_numerator->Clone("h_ratio_temp");
    h_ratio->SetTitle(title + ";#rho [cm];Dane / MC sum");

    // 3. Bezpośrednie dzielenie bin po binie (wraz z poprawnym przeliczeniem błędów)
    h_mc_sum->Add(h_regeneration, -1.0); // Dodajemy Regenerację do mianownika, żeby mieć pełną sumę MC
    h_ratio->Add(h_mc_sum, -1.0); // Odejmujemy sumę MC od licznika, żeby mieć tylko sygnał w liczniku
    h_ratio->Divide(h_regeneration);

    // 4. Stylizacja histogramu ilorazu
    h_ratio->SetMarkerStyle(20); // Pełne kropki
    h_ratio->SetMarkerSize(1.0);
    h_ratio->SetMarkerColor(kRed); // Czerwone punkty dla danych/MC
    h_ratio->SetLineColor(kBlack); // Czarne słupki błędów
    h_ratio->SetLineWidth(1);

    // Opcjonalnie: zawężamy oś Y wokół jedynki, żeby fluktuacje na krawędziach nie popsuły skali
    h_ratio->SetMinimum(0.0);
    h_ratio->SetMaximum(7.0);

    // 5. Rysowanie na Canvasie
    TCanvas *c_ratio = new TCanvas("c_ratio_temp", "Ratio Canvas", 800, 600);
    c_ratio->cd();
    gPad->SetGridy(); // Włączamy siatkę poziomą dla lepszej czytelności

    h_ratio->Draw("E1"); // Opcja "E1" rysuje punkty ze słupkami błędów

    // 6. Dodanie linii referencyjnej Y = 1.0
    double x_min = h_ratio->GetXaxis()->GetXmin();
    double x_max = h_ratio->GetXaxis()->GetXmax();
    TLine *line1 = new TLine(x_min, 1.0, x_max, 1.0);
    line1->SetLineColor(kGray + 2);
    line1->SetLineStyle(2); // Linia przerywana
    line1->SetLineWidth(2);
    line1->Draw("SAME");

    // 7. Zapis do pliku i sprzątanie pamięci
    c_ratio->SaveAs(filename);

    delete line1;
    delete c_ratio;
    delete h_ratio; // Usuwamy kopię, plik graficzny jest już bezpieczny na dysku
    std::cout << "[+] Wykres ilorazu zapisany do: " << filename << std::endl;
  }

  TH2 *plot_ratio_2D_and_save(TH2 *h_data, TH2 *h_mc_sum, TH2 *h_regeneration, TString title, TString filename)
  {
    // 1. Sprawdzenie czy wszystkie histogramy 2D istnieją
    if (!h_data || !h_mc_sum || !h_regeneration)
    {
      std::cout << "[-] Błąd: Któryś z histogramów 2D nie istnieje!" << std::endl;
      return nullptr;
    }

    // 2. Krok 1: h_mc_bkg = h_mc_sum - h_regeneration
    //    Odejmujemy Regenerację od sumy MC, żeby zostało samo tło (bez sygnału regeneracji)
    TH2D *h_mc_bkg = (TH2D *)h_mc_sum->Clone("h_mc_bkg_temp");
    h_mc_bkg->SetDirectory(0);
    h_mc_bkg->Add(h_regeneration, -1.0);

    // 3. Krok 2: h_data_signal = h_data - h_mc_bkg
    //    Odejmujemy tło MC od danych, żeby wyizolować sygnał regeneracji w danych
    TH2D *h_data_signal = (TH2D *)h_data->Clone("h_data_signal_temp");
    h_data_signal->SetDirectory(0);
    h_data_signal->Add(h_mc_bkg, -1.0);

    // 4. Krok 3: ratio = h_data_signal / h_regeneration
    //    Dzielimy oczyszczony sygnał w danych przez MC Regenerację
    TH2D *h_ratio_2D = (TH2D *)h_data_signal->Clone("h_ratio_2D_temp");
    h_ratio_2D->SetDirectory(0);

    TString x_title = h_data->GetXaxis()->GetTitle();
    TString y_title = h_data->GetYaxis()->GetTitle();
    h_ratio_2D->SetTitle(title + ";" + x_title + ";" + y_title);

    // Opcja "B" (Bayesowska) dba o poprawne błędy przy zerach w mianowniku
    h_ratio_2D->Divide(h_data_signal, h_regeneration, 1.0, 1.0, "B");

    // 5. Stylizacja i ustawienie skali kolorów
    h_ratio_2D->SetMinimum(0.0);
    h_ratio_2D->SetMaximum(7.0);
    h_ratio_2D->SetStats(kFALSE);

    // 6. Rysowanie na Canvasie
    TCanvas *c_ratio_2D = new TCanvas("c_ratio_2D_temp", "2D Ratio Canvas", 900, 800);
    c_ratio_2D->SetRightMargin(0.15);
    c_ratio_2D->cd();
    c_ratio_2D->SetLogz(0);
    gStyle->SetPalette(kRainBow);
    h_ratio_2D->Draw("COLZ");

    // 7. Zapis do pliku i czyszczenie pamięci
    c_ratio_2D->SaveAs(filename);
    delete c_ratio_2D;
    delete h_mc_bkg;
    delete h_data_signal;

    std::cout << "[+] Macierz poprawek 2D zapisana do: " << filename << std::endl;
    return h_ratio_2D;
  }

  double get_weight_2D(TH2* h_corr_matrix, double r_val, double dt_val)
  {
    // 1. Zabezpieczenie: Sprawdzamy, czy macierz w ogóle istnieje w pamięci
    if (!h_corr_matrix) {
      // Jeśli macierzy nie ma, zwracamy 1.0 (neutralna waga - nie zmienia wyniku)
      return 1.0; 
    }

    // 2. ROOT automatycznie znajduje globalny numer binu (X, Y) dla podanych wartości
    int bin2D = h_corr_matrix->FindBin(dt_val, r_val);

    // 3. Pobieramy czynnik poprawkowy (wagę) z tego binu
    double weight = h_corr_matrix->GetBinContent(bin2D);

    // 4. Fizyczne filtry bezpieczeństwa dla krawędzi macierzy (Low Statistics)
    // Jeśli waga jest zerowa, ujemna lub absurdalnie wysoka, zwracamy 1.0
    if (weight <= 0.0 || weight > 15.0) {
      return 1.0;
    }

    // 5. Zwracamy poprawnie wyznaczoną wagę
    return weight;
  }

  void initialize_matrices_from_file(TString input_filename, std::unordered_map<std::string, TH2 *> &h_corr_matrices)
  {
    // 2. Otwieramy plik ROOT w trybie tylko do odczytu ("READ")
    TFile* f_input = TFile::Open(input_filename, "READ");
    
    if (!f_input || f_input->IsZombie()) {
      std::cout << "[-] Błąd: Nie można otworzyć pliku " << input_filename << std::endl;
      return;
    }
    std::cout << "[+] Pomyślnie otwarto plik z wagami: " << input_filename << std::endl;

    // 3. Wczytujemy macierze za pomocą metody Get() 
    // WAŻNE: Musisz podać dokładnie taką samą nazwę (string), z jaką macierz została zapisana!
    h_corr_matrices["spherical"] = (TH2D*)f_input->Get("h_ratio_R_dt");
    h_corr_matrices["cylindrical"] = (TH2D*)f_input->Get("h_ratio_rho_dt");

    // 4. KRYTYCZNY KROK W ROOT: Odwiązanie histogramów od pliku (SetDirectory(0))
    // Jeśli tego nie zrobisz, w momencie zamknięcia pliku (f_input->Close()), 
    // Twoje macierze zostaną automatycznie usunięte z pamięci RAM!
    for (auto const& pair : h_corr_matrices) {
      std::string key = pair.first;
      TH2* matrix = pair.second;

      if (matrix) {
        matrix->SetDirectory(0); 
        std::cout << "[+] Załadowano macierz: " << key << " (Biny: " << matrix->GetNbinsX() << "x" << matrix->GetNbinsY() << ")" << std::endl;
      } else {
        std::cout << "[-] Ostrzeżenie: Nie znaleziono macierzy dla klucza: " << key << " w pliku ROOT!" << std::endl;
      }
    }

    // 5. Bezpiecznie zamykamy plik, macierze zostają w RAM-ie dzięki SetDirectory(0)
    f_input->Close();
    delete f_input;
  }

  void subtract_background_and_make_ratio(TH1* h_data_raw, TH1 *h_mc_sum, TH2 *h_data_dt_r, TH2 *h_mc_dt_r, double fit_min, double fit_max, TString output_root_file, TString output_png_filename)
  {
    if (!h_data_raw || !h_mc_sum) {
      std::cout << "[-] Blad: Ktorys z histogramow nie istnieje w subtract_background_and_make_ratio!" << std::endl;
      return;
    }

    // 1. Definiujemy funkcję fitu (eksponenta dla tla gazowego)
    TF1* f_bg = new TF1("f_bg", "expo", fit_min, fit_max);
    f_bg->SetParameter(1, -0.05); // Przybliżenie startowe (opadające tło)
    
    // Uruchamiamy dopasowanie tylko w zadanym zakresie tła ("R" - Range, "Q" - Quiet)
    h_data_raw->Fit(f_bg, "RQ");

    // 2. Tworzymy histogram dla SAMEGO WYMODELOWANEGO TŁA (w pełnym zakresie osi)
    TH1D* h_background_model = (TH1D*)h_data_raw->Clone("h_background_1D_model");
    h_background_model->SetDirectory(0);
    h_background_model->SetTitle("Wymodelowane tlo z fitu (Eksponenta)");
    h_background_model->Reset(); // Czyszczenie zawartości, zachowujemy tylko osie

    // 3. Tworzymy kopię Danych, z której odejmiemy tło
    TH1D* h_data_subtracted = (TH1D*)h_data_raw->Clone("h_data_1D_subtracted");
    h_data_subtracted->SetDirectory(0);
    h_data_subtracted->SetTitle("Dane po odjeciu tla");

    // 4. ODEJMOWANIE TŁA: Przechodzimy bin po binie i wypełniamy histogramy pośrednie
    for (int i = 1; i <= h_data_subtracted->GetNbinsX(); ++i) {
      double bin_center = h_data_subtracted->GetBinCenter(i);
      double raw_content = h_data_raw->GetBinContent(i);
      double raw_error = h_data_raw->GetBinError(i);

      // Wartość wykładnicza tła dla środka tego binu
      double bg_value = f_bg->Eval(bin_center);
      h_background_model->SetBinContent(i, bg_value);

      double new_content = raw_content - bg_value;
      if (new_content < 0) new_content = 0.0; // Bezpiecznik przed ujemną statystyką

      h_data_subtracted->SetBinContent(i, new_content);
      h_data_subtracted->SetBinError(i, raw_error); // Błąd statystyczny z surowych danych
    }

    // 5. TWORZENIE RATIO: Dzielimy oczyszczone Dane przez MC sum
    TH1D* h_ratio = (TH1D*)h_data_subtracted->Clone("h_ratio_subtracted_final");
    h_ratio->SetDirectory(0);
    h_ratio->SetTitle("Wagi czystego sygnalu (Dane bez tla / MC sum);#rho [cm];Stosunek");
    
    // Opcja "B" dla poprawnych błędów skorelowanych
    h_ratio->Divide(h_data_subtracted, h_mc_sum, 1.0, 1.0, "B");
    h_ratio->SetMinimum(0.0);
    h_ratio->SetMaximum(7.0);

    // 6. ZAPIS DO PLIKU .ROOT (Wszystkie kroki pośrednie lądują tutaj)
    TFile* f_out = TFile::Open(output_root_file, "RECREATE");
    if (f_out && !f_out->IsZombie()) {
      h_data_raw->Write("h_data_1D_raw");                 // 1. Surowe dane wejściowe
      f_bg->Write("f_background_fit_function");           // 2. Sama funkcja fitu
      h_background_model->Write("h_background_1D_model"); // 3. Wygenerowane na bazie fitu tło 1D
      h_data_subtracted->Write("h_data_1D_subtracted");   // 4. Czyste dane po odjęciu tła
      h_mc_sum->Write("h_mc_1D_sum");                     // 5. Wejściowe MC
      h_ratio->Write("h_ratio_subtracted_final");         // 6. Gotowe wagi
      
      f_out->Close();
      delete f_out;
      std::cout << "[+] Zapisano histogramy posrednie do pliku ROOT: " << output_root_file << std::endl;
    } else {
      std::cerr << "[-] Blad: Nie udalo sie utworzyc pliku ROOT: " << output_root_file << std::endl;
    }

    // 7. RYSOWANIE I WERYFIKACJA GRAFICZNA (.png)
    TCanvas* c_verify = new TCanvas("c_verify_temp", "Weryfikacja fitu i odejmowania tla", 1200, 600);
    c_verify->Divide(2, 1);

    c_verify->cd(1);
    h_data_raw->Draw("E");
    f_bg->SetLineColor(kRed);
    f_bg->SetLineWidth(2);
    f_bg->Draw("same"); // Pokazuje czerwoną krzywą tła na surowych danych

    c_verify->cd(2);
    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerColor(kBlue);
    h_ratio->Draw("E1"); // Pokazuje finalne wagi

    // Linia referencyjna Y = 1.0
    double x_min = h_ratio->GetXaxis()->GetXmin();
    double x_max = h_ratio->GetXaxis()->GetXmax();
    TLine* line1 = new TLine(x_min, 1.0, x_max, 1.0);
    line1->SetLineColor(kGray + 2);
    line1->SetLineStyle(2);
    line1->Draw("SAME");

    c_verify->SaveAs(output_png_filename);

    // 8. SPRZĄTANIE PAMIĘCI
    long old_level = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kSysError;

    delete line1;
    delete c_verify;
    delete f_bg;
    delete h_background_model;
    delete h_data_subtracted;

    gErrorIgnoreLevel = old_level;
    std::cout << "[+] Sukces: Wykres kontrolny zapisany do: " << output_png_filename << std::endl;
  }

 TH2* subtract_1D_fit_from_2D_matrix_and_ratio(
      TH2* h_data_raw, TH2* h_mc_sum, TH2* h_fit_reference,
      double rho_fit_min, double rho_fit_max, 
      double ratio_x_min, double ratio_x_max,
      double ratio_y_min, double ratio_y_max,
      TString output_root_file, TString output_png_filename)
  {
    if (!h_data_raw || !h_mc_sum || !h_fit_reference) {
      std::cout << "[-] Blad: Brak wymaganych histogramow 2D w subtract_1D_fit_from_2D_matrix_and_ratio!" << std::endl;
      return nullptr;
    }

    // 1. Rzutujemy wybrany HISTOGRAM REFERENCYJNY na os Y, aby wykonac stabilny fit 1D
    TH1D* h_ref_rho_1D = h_fit_reference->ProjectionY("h_ref_rho_1D_temp", 1, h_fit_reference->GetNbinsX());

    // 2. Fitujemy tlo eksponenta na rozkladzie 1D
    TF1* f_bg_1D = new TF1("f_bg_1D", "expo", rho_fit_min, rho_fit_max);
    f_bg_1D->SetParameter(1, -0.05);
    f_bg_1D->SetParLimits(1, -5.0, 0.0); 
    h_ref_rho_1D->Fit(f_bg_1D, "RQ");

    double global_slope = f_bg_1D->GetParameter(1);

    // 3. Przygotowujemy histogramy posrednie w 2D do zapisu
    TH2D* h_background_2D_model = (TH2D*)h_data_raw->Clone("h_background_2D_model");
    h_background_2D_model->SetDirectory(0);
    h_background_2D_model->SetTitle("Wymodelowane tlo 2D (z fitu 1D ref);#Delta t [ns];#rho [cm]");
    h_background_2D_model->Reset();

    TH2D* h_data_2D_subtracted = (TH2D*)h_data_raw->Clone("h_data_2D_subtracted");
    h_data_2D_subtracted->SetDirectory(0);
    h_data_2D_subtracted->SetTitle("Dane 2D po odjeciu tla;#Delta t [ns];#rho [cm]");

    // 4. SKALOWANIE I ODEJMOWANIE TLA W 2D (Pasek po pasku Delta t)
    for (int bx = 1; bx <= h_data_raw->GetNbinsX(); ++bx) {
      TH1D* h_slice_rho = h_data_raw->ProjectionY("h_slice_rho_temp", bx, bx);
      
      TF1* f_slice_bg = new TF1("f_slice_bg", "expo", rho_fit_min, rho_fit_max);
      f_slice_bg->FixParameter(1, global_slope);
      h_slice_rho->Fit(f_slice_bg, "RQ");

      for (int by = 1; by <= h_data_raw->GetNbinsY(); ++by) {
        double rho_val = h_data_raw->GetYaxis()->GetBinCenter(by);
        double raw_content = h_data_raw->GetBinContent(bx, by);
        double raw_error   = h_data_raw->GetBinError(bx, by);

        double bg_val = f_slice_bg->Eval(rho_val);
        h_background_2D_model->SetBinContent(bx, by, bg_val);

        double new_content = raw_content - bg_val;
        if (new_content < 0.0) new_content = 0.0; 

        h_data_2D_subtracted->SetBinContent(bx, by, new_content);
        h_data_2D_subtracted->SetBinError(bx, by, raw_error);
      }

      delete h_slice_rho;
      delete f_slice_bg;
    }

    // 5. TWORZENIE FINALNEJ MACIERZY RATIO 2D (Z bezpiecznikiem dla MC == 0)
    TH2D* h_ratio_2D = (TH2D*)h_data_2D_subtracted->Clone("h_ratio_rho_dt");
    h_ratio_2D->SetDirectory(0);
    h_ratio_2D->SetTitle("Finalna macierz wag 2D (Dane bez tla / MC sum);#Delta t [ns];#rho [cm]");
    h_ratio_2D->Reset(); // Resetujemy zawartość, napełnimy ją ręcznie bin po binie

    for (int bx = 1; bx <= h_data_2D_subtracted->GetNbinsX(); ++bx) {
      double x_val = h_data_2D_subtracted->GetXaxis()->GetBinCenter(bx);
      for (int by = 1; by <= h_data_2D_subtracted->GetNbinsY(); ++by) {
        double y_val = h_data_2D_subtracted->GetYaxis()->GetBinCenter(by);
        
        // Sprawdzamy czy bin mieści się w zdefiniowanym przez Ciebie oknie analizy
        if (x_val < ratio_x_min || x_val > ratio_x_max || y_val < ratio_y_min || y_val > ratio_y_max) {
          h_ratio_2D->SetBinContent(bx, by, -1.0);
          h_ratio_2D->SetBinError(bx, by, 0.0);
          continue;
        }

        double data_val = h_data_2D_subtracted->GetBinContent(bx, by);
        double data_err = h_data_2D_subtracted->GetBinError(bx, by);
        double mc_val   = h_mc_sum->GetBinContent(bx, by);
        double mc_err   = h_mc_sum->GetBinError(bx, by);

        // KRYTYCZNY BEZPIECZNIK: Jeśli w MC nie ma zliczeń, ustawiamy wagę na zero
        if (mc_val <= 0.0) {
          h_ratio_2D->SetBinContent(bx, by, -1.0);
          h_ratio_2D->SetBinError(bx, by, 0.0);
        } 
        else {
          // Standardowe dzielenie
          double ratio = data_val / mc_val;
          
          // Klasyczne propagowanie błędów dla niezależnych wartości (zamiennik opcji "B")
          double ratio_err = 0.0;
          if (data_val > 0.0) {
            ratio_err = ratio * std::sqrt(std::pow(data_err/data_val, 2) + std::pow(mc_err/mc_val, 2));
          }

          h_ratio_2D->SetBinContent(bx, by, ratio);
          h_ratio_2D->SetBinError(bx, by, ratio_err);
        }
      }
    }

    h_ratio_2D->SetMinimum(0.0);
    h_ratio_2D->SetMaximum(15.0);
    h_ratio_2D->SetStats(kFALSE);

    // 6. ZAPIS WSZYSTKICH KROKÓW POŚREDNICH DO PLIKU .ROOT
    TFile* f_out = TFile::Open(output_root_file, "RECREATE");
    if (f_out && !f_out->IsZombie()) {
      h_data_raw->Write("h_data_2D_raw");                     
      h_fit_reference->Write("h_fit_2D_reference_used");      
      h_ref_rho_1D->Write("h_data_1D_ref_projection");       
      f_bg_1D->Write("f_background_1D_fit");                  
      h_background_2D_model->Write("h_background_2D_model");   
      h_data_2D_subtracted->Write("h_data_2D_subtracted");   
      h_mc_sum->Write("h_mc_2D_sum");                         
      h_ratio_2D->Write("h_ratio_rho_dt");                    
      
      f_out->Close();
      delete f_out;
      std::cout << "[+] Zapisano pelna historie 2D do pliku: " << output_root_file << std::endl;
    }

    // 7. SZYBKI ZAPIS GRAFICZNY PANELU 2D (.png)
    TCanvas* c_view = new TCanvas("c_view_temp", "Kontrola 1D Fit -> 2D Subtraction", 1600, 1200);
    c_view->Divide(2, 2);
    gStyle->SetPalette(kRainBow);

    c_view->cd(1); gPad->SetRightMargin(0.15); h_data_raw->Draw("COLZ");
    c_view->cd(2); gPad->SetRightMargin(0.15); h_background_2D_model->Draw("COLZ");
    c_view->cd(3); gPad->SetRightMargin(0.15); h_data_2D_subtracted->Draw("COLZ");
    c_view->cd(4); gPad->SetRightMargin(0.15); h_ratio_2D->Draw("COLZ");
    
    c_view->SaveAs(output_png_filename);

    // 8. SPRZĄTANIE RAMu
    long old_level = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kSysError;
    delete c_view;
    delete f_bg_1D;
    delete h_ref_rho_1D;
    delete h_background_2D_model;
    delete h_data_2D_subtracted;
    gErrorIgnoreLevel = old_level;

    return h_ratio_2D;
  }
}

#endif