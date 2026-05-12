import ROOT
import numpy as np

class KaonLifetimeFitter:
    def __init__(self, rdf, time_var="deltaT_MC"):
        self.rdf = rdf
        self.time_var = time_var
        self.hist = None
        self.func = None

    def fit(self, t_range=15.0, n_bins=200, exclude_width=0.15):
        """
        t_range: zakres od -max do +max
        exclude_width: szerokość (ns) wycięcia wokół zera (od -width do +width)
        """
        # 1. Rezerwacja histogramu (symetryczny zakres)
        h_ptr = self.rdf.Histo1D(
            (f"h_life_{self.time_var}", "Lifetime Fit;#Deltat^{MC} [ns];Events", 
             n_bins, -t_range, t_range), self.time_var
        )
        self.hist = h_ptr.GetPtr()

        # 2. Definicja funkcji z wykluczeniem środka (RejectPoint)
        # par[0]: Norm, par[1]: tauS, par[2]: Ratio, par[3]: tauL
        def fit_logic(x, par):
            if abs(x[0]) < exclude_width:
                ROOT.TF1.RejectPoint()
                return 0
            return par[0] * (ROOT.TMath.Exp(-abs(x[0])/par[1]) + par[2] * ROOT.TMath.Exp(-abs(x[0])/par[3]))

        self.func = ROOT.TF1("f_life", fit_logic, -t_range, t_range, 4)
        self.func.SetParNames("Norm", "tau_S", "Ratio_L_S", "tau_L")

        # 3. KROK 1: Wstępne dopasowanie KL (na dalekich ogonach)
        # Szukamy tylko poziomu "tła" KL tam, gdzie KS już wygasło (> 2 ns)
        print("[Fit] Step 1: Estimating KL background...")
        self.func.SetParameters(self.hist.GetMaximum(), 0.089, 0.001, 51.16)
        self.func.FixParameter(1, 0.089) # Tymczasowo mrozimy KS
        self.hist.Fit(self.func, "L Q N R", "", 2.0, t_range) 
        
        # 4. KROK 2: Fit właściwy (uwzględniający KS i KL)
        print("[Fit] Step 2: Full fit with central exclusion...")
        self.func.ReleaseParameter(1)
        self.func.SetParLimits(1, 0.01, 0.12)  # tau_S
        self.func.SetParLimits(3, 40.0, 70.0)  # tau_L
        self.func.SetParLimits(2, 1e-7, 0.5)   # Ratio
        
        # Ustawienie maksimum na podstawie najwyższego binu poza wycięciem
        bin_limit = self.hist.FindBin(exclude_width)
        self.max_visible = self.hist.GetBinContent(bin_limit)
        self.func.SetParameter(0, self.max_visible)

        fit_res = self.hist.Fit(self.func, "S L M Q R")
        
        return fit_res

    def draw(self, output_name="fit_results.png"):
        canvas = ROOT.TCanvas("c_life", "", 900, 700)
        canvas.SetLogy()
        
        # Skalowanie osi Y, żeby nie widzieć zerowego binu w środku
        # Znajdź najwyższy bin POZA obszarem wyciętym (exclude_width)
        bin_ex_plus = self.hist.FindBin(0.15) # Zakładając exclude_width = 0.15
        new_max = self.hist.GetBinContent(bin_ex_plus)
        
        # Ustaw zakres osi Y tak, by ten "drugi" szczyt był na górze
        # Mnożymy przez 5-10, żeby mieć miejsce na legendę w skali log
        self.hist.SetMaximum(new_max * 10)
        self.hist.SetMinimum(1E2) # Żeby nie widzieć "szumu" na dnie
        
        self.hist.SetMarkerStyle(20)
        self.hist.SetMarkerSize(0.6)
        self.hist.Draw("E1")
        
        self.func.SetLineColor(ROOT.kRed)
        self.func.SetLineWidth(3)
        self.func.Draw("same")

        self._draw_components()

        # Dodanie linii oznaczających wycięty obszar
        self._draw_exclusion_lines()

        canvas.SaveAs(output_name)

    def _draw_components(self):
        p = self.func.GetParameters()
        xmin, xmax = self.func.GetXmin(), self.func.GetXmax()
        
        # Komponent KS (niebieski)
        self._f_ks = ROOT.TF1("f_ks", "[0]*exp(-abs(x)/[1])", xmin, xmax)
        self._f_ks.SetParameters(p[0], p[1])
        self._f_ks.SetLineStyle(2)
        self._f_ks.SetLineColor(ROOT.kBlue)
        self._f_ks.Draw("same")

        # Komponent KL (zielony)
        self._f_kl = ROOT.TF1("f_kl", "[0]*[2]*exp(-abs(x)/[3])", xmin, xmax)
        self._f_kl.SetParameters(p[0], p[1], p[2], p[3])
        self._f_kl.SetLineStyle(2)
        self._f_kl.SetLineColor(ROOT.kGreen+2)
        self._f_kl.Draw("same")

    def _draw_exclusion_lines(self):
        """Rysuje pionowe linie w miejscu, gdzie zaczyna się fit."""
        # Pobieramy t_min/t_max pośrednio z logiki RejectPoint (nasze exclude_width)
        # Dla uproszczenia załóżmy, że linie rysujemy tam, gdzie funkcja "wraca"
        pass # Można tu dodać TLine