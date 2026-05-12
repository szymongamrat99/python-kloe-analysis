import ROOT
import numpy as np

class EllipseFitter:
    def __init__(self, rdf, var_u, var_v):
        self.rdf = rdf
        self.u = var_u
        self.v = var_v
        self.params = {}

    def fit(self, n_iterations=3, sigma_cut=3.0):
      """
      Oblicza parametry elipsy z iteracyjnym usuwaniem wartości odstających (robust).
      n_iterations: liczba powtórzeń (zazwyczaj 2-3 wystarczą)
      sigma_cut: próg odległości Mahalanobisa, powyżej którego punkt jest uznawany za ogon
      """
      # Pobieramy dane raz na początku
      cols = self.rdf.AsNumpy(columns=[self.u, self.v])
      u_data = cols[self.u]
      v_data = cols[self.v]
      
      # Inicjalizujemy maskę (na początku bierzemy wszystkie punkty)
      mask = np.ones(len(u_data), dtype=bool)

      for i in range(n_iterations):
          u_curr = u_data[mask]
          v_curr = v_data[mask]

          if len(u_curr) < 2:
              print(f"Warning: Too few points left after iteration {i}")
              break

          # 1. Obliczamy parametry dla aktualnego zbioru (rdzenia)
          u0 = np.mean(u_curr)
          v0 = np.mean(v_curr)
          cov = np.cov(u_curr, v_curr)
          
          su = np.sqrt(cov[0, 0])
          sv = np.sqrt(cov[1, 1])
          # Zabezpieczenie przed dzieleniem przez zero
          rho = cov[0, 1] / (su * sv) if (su * sv) > 0 else 0

          self.params = {
              "u0": float(u0),
              "v0": float(v0),
              "su": float(su),
              "sv": float(sv),
              "rho": float(rho)
          }

          # 2. Obliczamy odległości Mahalanobisa dla WSZYSTKICH punktów względem obecnej elipsy
          # d^2 = 1/(1-rho^2) * [ (du/su)^2 + (dv/sv)^2 - 2*rho*du*dv/(su*sv) ]
          du = u_data - u0
          dv = v_data - v0
          
          # Współczynnik normalizacyjny dla korelacji
          inv_rho_sq = 1.0 / (1.0 - rho**2)
          dist_sq = inv_rho_sq * (
              (du / su)**2 + 
              (dv / sv)**2 - 
              2 * rho * (du * dv) / (su * sv)
          )

          # 3. Aktualizujemy maskę na kolejną iterację (wycinamy ogony)
          mask = dist_sq < (sigma_cut**2)
          
          n_removed = len(u_data) - np.sum(mask)
          # print(f"Iteration {i+1}: Removed {n_removed} outliers beyond {sigma_cut} sigma")

      # Zaokrąglamy wyniki na samym końcu dla estetyki w JSON
      for key in self.params:
          self.params[key] = round(self.params[key], 4)

      return self.params

    def draw_fit(self, n_sigma=1.0, output_path="ellipse_fit.png"):
        if not self.params:
            self.fit()

        p = self.params
        
        # 1. Pobieramy histogram (Lazy)
        h2_ptr = self.rdf.Histo2D(
            (f"h2_{self.u}_{self.v}", ";m^{inv}_{4#gamma} [MeV/c^{2}];m^{inv}_{6#gamma} [MeV/c^{2}]", 
             201, -400.0, 500.0,
             201, -400.0, 500.0),
            self.u, self.v
        )
        
        # Wyzwalany jest RunGraphs lub GetPtr
        h2 = h2_ptr.GetPtr()

        c = ROOT.TCanvas("c_fit", "Ellipse Fit", 750, 750)
        h2.GetXaxis().SetRangeUser(-400.0, 500.0)
        h2.GetYaxis().SetRangeUser(-400.0, 500.0)
        h2.Draw("COLZ")

        # --- POPRAWIONA MATEMATYKA ELIPSY ---
        # Macierz kowariancji
        cov_uv = p['rho'] * p['su'] * p['sv']
        
        # Obliczamy wartości własne (Eigenvalues) macierzy 2x2
        # lambda = (tr(V) +/- sqrt(tr(V)^2 - 4*det(V))) / 2
        tr_v = p['su']**2 + p['sv']**2
        det_v = (p['su']**2 * p['sv']**2) - (cov_uv**2)
        
        # Zabezpieczenie przed ujemnym wyznacznikiem (numeryczne)
        delta = np.sqrt(max(0, tr_v**2 - 4*det_v))
        lambda1 = (tr_v + delta) / 2
        lambda2 = (tr_v - delta) / 2
        
        # Półosie to pierwiastki z wartości własnych
        a = np.sqrt(lambda1) * n_sigma
        b = np.sqrt(lambda2) * n_sigma

        # Kąt obrotu (Angle of the eigenvector)
        # theta = 0.5 * arctan(2 * cov(u,v) / (var(u) - var(v)))
        angle = 0.5 * np.arctan2(2 * cov_uv, p['su']**2 - p['sv']**2)
        angle_deg = np.degrees(angle)

        # 3. Rysowanie
        el = ROOT.TEllipse(p['u0'], p['v0'], a, b, 0, 360, angle_deg)
        el.SetFillStyle(0)
        el.SetLineColor(ROOT.kBlack)
        el.SetLineWidth(3)
        el.Draw("same")

        c.SaveAs(output_path)
        return c, h2, el

    def generate_json_snippet(self):
        """Generuje gotowy tekst do wklejenia do JSONa."""
        if not self.params:
            self.fit()
        
        p = self.params
        var_str = f"getCorrelatedEllipseDist({self.u}, {self.v}, {p['u0']}, {p['v0']}, {p['su']}, {p['sv']}, {p['rho']})"
        
        return {
            "variable": var_str,
            "value": 1.0,
            "operator": "<",
            "description": f"Optimized ellipse for {self.u} vs {self.v}"
        }