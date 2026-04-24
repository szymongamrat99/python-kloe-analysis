import ROOT
from lib.hist.hist_model_loader import HistModel


class MultiChannelHist:
    def __init__(self, rdf, channName, channColor, hist_model=None,
                 var=None, model=None, base_filter="", global_filter="", scale_map={}, config=None):
        """
        rdf           - ROOT.RDataFrame (or filtered node)
        channName     - dict {mctruth_id: channel_name}
        channColor    - dict {channel_name: ROOT color int}
        hist_model    - HistModel instance (preferred, from HistModelLoader)
        var           - column name (fallback if hist_model not given)
        model         - tuple (prefix, title, nbins, xmin, xmax) (fallback)
        base_filter   - per-histogram filter (overrides hist_model.base_filter)
        global_filter - filter applied to ALL histograms (combined with base_filter)
        """
        self.channName = channName
        self.channColor = channColor
        self.scale_map = scale_map
        self.histos = {}

        # Resolve from HistModel or legacy tuple
        if hist_model is not None:
            self.var = hist_model.var
            self.base_filter = base_filter or hist_model.base_filter
            self._hist_model = hist_model
        else:
            self.var = var
            self.base_filter = base_filter
            prefix, title, nbins, xmin, xmax = model
            self._hist_model = HistModel(prefix, {
                "var": var, "title": title,
                "nbins": nbins, "xmin": xmin, "xmax": xmax,
            })

        self.global_filter = global_filter

        # Combine global + per-histogram filters
        filters = [f for f in [global_filter, self.base_filter] if f]
        combined_filter = " && ".join(f"({f})" for f in filters) if filters else ""

        rdf_mc = rdf.Filter("mcflag == 1 && mctruth != 0 && mctruth != -1")
        rdf_data = rdf.Filter("mcflag == 0")

        # Count per-channel BEFORE cuts (for efficiency calculation)
        self._counts_before = {}
        for key, chann in channName.items():
            if chann == "MC sum" or chann == "Data":
                continue
            self._counts_before[chann] = rdf_mc.Filter(f"mctruth == {key}").Count()

        if combined_filter:
            rdf_mc = rdf_mc.Filter(combined_filter)
            rdf_data = rdf_data.Filter(combined_filter)

        for key, chann in channName.items():
            if chann == "MC sum":
                continue

            hmodel = self._hist_model.make_root_model(suffix=chann)

            # --- NOWA LOGIKA ZABEZPIECZAJĄCA ---
            # Jeśli model to 'resolution', akceptujemy tylko MC (mctruth 0 lub 1)
            is_resolution = getattr(self._hist_model, 'resolution', False)
            
            if is_resolution and chann == "Data":
                # Tworzymy pusty histogram dla 'Data', żeby uniknąć błędów KeyErrro w draw()
                # ale nie przypisujemy mu żadnych danych z rdf_data
                self.histos[chann] = self._hist_model.make_empty_th(suffix=chann)
                continue
            # -----------------------------------

            if chann == "Data":
                self.histos[chann] = rdf_data.Histo1D(hmodel, self.var)
            elif chann == "Signal" and self._hist_model.weight:
                # Weighted Signal: fill with weight, then rescale to original entries
                rdf_sig = rdf_mc.Filter(f"mctruth == {key}")
                self.histos[chann] = rdf_sig.Histo1D(hmodel, self.var, self._hist_model.weight)
                self._signal_key = key
                self._signal_rdf = rdf_sig  # keep ref for entry count
            elif chann == "Regeneration":
                weight_far_left = self.scale_map["Regeneration"][0]
                weight_near_left = self.scale_map["Regeneration"][1]
                weight_near_right = self.scale_map["Regeneration"][2]
                weight_far_right = self.scale_map["Regeneration"][3]

                split_histos = []
                split_histos.append(rdf_mc.Define("weight_col", f"{weight_far_left}").Filter(f"mctruth == {key} && deltaT <= -30").Histo1D(hmodel, self.var, "weight_col"))
                split_histos.append(rdf_mc.Define("weight_col", f"{weight_near_left}").Filter(f"mctruth == {key} && deltaT > -30 && deltaT <= 0").Histo1D(hmodel, self.var, "weight_col"))
                split_histos.append(rdf_mc.Define("weight_col", f"{weight_near_right}").Filter(f"mctruth == {key} && deltaT > 0 && deltaT <= 30").Histo1D(hmodel, self.var, "weight_col"))
                split_histos.append(rdf_mc.Define("weight_col", f"{weight_far_right}").Filter(f"mctruth == {key} && deltaT > 30").Histo1D(hmodel, self.var, "weight_col"))

                split_histos[0].Add(split_histos[1].GetPtr())
                split_histos[0].Add(split_histos[2].GetPtr())
                split_histos[0].Add(split_histos[3].GetPtr())

                self.histos[chann] = split_histos[0]
            else:
                weight = self.scale_map.get(chann, 1.0)
                self.histos[chann] = rdf_mc.Define("weight_col", f"{weight}").Filter(f"mctruth == {key}").Histo1D(hmodel, self.var, "weight_col")

        self.histos["MC sum"] = self._hist_model.make_empty_th(suffix="MC sum")
        self._mc_sum_built = False

    def build(self):
        """Trigger RDF evaluation, build MC sum, and scale. Call before save/draw."""
        self._rescale_weighted_signal()
        self._build_mc_sum()
        return self

    def _resolution_fitter(self):
        """Custom fitter for resolution histograms."""

        # Get the resolution histogram for Signal
        h = self.histos["MC sum"]

        if not h:
          return (None,) * 9

        hp = h.GetPtr() if hasattr(h, 'GetPtr') else h

        # Fit with a Gaussian (or any other suitable function)
        double_gaus = ROOT.TF1("double_gaus", ROOT.double_gaus, hp.GetXaxis().GetXmin(), hp.GetXaxis().GetXmax(), 5)
        double_gaus.SetParameters(hp.GetMaximum(), hp.GetRMS(), hp.GetMaximum()/2, hp.GetRMS()*2, hp.GetMean())
        double_gaus.SetParLimits(1, 0.01 * hp.GetRMS(), 10 * hp.GetRMS())  # sigma_core limits
        double_gaus.SetParLimits(3, 0.01 * hp.GetRMS(), 20 * hp.GetRMS())  # sigma_tail limits

        double_gaus.SetParNames("A_core", "sigma_core", "A_tail", "sigma_tail", "mean")
        fit_result = hp.Fit(double_gaus, "S")
        if fit_result.IsValid():
            chi2 = fit_result.Chi2()
            ndf = fit_result.Ndf()
            mean = fit_result.Parameter(4)
            mean_err = fit_result.ParError(4)
            sigma_core = fit_result.Parameter(1)
            sigma_core_err = fit_result.ParError(1)
            sigma_tail = fit_result.Parameter(3)
            sigma_tail_err = fit_result.ParError(3)
            print(f"Fitted resolution: mean={mean:.3f} ± {mean_err:.3f}, sigma_core={sigma_core:.3f} ± {sigma_core_err:.3f}, sigma_tail={sigma_tail:.3f} ± {sigma_tail_err:.3f}, chi2/ndf={chi2}/{ndf}")

            # 5. Bezpieczne zwracanie wyników
        if fit_result and fit_result.IsValid():
          return (
            fit_result.Parameter(1), # sigma_core
            fit_result.ParError(1),  # sigma_core_err
            fit_result.Parameter(3), # sigma_tail
            fit_result.ParError(3),  # sigma_tail_err
            fit_result.Parameter(4), # mean
            fit_result.ParError(4),  # mean_err
            fit_result.Chi2(),       # chi2
            fit_result.Ndf(),        # ndf
            double_gaus
          )
        else:
          print(f"Warning: Fit failed for {hp.GetName()}. Returning null results.")
          return (None,) * 9

    def _rescale_weighted_signal(self):
        """If Signal was weighted, rescale so integral matches original entry count."""
        if not getattr(self, '_signal_rdf', None):
            return
        h = self.histos["Signal"]
        hp = h.GetPtr() if hasattr(h, 'GetPtr') else h
        integral = hp.Integral(0, hp.GetNbinsX() + 1)
        if integral <= 0:
            return
        original_entries = self._signal_rdf.Count().GetValue()
        if original_entries > 0:
            weight = self.scale_map.get("Signal", 1.0)
            hp.Scale(weight * original_entries / integral)

    def _build_mc_sum(self):
        if self._mc_sum_built:
            return
        mc_sum = self.histos["MC sum"]
        mc_sum.Reset()
        for chann, hist in self.histos.items():
            if chann != "Data" and chann != "MC sum":
                mc_sum.Add(hist.GetPtr() if hasattr(hist, 'GetPtr') else hist)
        self._mc_sum_built = True

    def _scale_mc_to_data(self):
        if "Data" not in self.histos or "MC sum" not in self.histos:
            return

        data_entries = self.histos["Data"].GetEntries()
        mc_entries = self.histos["MC sum"].GetEntries()

        if mc_entries == 0:
            return

        scale = data_entries / mc_entries
        for chann, hist in self.histos.items():
            if chann != "Data":
                h = hist.GetPtr() if hasattr(hist, 'GetPtr') else hist
                h.Scale(scale)

    def draw(self, canvas_name=None, width=850, height=850,
             scale_mc=True, save_as=None):

        self._build_mc_sum()
        if scale_mc:
            self._scale_mc_to_data()

        hm = self._hist_model

        if canvas_name is None:
            canvas_name = f"c_{hm.name}"

        canva = ROOT.TCanvas(canvas_name, "", width, height)

        # Stats box
        ROOT.gStyle.SetOptStat(1 if hm.show_stats else 0)

        # Log scales with safety checks
        if hm.log_x:
            if hm.xmin > 0:
                canva.SetLogx(1)
            else:
                print(f"[warn] {hm.name}: log_x requested but xmin={hm.xmin} <= 0, fallback to linear")
        if hm.log_y:
            canva.SetLogy(1)

        # Find channel with highest peak for drawing first
        peaks = {}
        for chann, hist in self.histos.items():
            h = hist.GetPtr() if hasattr(hist, 'GetPtr') else hist
            peaks[chann] = h.GetMaximum()

        draw_order = sorted(peaks, key=peaks.get, reverse=True)

        # Legend from model config
        lp = hm.legend_pos if hasattr(hm, 'legend_pos') else [0.65, 0.65, 0.88, 0.88]
        legend = ROOT.TLegend(*lp)
        legend.SetBorderSize(0)

        first = True
        for chann in draw_order:
            hist = self.histos[chann]
            h = hist.GetPtr() if hasattr(hist, 'GetPtr') else hist
            color = self.channColor.get(chann, ROOT.kBlack)
            h.SetLineColor(color)

            is_data = (chann == "Data")

            if first:
                h.Draw("PE1" if is_data else "HIST")

                # Y-axis range
                if hm.ymax_display is not None:
                    h.SetMaximum(hm.ymax_display)
                if hm.ymin is not None and hm.dim != 2:
                    if hm.log_y and hm.ymin <= 0:
                        print(f"[warn] {hm.name}: log_y with ymin={hm.ymin} <= 0, setting ymin=0.1")
                        h.SetMinimum(0.1)
                    else:
                        h.SetMinimum(hm.ymin)
                elif hm.log_y:
                    h.SetMinimum(0.1)

                if hm.resolution and chann == "MC sum":
                    result = self._resolution_fitter()

                    if result[0] is not None:

                      if result[0] < result[2]:
                          sigma_core = result[0]
                          sigma_tail = result[2]
                          sigma_core_err = result[1]
                          sigma_tail_err = result[3]
                      else:
                          sigma_core = result[2]
                          sigma_tail = result[0]
                          sigma_core_err = result[3]
                          sigma_tail_err = result[1]

                      mean = result[4]
                      mean_err = result[5]

                      chi2 = result[6]
                      ndf = result[7]

                      # Tworzymy tabelkę z parametrami (TPaveText)
                      # Współrzędne: x1, y1, x2, y2 (w jednostkach NDC: 0-1)
                      # Umieszczamy ją poniżej legendy (legenda zazwyczaj kończy się na ~0.65)
                      pave = ROOT.TPaveText(0.65, 0.40, 0.88, 0.63, "NDC")
                      pave.SetFillColor(0)
                      pave.SetTextAlign(12) # wyrównanie do lewej
                      pave.SetTextFont(42)
                      pave.SetTextSize(0.025)
                      pave.SetBorderSize(1)

                      pave.AddText(f"Mean: {mean:.3f} #pm {mean_err:.3f}")
                      pave.AddText(f"#sigma_{{core}}: {sigma_core:.3f} #pm {sigma_core_err:.3f}")
                      pave.AddText(f"#sigma_{{tail}}: {sigma_tail:.3f} #pm {sigma_tail_err:.3f}")
                      pave.AddText(f"#chi^{{2}} / ndf: {chi2:.1f} / {ndf}")
                      
                      pave.Draw()

                      result[-1].SetLineColor(ROOT.kRed)
                      result[-1].SetLineWidth(2)
                      result[-1].Draw("SAME")
                    else:
                      print(f"Resolution fit failed for {chann}. No curve will be drawn.")

                first = False
            else:
                h.Draw("PE1SAME" if is_data else "HISTSAME")

                if hm.resolution and chann == "MC sum":
                    result = self._resolution_fitter()

                    if result[0] is not None:

                      if result[0] < result[2]:
                          sigma_core = result[0]
                          sigma_tail = result[2]
                          sigma_core_err = result[1]
                          sigma_tail_err = result[3]
                      else:
                          sigma_core = result[2]
                          sigma_tail = result[0]
                          sigma_core_err = result[3]
                          sigma_tail_err = result[1]

                      mean = result[4]
                      mean_err = result[5]

                      chi2 = result[6]
                      ndf = result[7]

                      # Tworzymy tabelkę z parametrami (TPaveText)
                      # Współrzędne: x1, y1, x2, y2 (w jednostkach NDC: 0-1)
                      # Umieszczamy ją poniżej legendy (legenda zazwyczaj kończy się na ~0.65)
                      pave = ROOT.TPaveText(0.65, 0.40, 0.88, 0.63, "NDC")
                      pave.SetFillColor(0)
                      pave.SetTextAlign(12) # wyrównanie do lewej
                      pave.SetTextFont(42)
                      pave.SetTextSize(0.025)
                      pave.SetBorderSize(1)

                      pave.AddText(f"Mean: {mean:.3f} #pm {mean_err:.3f}")
                      pave.AddText(f"#sigma_{{core}}: {sigma_core:.3f} #pm {sigma_core_err:.3f}")
                      pave.AddText(f"#sigma_{{tail}}: {sigma_tail:.3f} #pm {sigma_tail_err:.3f}")
                      pave.AddText(f"#chi^{{2}} / ndf: {chi2:.1f} / {ndf}")
                      
                      pave.Draw()

                      result[-1].SetLineColor(ROOT.kRed)
                      result[-1].SetLineWidth(2)
                      result[-1].Draw("SAME")
                    else:
                      print(f"Resolution fit failed for {chann}. No curve will be drawn.")

            legend.AddEntry(h, f"{chann} ({int(h.Integral(0, h.GetNbinsX() + 1))})", "lp" if is_data else "l")

        legend.Draw()

        if save_as:
            canva.SaveAs(save_as)

        return canva, legend
