import ROOT


class MultiChannelHist:
    def __init__(self, rdf, channName, channColor, var, model, base_filter=""):
        """
        rdf        - ROOT.RDataFrame (or filtered node)
        channName  - dict {mctruth_id: channel_name}
        channColor - dict {channel_name: ROOT color int}
        var        - column name to histogram (e.g. "Chi2SignalKinFit")
        model      - tuple (name_prefix, title, nbins, xmin, xmax)
        base_filter - optional additional filter applied to all channels
        """
        self.channName = channName
        self.channColor = channColor
        self.var = var
        self.model = model
        self.histos = {}

        rdf_mc = rdf.Filter("mcflag == 1 && mctruth != 0 && mctruth != -1")
        rdf_data = rdf.Filter("mcflag == 0")

        if base_filter:
            rdf_mc = rdf_mc.Filter(base_filter)
            rdf_data = rdf_data.Filter(base_filter)

        prefix, title, nbins, xmin, xmax = model

        for key, chann in channName.items():
            if chann == "MC sum":
                continue

            hmodel = ROOT.RDF.TH1DModel(
                f"{prefix}_{chann}", title, nbins, xmin, xmax
            )

            if chann == "Data":
                self.histos[chann] = rdf_data.Histo1D(hmodel, var)
            else:
                self.histos[chann] = rdf_mc.Filter(f"mctruth == {key}").Histo1D(hmodel, var)

        # MC sum — created manually after triggering lazy evaluation
        self.histos["MC sum"] = ROOT.TH1D(f"{prefix}_MC sum", title, nbins, xmin, xmax)

    def _build_mc_sum(self):
        mc_sum = self.histos["MC sum"]
        mc_sum.Reset()
        for chann, hist in self.histos.items():
            if chann != "Data" and chann != "MC sum":
                mc_sum.Add(hist.GetPtr() if hasattr(hist, 'GetPtr') else hist)

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

    def draw(self, canvas_name="c1", width=850, height=850,
             scale_mc=True, legend_pos=(0.65, 0.65, 0.88, 0.88), save_as=None):

        self._build_mc_sum()
        if scale_mc:
            self._scale_mc_to_data()

        canva = ROOT.TCanvas(canvas_name, "", width, height)

        # Find channel with highest peak for drawing first
        peaks = {}
        for chann, hist in self.histos.items():
            h = hist.GetPtr() if hasattr(hist, 'GetPtr') else hist
            peaks[chann] = h.GetMaximum()

        draw_order = sorted(peaks, key=peaks.get, reverse=True)

        legend = ROOT.TLegend(*legend_pos)
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
                first = False
            else:
                h.Draw("PE1SAME" if is_data else "HISTSAME")

            legend.AddEntry(h, chann, "lp" if is_data else "l")

        legend.Draw()

        if save_as:
            canva.SaveAs(save_as)

        return canva, legend
