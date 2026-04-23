import ROOT
from lib.hist.hist_model_loader import HistModel


class MultiChannelHist2D:
    """Per-channel 2D histograms (one per MC channel + Data), drawn with COLZ."""

    def __init__(self, rdf, channName, channColor, hist_model=None,
                 global_filter=""):
        self.channName = channName
        self.channColor = channColor
        self.histos = {}
        self._hist_model = hist_model

        self.global_filter = global_filter

        # Combine global + per-histogram filters
        filters = [f for f in [global_filter, hist_model.base_filter] if f]
        combined_filter = " && ".join(f"({f})" for f in filters) if filters else ""

        rdf_mc = rdf.Filter("mcflag == 1 && mctruth != 0 && mctruth != -1")
        rdf_data = rdf.Filter("mcflag == 0")

        if combined_filter:
            rdf_mc = rdf_mc.Filter(combined_filter)
            rdf_data = rdf_data.Filter(combined_filter)

        var_x = hist_model.var
        var_y = hist_model.var_y
        weight = hist_model.weight

        for key, chann in channName.items():
            if chann == "MC sum":
                continue

            hmodel = hist_model.make_root_model(suffix=chann)

            if chann == "Data":
                    self.histos[chann] = rdf_data.Histo2D(hmodel, var_x, var_y)
            else:
                rdf_ch = rdf_mc.Filter(f"mctruth == {key}")
                if weight and key == 1:
                    self.histos[chann] = rdf_ch.Histo2D(hmodel, var_x, var_y, weight)
                else:
                    self.histos[chann] = rdf_ch.Histo2D(hmodel, var_x, var_y)

        self.histos["Signal"].Scale(self.histos["Signal"].GetEntries() / self.histos["Signal"].Integral(0, self.histos["Signal"].GetNbinsX() + 1, 0, self.histos["Signal"].GetNbinsY() + 1))

        self._built = False

    def build(self):
        """Trigger lazy RDF evaluation."""
        for chann, h in self.histos.items():
            if hasattr(h, 'GetPtr'):
                h.GetPtr()  # force evaluation
        self._built = True
        return self

    def draw(self, canvas_name=None, width=850, height=850, save_as=None):
        """Draw one canvas per channel, always COLZ."""
        hm = self._hist_model
        canvases = []

        for chann, hist in self.histos.items():
            h = hist.GetPtr() if hasattr(hist, 'GetPtr') else hist

            cname = f"c_{hm.name}_{chann}".replace(" ", "_")
            if canvas_name:
                cname = f"{canvas_name}_{chann}".replace(" ", "_")

            c = ROOT.TCanvas(cname, f"{hm.name} — {chann}", width, height)
            ROOT.gStyle.SetOptStat(0)

            if hm.log_z:
                c.SetLogz(1)

            h.Draw("COLZ")

            if save_as:
                # Insert channel name before extension
                base, ext = save_as.rsplit(".", 1) if "." in save_as else (save_as, "svg")
                chann_safe = chann.replace(" ", "_")
                path = f"{base}_{chann_safe}.{ext}"
                c.SaveAs(path)

            canvases.append(c)

        return canvases
