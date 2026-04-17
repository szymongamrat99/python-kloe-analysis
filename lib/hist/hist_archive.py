import ROOT
import json
import os
import hashlib
from datetime import datetime


class HistArchive:
    """Archive histograms to dated ROOT files with JSON metadata.

    Structure:
        output/
          2026-04-17/
            run_abc123.root        # histograms
            run_abc123_meta.json   # full metadata (cuts, columns, models)

    Usage:
        archive = HistArchive("output/", columns={"col": "expr"}, global_filter="...")
        archive.add("chi2", mch, hist_model=hm)
        archive.add("mass", mch2, hist_model=hm2)
        archive.save()  # writes .root + .json
    """

    def __init__(self, base_dir, columns=None, global_filter=""):
        self.base_dir = base_dir
        self.global_filter = global_filter
        self.columns = columns or {}
        self._entries = []
        self._run_id = self._make_run_id()
        self._date_str = datetime.now().strftime("%Y-%m-%d")

    def _make_run_id(self):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        raw = f"{timestamp}_{id(self)}"
        short_hash = hashlib.md5(raw.encode()).hexdigest()[:8]
        return f"run_{timestamp}_{short_hash}"

    def add(self, name, multi_channel_hist, hist_model=None):
        """Register a histogram group for archiving."""
        self._entries.append({
            "name": name,
            "mch": multi_channel_hist,
            "hist_model": hist_model,
        })
        return self

    def save(self):
        """Write ROOT file with histograms and JSON with metadata."""
        out_dir = os.path.join(self.base_dir, self._date_str)
        os.makedirs(out_dir, exist_ok=True)

        root_path = os.path.join(out_dir, f"{self._run_id}.root")
        json_path = os.path.join(out_dir, f"{self._run_id}_meta.json")

        # Build full metadata
        metadata = {
            "run_id": self._run_id,
            "date": self._date_str,
            "timestamp": datetime.now().isoformat(),
            "global_filter": self.global_filter,
            "columns": self.columns,
            "histograms": {},
        }

        # Write ROOT file
        tfile = ROOT.TFile.Open(root_path, "RECREATE")
        if not tfile or tfile.IsZombie():
            raise IOError(f"Cannot create archive file: {root_path}")

        for entry in self._entries:
            name = entry["name"]
            mch = entry["mch"]
            hm = entry["hist_model"]

            # Ensure MC sum is built
            mch._build_mc_sum()

            # Create a subdirectory per histogram group
            tdir = tfile.mkdir(name)
            tdir.cd()

            hist_meta = {
                "per_hist_filter": hm.base_filter if hm else "",
                "combined_filter": self._combine_filters(hm),
                "variable": hm.var if hm else "",
            }

            if hm:
                hist_meta["model"] = {
                    "title": hm.title,
                    "nbins": hm.nbins,
                    "xmin": hm.xmin,
                    "xmax": hm.xmax,
                    "dim": hm.dim,
                }
                if hm.dim == 2:
                    hist_meta["model"].update({
                        "var_y": hm.var_y,
                        "nbins_y": hm.nbins_y,
                        "ymin": hm.ymin,
                        "ymax": hm.ymax,
                    })

            channels_info = {}
            for chann, hist in mch.histos.items():
                h = hist.GetPtr() if hasattr(hist, 'GetPtr') else hist
                clone = h.Clone(f"{name}_{chann}")
                clone.SetDirectory(tdir)
                clone.Write()
                channels_info[chann] = {
                    "entries": h.GetEntries(),
                    "mean": h.GetMean(),
                    "rms": h.GetRMS(),
                }

            hist_meta["channels"] = channels_info

            # Purity and efficiency (Signal channel)
            signal_after = channels_info.get("Signal", {}).get("entries", 0)
            bkg_after = sum(
                info["entries"] for ch, info in channels_info.items()
                if ch not in ("Signal", "Data", "MC sum")
            )

            total_after = signal_after + bkg_after
            hist_meta["purity"] = signal_after / total_after if total_after > 0 else 0.0

            # Efficiency: Signal_after / Signal_before
            counts_before = getattr(mch, '_counts_before', {})
            signal_before_count = counts_before.get("Signal", None)
            if signal_before_count is not None:
                sb = signal_before_count.GetValue() if hasattr(signal_before_count, 'GetValue') else signal_before_count
                hist_meta["efficiency"] = signal_after / sb if sb > 0 else 0.0
                hist_meta["signal_before"] = sb
            else:
                hist_meta["efficiency"] = None
                hist_meta["signal_before"] = None

            hist_meta["signal_after"] = signal_after
            hist_meta["background_after"] = bkg_after

            metadata["histograms"][name] = hist_meta

        tfile.Close()

        # Write JSON
        with open(json_path, 'w') as f:
            json.dump(metadata, f, indent=2)

        print(f"[archive] Saved: {root_path}")
        print(f"[archive] Meta:  {json_path}")

        return root_path, json_path

    def _combine_filters(self, hist_model):
        filters = []
        if self.global_filter:
            filters.append(self.global_filter)
        if hist_model and hist_model.base_filter:
            filters.append(hist_model.base_filter)
        return " && ".join(f"({f})" for f in filters) if filters else ""
