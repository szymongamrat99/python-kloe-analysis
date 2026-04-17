import ROOT
import json
import os
import hashlib


class HistCache:
    """Cache MultiChannelHist results in a ROOT file with metadata.

    Metadata stored per histogram:
        - filter (global + per-histogram)
        - column definitions used
        - histogram model config

    Usage:
        cache = HistCache("output/cache.root", columns_json="config/columns.json")
        # If cache is valid, loads from file; otherwise returns None
        histos = cache.load("chi2", global_filter="...", hist_model=hm)
        if histos is None:
            mch = MultiChannelHist(rdf, ...)
            cache.save("chi2", mch, global_filter="...", hist_model=hm)
    """

    def __init__(self, cache_path, columns_json=None):
        self.cache_path = cache_path
        self._columns = {}
        if columns_json and os.path.exists(columns_json):
            with open(columns_json, 'r') as f:
                self._columns = json.load(f)

    def _metadata_key(self, name):
        return f"{name}__metadata"

    def _build_metadata(self, name, global_filter, hist_model):
        meta = {
            "name": name,
            "global_filter": global_filter or "",
            "hist_filter": hist_model.base_filter if hist_model else "",
            "var": hist_model.var if hist_model else "",
            "columns": self._columns,
        }
        if hist_model:
            meta["model"] = {
                "title": hist_model.title,
                "nbins": hist_model.nbins,
                "xmin": hist_model.xmin,
                "xmax": hist_model.xmax,
                "dim": hist_model.dim,
            }
            if hist_model.dim == 2:
                meta["model"].update({
                    "var_y": hist_model.var_y,
                    "nbins_y": hist_model.nbins_y,
                    "ymin": hist_model.ymin,
                    "ymax": hist_model.ymax,
                })
        return meta

    def _meta_hash(self, meta):
        raw = json.dumps(meta, sort_keys=True)
        return hashlib.sha256(raw.encode()).hexdigest()

    def _read_stored_hash(self, tfile, name):
        key = self._metadata_key(name)
        obj = tfile.Get(key)
        if not obj:
            return None
        return str(obj.GetTitle())

    def load(self, name, global_filter="", hist_model=None, chann_names=None):
        """Try to load cached histograms. Returns dict {chann: TH1D} or None if stale/missing."""
        if not os.path.exists(self.cache_path):
            return None

        meta = self._build_metadata(name, global_filter, hist_model)
        current_hash = self._meta_hash(meta)

        tfile = ROOT.TFile.Open(self.cache_path, "READ")
        if not tfile or tfile.IsZombie():
            return None

        stored_hash = self._read_stored_hash(tfile, name)
        if stored_hash != current_hash:
            tfile.Close()
            return None

        histos = {}
        if chann_names is None:
            chann_names = []

        for chann in chann_names:
            hist_name = f"{name}_{chann}"
            h = tfile.Get(hist_name)
            if not h:
                tfile.Close()
                return None
            h.SetDirectory(0)
            histos[chann] = h

        tfile.Close()
        return histos

    def save(self, name, multi_channel_hist, global_filter="", hist_model=None):
        """Save all histograms from a MultiChannelHist + metadata to cache file."""
        meta = self._build_metadata(name, global_filter, hist_model)
        current_hash = self._meta_hash(meta)
        meta_json = json.dumps(meta, sort_keys=True, indent=2)

        cache_dir = os.path.dirname(self.cache_path)
        if cache_dir and not os.path.exists(cache_dir):
            os.makedirs(cache_dir)

        # Open in UPDATE mode to preserve other histograms
        if os.path.exists(self.cache_path):
            tfile = ROOT.TFile.Open(self.cache_path, "UPDATE")
        else:
            tfile = ROOT.TFile.Open(self.cache_path, "RECREATE")

        if not tfile or tfile.IsZombie():
            raise IOError(f"Cannot open cache file: {self.cache_path}")

        for chann, hist in multi_channel_hist.histos.items():
            h = hist.GetPtr() if hasattr(hist, 'GetPtr') else hist
            clone = h.Clone(f"{name}_{chann}")
            clone.SetDirectory(tfile)
            clone.Write("", ROOT.TObject.kOverwrite)

        # Store hash as TNamed
        hash_obj = ROOT.TNamed(self._metadata_key(name), current_hash)
        hash_obj.Write("", ROOT.TObject.kOverwrite)

        # Store full metadata as TNamed
        meta_obj = ROOT.TNamed(f"{name}__metadata_json", meta_json)
        meta_obj.Write("", ROOT.TObject.kOverwrite)

        tfile.Close()

    def list_cached(self):
        """List all histogram groups stored in the cache."""
        if not os.path.exists(self.cache_path):
            return []

        tfile = ROOT.TFile.Open(self.cache_path, "READ")
        if not tfile or tfile.IsZombie():
            return []

        names = set()
        for key in tfile.GetListOfKeys():
            kname = key.GetName()
            if kname.endswith("__metadata_json"):
                names.add(kname.replace("__metadata_json", ""))

        tfile.Close()
        return sorted(names)

    def get_metadata(self, name):
        """Read stored metadata JSON for a histogram group."""
        if not os.path.exists(self.cache_path):
            return None

        tfile = ROOT.TFile.Open(self.cache_path, "READ")
        if not tfile or tfile.IsZombie():
            return None

        obj = tfile.Get(f"{name}__metadata_json")
        if not obj:
            tfile.Close()
            return None

        meta = json.loads(str(obj.GetTitle()))
        tfile.Close()
        return meta
