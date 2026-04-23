import ROOT
import json


class HistModel:
    """Wrapper around a loaded histogram model definition."""

    def __init__(self, name, config):
        self.name = name
        self.var = config["var"]
        self.title = config.get("title", "")
        self.nbins = config["nbins"]
        self.xmin = config["xmin"]
        self.xmax = config["xmax"]
        self.dim = config.get("dim", 1)
        self.base_filter = config.get("filter", "")

        # Display options
        self.ymin = config.get("ymin_display", None)
        self.ymax_display = config.get("ymax_display", None)
        self.log_x = config.get("log_x", False)
        self.log_y = config.get("log_y", False)
        self.legend_pos = config.get("legend_pos", [0.65, 0.65, 0.88, 0.88])
        self.show_stats = config.get("show_stats", False)

        # Weighting
        self.weight = config.get("weight", None)
        
        # Is Resolution plot?
        self.resolution = config.get("resolution", False)

        # 2D-specific
        if self.dim == 2:
            self.var_y = config["var_y"]
            self.nbins_y = config["nbins_y"]
            self.ymin = config["ymin"]
            self.ymax = config["ymax"]
            self.log_z = config.get("log_z", False)

    def make_root_model(self, suffix=""):
        """Return a tuple that ROOT's Pythonization converts to TH1DModel/TH2DModel."""
        full_name = f"{self.name}_{suffix}" if suffix else self.name

        if self.dim == 2:
            return (
                full_name, self.title,
                self.nbins, self.xmin, self.xmax,
                self.nbins_y, self.ymin, self.ymax,
            )
        else:
            return (
                full_name, self.title,
                self.nbins, self.xmin, self.xmax,
            )

    def make_empty_th(self, suffix=""):
        """Return an empty TH1D or TH2D (for MC sum etc.)."""
        full_name = f"{self.name}_{suffix}" if suffix else self.name

        if self.dim == 2:
            return ROOT.TH2D(
                full_name, self.title,
                self.nbins, self.xmin, self.xmax,
                self.nbins_y, self.ymin, self.ymax,
            )
        else:
            return ROOT.TH1D(
                full_name, self.title,
                self.nbins, self.xmin, self.xmax,
            )

    def __repr__(self):
        if self.dim == 2:
            return f"HistModel({self.name}, {self.var} vs {self.var_y}, {self.dim}D)"
        return f"HistModel({self.name}, {self.var}, {self.dim}D)"


class HistModelLoader:
    """Load histogram model definitions from a JSON file.

    JSON format:
        {
            "chi2": {
                "var": "Chi2SignalKinFit",
                "title": "Chi2; #chi^{2}; Entries",
                "nbins": 100, "xmin": 0, "xmax": 100,
                "filter": "abs(minv4gam - 497.6) < 76"
            },
            "mass_2d": {
                "dim": 2,
                "var": "pi01Fit[5]", "var_y": "pi02Fit[5]",
                "title": "Pi0 masses; m1 [MeV]; m2 [MeV]",
                "nbins": 100, "xmin": 50, "xmax": 250,
                "nbins_y": 100, "ymin": 50, "ymax": 250
            }
        }
    """

    def __init__(self, path):
        self.path = path
        self._models = {}
        self._load(path)

    def _load(self, path):
        with open(path, 'r') as f:
            data = json.load(f)

        for name, cfg in data.items():
            self._models[name] = HistModel(name, cfg)

    def get(self, name):
        """Return a single HistModel by name."""
        if name not in self._models:
            raise KeyError(f"Histogram model '{name}' not found. Available: {list(self._models.keys())}")
        return self._models[name]

    def all(self):
        """Return dict of all loaded models."""
        return dict(self._models)

    def names(self):
        return list(self._models.keys())

    def __getitem__(self, name):
        return self.get(name)

    def __repr__(self):
        return f"HistModelLoader({self.names()})"
