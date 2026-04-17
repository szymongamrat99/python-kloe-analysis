import ROOT
import os
import math as mt
from lib.utils.initial_declarations import InitialDeclarations
from lib.utils.utils import load_files_into_chain, get_chann_dicts
from lib.hist.multi_channel_hist import MultiChannelHist
from lib.hist.hist_model_loader import HistModelLoader
from lib.hist.hist_cache import HistCache
from lib.hist.hist_archive import HistArchive
from lib.rdf.column_definer import ColumnDefiner

# Disable ROOT logon file loading
ROOT.PyConfig.DisableRootLogon = True
ROOT.gROOT.SetBatch(True)

# Load of the initial declarations from the config file
#-----------------------------------------
path_to_config = os.path.join(os.path.dirname(__file__), "..", "config", "config.json")
try:
    InitialDeclarations(path_to_config)
except Exception as e:
    print(f"Error loading initial declarations: {e}")
#-----------------------------------------
# Initial root config
ROOT.gROOT.SetBatch(True)

ROOT.KLOE.setGlobalStyle() # Set the global style for KLOE analysis

logger = ROOT.ErrorHandling.ErrorLogs("../log/") # Initialize the logger for error handling
ROOT.Utils.InitializeVariables(logger) # Initialize physics constants
pm00 = ROOT.KLOE.pm00() # General pm00 class

channName, channColor = get_chann_dicts() # Chann name and color dicts to python
# ----------------------------------------
# Set up the TChain and load files into it
chain = ROOT.TChain("h1")
load_files_into_chain(config_path=path_to_config, chain=chain) 
# ----------------------------------------
# Path to image folder
img_dir = os.path.join(os.path.dirname(__file__), "..", "img")
if not os.path.exists(img_dir):
    os.makedirs(img_dir)
# ----------------------------------------

rdf = ROOT.RDataFrame(chain)

# Declare C++ helper wrappers for functions with pointer-based signatures
ROOT.gInterpreter.Declare("""
double wrap_double_exponential(double dt) {
    double x[1] = {dt};
    double par[1] = {0};
    return double_exponential(x, par);
};

double wrap_interf_function(double dt) {
    double x[1] = {dt};
    double par[2] = {PhysicsConstants::Re, PhysicsConstants::Im_nonCPT};
    return interf_function(x, par);
};
""")

columns_config = os.path.join(os.path.dirname(__file__), "..", "config", "columns.json")
rdf = ColumnDefiner(rdf).from_json(columns_config).apply()

# Load histogram models from config
hist_models = HistModelLoader(os.path.join(os.path.dirname(__file__), "..", "config", "histograms.json"))

# Global filter applied to all histograms
global_filter = f"abs(minv4gam - {ROOT.PhysicsConstants.mK0}) < 76"

# Cache file
cache_path = os.path.join(os.path.dirname(__file__), "..", "output", "hist_cache.root")
cache = HistCache(cache_path, columns_json=columns_config)

chann_list = list(channName.values())

# Archive for dated ROOT files + JSON metadata
output_dir = os.path.join(os.path.dirname(__file__), "..", "output")
archive = HistArchive(output_dir, columns=ColumnDefiner._load_json(columns_config),
                      global_filter=global_filter)

for name, hm in hist_models.all().items():
    if hm.dim != 1:
        continue

    # Try loading from cache
    cached = cache.load(name, global_filter=global_filter,
                        hist_model=hm, chann_names=chann_list)

    if cached is not None:
        print(f"[cache hit] {name}")
        mch = MultiChannelHist.__new__(MultiChannelHist)
        mch.channName = channName
        mch.channColor = channColor
        mch.histos = cached
        mch._hist_model = hm
        mch._mc_sum_built = True
    else:
        print(f"[cache miss] {name} — computing...")
        mch = MultiChannelHist(rdf, channName, channColor,
                               hist_model=hm, global_filter=global_filter)
        mch.build()  # trigger RDF evaluation + MC sum
        cache.save(name, mch, global_filter=global_filter, hist_model=hm)

    archive.add(name, mch, hist_model=hm)
    mch.draw(save_as=os.path.join(img_dir, f"{name}.svg"))

archive.save()