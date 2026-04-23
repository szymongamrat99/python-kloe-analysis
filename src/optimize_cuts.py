# Disable ROOT logon file loading — MUST be before import ROOT
import ROOT
ROOT.PyConfig.DisableRootLogon = True

import os
import sys
import math as mt
from lib.utils.initial_declarations import InitialDeclarations
from lib.rdf.column_definer import ColumnDefiner
from lib.utils.utils import load_files_into_chain, get_chann_dicts

from lib.hist.hist_model_loader import HistModelLoader
from lib.hist.hist_model_loader import HistModel

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

path_to_log = os.path.join(os.path.dirname(__file__), "..", "log/")
logger = ROOT.ErrorHandling.ErrorLogs(str(path_to_log)) # Initialize the logger for error handling
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

ROOT.gInterpreter.Declare("""
// Compute angle (in degrees) between two 3-vectors given as array components
double vec3Angle(double x1, double y1, double z1,
                 double x2, double y2, double z2) {
    TVector3 v1(x1, y1, z1);
    TVector3 v2(x2, y2, z2);
    return v1.Angle(v2);
}
""")

columns_config = os.path.join(os.path.dirname(__file__), "..", "config", "columns.json")
rdf = ColumnDefiner(rdf).from_json(columns_config).apply()

# Load histogram models from config
hist_models = HistModelLoader(os.path.join(os.path.dirname(__file__), "..", "config", "cut_optimization.json"))



histograms1D = {}
histograms2D = {}

for name, model in hist_models.items():
    if model.dim == 2:
      histograms2D[name] = rdf.Histo2D(model)
    else:
      histograms1D[name] = rdf.Histo1D(model)