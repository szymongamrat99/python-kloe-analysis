import ROOT
import os
import math as mt
from lib.utils.initial_declarations import InitialDeclarations
from lib.utils.utils import load_files_into_chain, get_chann_dicts
from lib.hist.multi_channel_hist import MultiChannelHist
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
ROOT.EnableImplicitMT()

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

cols = ColumnDefiner(rdf)
cols.add("combined_Pi0",
         "sqrt(pow(pi01Fit[5] - PhysicsConstants::mPi0, 2)"
         " + pow(pi02Fit[5] - PhysicsConstants::mPi0, 2))")
cols.add("totalE", "pi01Fit[3] + pi02Fit[3]")

rdf = cols.apply() # Apply the column definitions to get a new RDF node with the new columns defined

chi2_hist = MultiChannelHist(
    rdf, channName, channColor,
    var="combined_Pi0",
    model=("hChi2", "Chi2 of the signal kinfit; #chi^{2} [-]; Entries", 100, 0, 10),
    base_filter=f"abs(minv4gam - {ROOT.PhysicsConstants.mK0}) < 76"
)

img_path = os.path.join(img_dir, "Chi2SignalKinFit.svg")
chi2_hist.draw(save_as=img_path)