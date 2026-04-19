# Disable ROOT logon file loading — MUST be before import ROOT
import ROOT
ROOT.PyConfig.DisableRootLogon = True
ROOT.gROOT.SetBatch(True)

import os
import math as mt
from lib.utils.initial_declarations import InitialDeclarations
from lib.utils.utils import load_files_into_chain, get_chann_dicts
from lib.hist.multi_channel_hist import MultiChannelHist
from lib.hist.multi_channel_hist_2d import MultiChannelHist2D
from lib.hist.hist_model_loader import HistModelLoader
from lib.hist.hist_cache import HistCache
from lib.hist.hist_archive import HistArchive
from lib.rdf.column_definer import ColumnDefiner
from lib.hist.efficiency_purity_calculator import EfficiencyPurityCalculator

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

chi2Cut = 30.0
mPiChCut = 1.2
mCombPi0Cut = 10.
mPiNeCut = 76.
qMissCut = 3.75

single_cuts = {
                f"Chi2SignalKinFit < {chi2Cut}": True,
                f"abs(Kchrec[5] - {ROOT.PhysicsConstants.mK0}) < {mPiChCut}": True,
                f"combined_Pi0 < {mCombPi0Cut}": True,
                f"abs(minv4gam - {ROOT.PhysicsConstants.mK0}) < {mPiNeCut}": True,
                f"Qmiss < {qMissCut}": True,
              }

global_filter = ""

for cut, enabled in single_cuts.items():
    if enabled:
        if global_filter:
            global_filter += " && "
        global_filter += f"({cut})"

global_filter = "(" + global_filter + ")" if global_filter else ""


# Efficiency
eff_pur_calculator = EfficiencyPurityCalculator(rdf, truth_condition="mctruth == 1", truth_condition_before_cut="mctruth == 0", bin_width=1.0, xmin=-30.0, xmax=30.0)
efficiency_hist = eff_pur_calculator.efficiency_histo_per_time_difference(selection=global_filter)

canvas_eff = ROOT.TCanvas("canvas_efficiency", "Efficiency vs #Deltat", 800, 800)
efficiency_hist.SetLineColor(ROOT.kBlack)
efficiency_hist.Draw("PE1")
ROOT.gPad.Update()
efficiency_hist.GetPaintedGraph().GetYaxis().SetRangeUser(0.0, 1.0)
canvas_eff.SaveAs(os.path.join(img_dir, "efficiency_vs_deltat.svg"))

histSelectedSignal, histSelected, histPurity = eff_pur_calculator.purity_histo_per_time_difference(selection=global_filter)
canvas_pur = ROOT.TCanvas("canvas_purity", "Purity vs #Deltat", 800, 800)
histPurity.SetLineColor(ROOT.kBlack)
histPurity.Draw("PE1")
ROOT.gPad.Update()
histPurity.GetPaintedGraph().GetYaxis().SetRangeUser(0.0, 1.0)
canvas_pur.SaveAs(os.path.join(img_dir, "purity_vs_deltat.svg"))

canvas_selsig = ROOT.TCanvas("canvas_selected_signal", "Selected Signal vs #Deltat", 800, 800)
histSelectedSignal.SetLineColor(ROOT.kOrange)
histSelectedSignal.Draw("HIST")
ROOT.gPad.Update()
histSelectedSignal.GetYaxis().SetRangeUser(0.0, 1.2 * histSelectedSignal.GetMaximum())
canvas_selsig.SaveAs(os.path.join(img_dir, "selected_signal_vs_deltat.svg"))

canvas_sel = ROOT.TCanvas("canvas_selected", "Selected vs #Deltat", 800, 800)
histSelected.SetLineColor(ROOT.kOrange)
histSelected.Draw("HIST")
ROOT.gPad.Update()
histSelected.GetYaxis().SetRangeUser(0.0, 1.2 * histSelected.GetMaximum())
canvas_sel.SaveAs(os.path.join(img_dir, "selected_vs_deltat.svg"))

# Cache file
cache_path = os.path.join(os.path.dirname(__file__), "..", "output", "hist_cache.root")
cache = HistCache(cache_path, columns_json=columns_config)

chann_list = list(channName.values())

# Archive for dated ROOT files + JSON metadata
output_dir = os.path.join(os.path.dirname(__file__), "..", "output")
archive = HistArchive(output_dir, columns=ColumnDefiner._load_json(columns_config),
                      global_filter=global_filter)

for name, hm in hist_models.all().items():
    if hm.dim > 2:
        continue

    if hm.dim == 2:
        # --- 2D histograms: one canvas per channel, COLZ ---
        print(f"[2D] {name} — computing...")
        dir_name = os.path.join(img_dir, name)
        os.makedirs(dir_name, exist_ok=True)
        mch2d = MultiChannelHist2D(rdf, channName, channColor,
                                   hist_model=hm, global_filter=global_filter)
        mch2d.build()
        mch2d.draw(save_as=os.path.join(dir_name, f"{name}.svg"))
        continue

    # --- 1D histograms ---
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