# Disable ROOT logon file loading — MUST be before import ROOT
import ROOT
ROOT.PyConfig.DisableRootLogon = True
ROOT.gROOT.SetBatch(True)

import os
import sys
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

path_to_log = os.path.join(os.path.dirname(__file__), "..", "log/")
logger = ROOT.ErrorHandling.ErrorLogs(str(path_to_log)) # Initialize the logger for error handling
ROOT.Utils.InitializeVariables(logger) # Initialize physics constants
pm00 = ROOT.KLOE.pm00() # General pm00 class

channName, channColor = get_chann_dicts() # Chann name and color dicts to python

physics_error_names = [
    "NO_VTX_WITH_TWO_TRACKS", "LESS_THAN_FOUR_NEUTRAL_CLUSTERS", "NO_VTX_WITH_OPPOSITE_TRACKS", "LESS_THAN_FOUR_CLUSTERS_WITH_GOOD_ENERGY", 
    "TRILATERATION_KIN_FIT"
]

ErrorCodes = ROOT.ErrorHandling.ErrorCodes  # adjust namespace
error_dict = {name: int(getattr(ErrorCodes, name)) for name in physics_error_names}
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
hist_models = HistModelLoader(os.path.join(os.path.dirname(__file__), "..", "config", "histograms.json"))

chi2Cut = 30.0
mPiChCut = 2.637
mCombPi0Cut = 10.
mPiNeCut = 76.
qMissCut = 3.75
bestErrorSixGammaCut = 1000.0

Kchrec_mass_mean = 497.605

# Cuts on pi0-plane
pi0Mass1_mean = 134.83240924168848
pi0Mass1_sigma = 3.402793221980809
pi0Mass2_mean = 134.87134080668446
pi0Mass2_sigma = 3.22758797811996
rho = -0.3399659203793453

u_cut = 8.085754304
v_cut = 11.514665004

u_column = f"((pi01Fit[5] - {pi0Mass1_mean}) + (pi02Fit[5] - {pi0Mass2_mean})) / sqrt(2)"
v_column = f"((pi01Fit[5] - {pi0Mass1_mean}) - (pi02Fit[5] - {pi0Mass2_mean})) / sqrt(2)"

# Cut for resolution improvement
radius_00_limit = 1.5
radius_pm_limit = 2.0
z_00_limit = 1.5
z_pm_limit = 1.5

fiducial_Volume_Close = f"radius_00 < {radius_00_limit} && z_00 < {z_00_limit} && radius_pm < {radius_pm_limit} && z_pm < {z_pm_limit}"

# Cut for omega-pi0 rejection
radius_00_pm_limit = 2.05
z_00_pm_limit = 2.45

fiducial_Volume_Omega = f"radius_00_pm < {radius_00_pm_limit} && z_00_pm < {z_00_pm_limit}"

# Angle between trk and Kchrec
angle_trk1_kchrec_limit = 0.8
angle_trk2_kchrec_limit = 0.8
# ---
#Rho pm
rho_limit = 0.8
# ---
# Definition of columns
rdf = rdf.Define("u", u_column).Define("v", v_column)

# Cut list
single_cuts = {
                f"Chi2SignalKinFit < {chi2Cut}": True,
                f"abs(Kchrec_mass_centered) < {mPiChCut}": True,
                f"abs(u) < {u_cut}": True,
                f"abs(v) < {v_cut}": True,
                f"(({fiducial_Volume_Close} && (abs(cosAngleTrk1Kch) < {angle_trk1_kchrec_limit} && abs(cosAngleTrk2Kch) < {angle_trk2_kchrec_limit})) || !({fiducial_Volume_Close}))": True,
                f"(({fiducial_Volume_Omega} && rho > {rho_limit}) || !({fiducial_Volume_Omega}))": True,
                f"abs(Knerec_mass_centered) < {mPiNeCut}": True,
                f"Qmiss < {qMissCut}": False,
                f"bestErrorSixGamma > {bestErrorSixGammaCut}": False
              }

# Build list of enabled cuts (preserving order)
enabled_cuts = [cut for cut, enabled in single_cuts.items() if enabled]


def build_global_filter(cuts):
    """Build a global filter string from a list of cut expressions."""
    if not cuts:
        return ""
    return "(" + " && ".join(f"({c})" for c in cuts) + ")"


def run_analysis(rdf, global_filter, img_dir, output_suffix=""):
    """Run the full analysis (efficiency, purity, histograms) with a given filter."""

    suffix_dir = os.path.join(img_dir, output_suffix) if output_suffix else img_dir
    os.makedirs(suffix_dir, exist_ok=True)

    # Save cut info
    if output_suffix:
        with open(os.path.join(suffix_dir, "cuts_applied.txt"), "w") as f:
            f.write(f"Global filter:\n{global_filter}\n")

    scale_map = {
      "Signal": 0.996761,
      "Omega": 1.07463,
      "Regeneration": [0.499798, 3.47319, 3.71238, 0.575271],
      "3pi0": 0.782048,
      "Semileptonic": 0.0,
      "Other": 3.488
    }

    # scale_map = {
    #   "Signal": 1.0,
    #   "Omega": 1.0,
    #   "Regeneration": [1.0, 1.0, 1.0, 1.0],
    #   "3pi0": 1.0,
    #   "Semileptonic": 1.0,
    #   "Other": 1.0
    # }

    # Efficiency
    eff_pur_calculator = EfficiencyPurityCalculator(rdf, truth_condition="mctruth == 1", truth_condition_before_cut="mctruth == 0", bin_width=1.0, xmin=-50.0, xmax=50.0)

    total, selected, eff_errors, all_good_clusters, single_bad_cluster = eff_pur_calculator.calculate_preselection_eff(err_codes=error_dict)

    print(f"Total true signal events: {total}")
    for err_name, count in selected.items():
        print(f"Selected events with error '{err_name}': {count} (efficiency: {eff_errors[err_name]:.4f})")
        print(f"  - All good clusters: {all_good_clusters[err_name]}")
        print(f"  - Single bad cluster: {single_bad_cluster[err_name]}")
        print()

    total_true, selected_true, true_selected, selected_count, efficiency, purity = eff_pur_calculator.calculate_efficiency_purity(selection=global_filter)

    print(f"Efficiency: {efficiency:.4f} ({selected_true}/{total_true})")
    print(f"Purity: {purity:.4f} ({true_selected}/{selected_count})")

    efficiency_hist = eff_pur_calculator.efficiency_histo_per_time_difference(selection=global_filter)

    canvas_eff = ROOT.TCanvas("canvas_efficiency", "Efficiency vs #Deltat", 800, 800)
    efficiency_hist.SetLineColor(ROOT.kBlack)
    efficiency_hist.Draw("PE1")
    ROOT.gPad.Update()
    efficiency_hist.GetPaintedGraph().GetYaxis().SetRangeUser(0.0, 1.0)
    canvas_eff.SaveAs(os.path.join(suffix_dir, "efficiency_vs_deltat.svg"))

    histSelectedSignal, histSelected, histPurity = eff_pur_calculator.purity_histo_per_time_difference(selection=global_filter)
    canvas_pur = ROOT.TCanvas("canvas_purity", "Purity vs #Deltat", 800, 800)
    histPurity.SetLineColor(ROOT.kBlack)
    histPurity.Draw("PE1")
    ROOT.gPad.Update()
    histPurity.GetPaintedGraph().GetYaxis().SetRangeUser(0.0, 1.0)
    canvas_pur.SaveAs(os.path.join(suffix_dir, "purity_vs_deltat.svg"))

    canvas_selsig = ROOT.TCanvas("canvas_selected_signal", "Selected Signal vs #Deltat", 800, 800)
    histSelectedSignal.SetLineColor(ROOT.kOrange)
    histSelectedSignal.Draw("HIST")
    ROOT.gPad.Update()
    histSelectedSignal.GetYaxis().SetRangeUser(0.0, 1.2 * histSelectedSignal.GetMaximum())
    canvas_selsig.SaveAs(os.path.join(suffix_dir, "selected_signal_vs_deltat.svg"))

    canvas_sel = ROOT.TCanvas("canvas_selected", "Selected vs #Deltat", 800, 800)
    histSelected.SetLineColor(ROOT.kOrange)
    histSelected.Draw("HIST")
    ROOT.gPad.Update()
    histSelected.GetYaxis().SetRangeUser(0.0, 1.2 * histSelected.GetMaximum())
    canvas_sel.SaveAs(os.path.join(suffix_dir, "selected_vs_deltat.svg"))

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
            print(f"[2D] {name} — computing...")
            dir_name = os.path.join(suffix_dir, name)
            os.makedirs(dir_name, exist_ok=True)
            mch2d = MultiChannelHist2D(rdf, channName, channColor,
                                       hist_model=hm, global_filter=global_filter)
            mch2d.build()
            mch2d.draw(save_as=os.path.join(dir_name, f"{name}.svg"))
            continue

        # --- 1D histograms ---
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
                                   hist_model=hm, global_filter=global_filter, scale_map=scale_map)
            mch.build()
            cache.save(name, mch, global_filter=global_filter, hist_model=hm)

        archive.add(name, mch, hist_model=hm)
        mch.draw(save_as=os.path.join(suffix_dir, f"{name}.svg"), scale_mc=False)

    archive.save()


# --- Mode selection ---
mode = int(sys.argv[1]) if len(sys.argv) > 1 else 2

if mode == 1:
    # Sequential cuts: run analysis adding one cut at a time
    print("=== MODE 1: Sequential cuts ===\n")
    for i in range(1, len(enabled_cuts) + 1):
        cuts_subset = enabled_cuts[:i]
        gf = build_global_filter(cuts_subset)
        folder_name = f"step_{i:02d}_cuts_{i}"
        print(f"\n{'='*60}")
        print(f"Step {i}/{len(enabled_cuts)}: applying cuts 1..{i}")
        print(f"Filter: {gf}")
        print(f"{'='*60}\n")
        run_analysis(rdf, gf, img_dir, output_suffix=folder_name)

elif mode == 2:
    # Combined filter: all enabled cuts at once (original behavior)
    print("=== MODE 2: Combined global filter ===\n")
    global_filter = build_global_filter(enabled_cuts)
    print(f"Filter: {global_filter}\n")
    run_analysis(rdf, global_filter, img_dir)

else:
    print(f"Unknown mode: {mode}. Use 1 (sequential) or 2 (combined).")
    sys.exit(1)