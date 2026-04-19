import ROOT
from lib.hist.hist_model_loader import HistModel


class EfficiencyPurityCalculator:
  def __init__(self, rdf: ROOT.RDataFrame, truth_condition: str, truth_condition_before_cut: str, bin_width=1.0, xmin=-100, xmax=100):
    self.time_diff_MC = "KaonChTimeCMMC - KaonNeTimeCMMC"
    self.time_diff_rec = "KaonChTimeCMSignalFit - KaonNeTimeCMSignalFit"

    define = "Redefine" if rdf.HasColumn("deltaT_MC") else "Define"
    self.rdf = getattr(rdf, define)("deltaT_MC", self.time_diff_MC)
    define = "Redefine" if self.rdf.HasColumn("deltaT") else "Define"
    self.rdf = getattr(self.rdf, define)("deltaT", self.time_diff_rec)

    mc_total = f"mcflag == 1 && (mctruth == -1 || {truth_condition} || {truth_condition_before_cut})"
    self.rdfSignalTot = self.rdf.Filter(mc_total)
    self.rdfMC = self.rdf.Filter("mcflag == 1 && mctruth != -1")
    self.rdfData = self.rdf.Filter("mcflag == 0")

    self.xmin = xmin
    self.xmax = xmax
    self.nbins = int((self.xmax - self.xmin) / bin_width)

    if self.nbins % 2 == 0:
        self.nbins += 1  # Ensure odd number of bins for symmetry around zero

    self.truth_condition = f"({truth_condition})" if truth_condition else None
    self.truth_condition_before_cut = f"({truth_condition_before_cut})" if truth_condition_before_cut else None

    self.combined_truth_condition = f"({self.truth_condition} || {self.truth_condition_before_cut})" if self.truth_condition_before_cut else self.truth_condition

  def efficiency_histo_per_time_difference(self, selection: str):

    histSignalTotal = self.rdfMC.Filter(self.combined_truth_condition).Histo1D(("selected_signal", ";#Deltat [#tau_{S}];Efficiency [-]", self.nbins, self.xmin, self.xmax), "deltaT_MC", "interference_weight")

    scale = histSignalTotal.GetEntries() / histSignalTotal.Integral(0, self.nbins + 1) if histSignalTotal.GetEntries() > 0 else 1.0
    histSignalTotal.Scale(scale)

    selection = selection + " && " + self.combined_truth_condition if selection else self.combined_truth_condition
    histSelected = self.rdfMC.Filter(selection).Histo1D(("selected_signal", ";#Deltat [#tau_{S}];Efficiency [-]", self.nbins, self.xmin, self.xmax), "deltaT_MC", "interference_weight")

    scale = histSelected.GetEntries() / histSelected.Integral(0, self.nbins + 1) if histSelected.GetEntries() > 0 else 1.0
    histSelected.Scale(scale)

    self.histEfficiency = ROOT.TEfficiency(histSelected.GetPtr(), histSignalTotal.GetPtr())
    self.histEfficiency.SetStatisticOption(ROOT.TEfficiency.kBUniform)

    self.histEfficiency.SetTitle(";#Deltat^{MC} [#tau_{S}];Efficiency [-]")
    self.histEfficiency.SetName("efficiency_per_time_difference")

    return self.histEfficiency

  def purity_histo_per_time_difference(self, selection: str):

    comb_selection = selection + " && " + self.truth_condition if selection else self.truth_condition

    bkg_selection = selection + " && " + "!" + self.combined_truth_condition if selection else "!" + self.combined_truth_condition 


    histSelectedSignal = self.rdfMC.Filter(comb_selection).Histo1D(("selected_signal_only", ";#Deltat [#tau_{S}];Purity [-]", self.nbins, self.xmin, self.xmax), "deltaT", "interference_weight")

    scale = histSelectedSignal.GetEntries() / histSelectedSignal.Integral(0, self.nbins + 1) if histSelectedSignal.GetEntries() > 0 else 1.0
    histSelectedSignal.Scale(scale)

    histSelected = self.rdfMC.Filter(bkg_selection).Histo1D(("selected_signal_bkg", ";#Deltat [#tau_{S}];Purity [-]", self.nbins, self.xmin, self.xmax), "deltaT")

    histSelected.Add(histSelectedSignal.GetPtr())

    self.histPurity = ROOT.TEfficiency(histSelectedSignal.GetPtr(), histSelected.GetPtr())
    self.histPurity.SetStatisticOption(ROOT.TEfficiency.kBUniform)

    self.histPurity.SetTitle(";#Deltat [#tau_{S}];Purity [-]")
    self.histPurity.SetName("purity_per_time_difference")

    return histSelectedSignal, histSelected, self.histPurity

  def calculate_efficiency_purity(self, selection: str) -> tuple:
    """
    Calculate efficiency and purity for a given selection and truth condition.
    Calculates only for MC events.

    Args:
        selection (str): The selection criteria to apply to the data.
        truth_condition (str): The condition that defines the "true" events.

    Returns:
        tuple: A tuple containing the efficiency and purity values.
    """
    
    comb_total_selection = selection + " && " + self.combined_truth_condition if selection else self.truth_condition

    comb_selection = selection + " && " + self.truth_condition if selection else self.truth_condition

    bkg_selection = selection + " && " + "!" + self.truth_condition_before_cut if selection else "!" + self.truth_condition_before_cut 

    # Count total true events
    total_true = self.rdfMC.Filter(self.combined_truth_condition).Count().GetValue()

    # Count selected events
    selected_true = self.rdfMC.Filter(comb_total_selection).Count().GetValue()

    # Count true events that are selected
    true_selected = self.rdfMC.Filter(comb_selection).Count().GetValue()

    selected = self.rdfMC.Filter(bkg_selection).Count().GetValue()

    # Calculate efficiency and purity
    efficiency = selected_true / total_true if total_true > 0 else 0
    purity = true_selected / selected if selected > 0 else 0

    return total_true, selected_true, true_selected, selected, efficiency, purity

  def calculate_preselection_eff(self, err_codes: dict[str, int]):

    total_true = self.rdfSignalTot.Count().GetValue()

    selection = ""

    all_good_clusters = {}
    single_bad_cluster = {}

    selected = {}
    efficiency = {}

    for name, code in err_codes.items():
      selection = f"{selection} && (errorcode != {code})" if selection else f"(errorcode != {code})"

      selected_true = self.rdfSignalTot.Filter(selection).Count().GetValue()

      efficiency[name] = selected_true / total_true if total_true > 0 else 0
      selected[name] = selected_true

      all_good_clusters[name] = self.rdfSignalTot.Filter(selection + f" && (goodClustersTriKinFitSize == 4)").Count().GetValue()
      single_bad_cluster[name] = self.rdfSignalTot.Filter(selection + f" && (goodClustersTriKinFitSize >= 3)").Count().GetValue()


    return total_true, selected, efficiency, all_good_clusters, single_bad_cluster