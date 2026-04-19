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

    self.histEfficiency.SetTitle(";#Deltat [#tau_{S}];Efficiency [-]")
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
    
    rdfMC = self.rdf.Filter("mcflag == 1")

    # Count total true events
    total_true = rdfMC.Filter(truth_condition).Count().GetValue()

    # Count selected events
    selected = rdfMC.Filter(selection).Count().GetValue()

    # Count true events that are selected
    true_selected = rdfMC.Filter(f"{selection} && {truth_condition}").Count().GetValue()

    # Calculate efficiency and purity
    efficiency = true_selected / total_true if total_true > 0 else 0
    purity = true_selected / selected if selected > 0 else 0

    return efficiency, purity