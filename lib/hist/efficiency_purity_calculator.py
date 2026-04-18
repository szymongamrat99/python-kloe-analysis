import ROOT


class EfficiencyPurityCalculator:
  def __init__(self, rdf: ROOT.RDataFrame, truth_condition: str, bin_width=1.0, xmin=-100, xmax=100):
    self.time_diff_MC = "KaonChTimeCMMC - KaonNeTimeCMMC"
    self.time_diff_rec = "KaonChTimeCMSignalFit - KaonNeTimeCMSignalFit"

    define = "Redefine" if rdf.HasColumn("deltaT_MC") else "Define"
    self.rdf = getattr(rdf, define)("deltaT_MC", self.time_diff_MC)
    define = "Redefine" if self.rdf.HasColumn("deltaT") else "Define"
    self.rdf = getattr(self.rdf, define)("deltaT", self.time_diff_rec)
    self.rdfMC = self.rdf.Filter("mcflag == 1")
    self.rdfData = self.rdf.Filter("mcflag == 0")

    self.xmin = xmin
    self.xmax = xmax
    self.nbins = int((self.xmax - self.xmin) / bin_width)

    if self.nbins % 2 == 0:
        self.nbins += 1  # Ensure odd number of bins for symmetry around zero

    self.hist_models = {}

    self.hist_models["total_signal"] = ("total_signal", ";#Deltat [#tau_{S}];Efficiency [-]", self.nbins, self.xmin, self.xmax)
    self.hist_models["selected_signal"] = ("selected_signal", ";#Deltat [#tau_{S}];Efficiency [-]", self.nbins, self.xmin, self.xmax)

    self.truth_condition = truth_condition

  def efficiency_histo_per_time_difference(self, selection: str, bin_width=1.0, xmin=-100, xmax=100):

    histSignalTotal = self.rdfMC.Filter(self.truth_condition).Histo1D(self.hist_models["total_signal"], "deltaT_MC")

    selection = selection + " && " + self.truth_condition if selection else self.truth_condition
    histSelected = self.rdfMC.Filter(selection).Histo1D(self.hist_models["selected_signal"], "deltaT_MC")

    self.histEfficiency = ROOT.TEfficiency(histSelected.GetPtr(), histSignalTotal.GetPtr())

    self.histEfficiency.SetTitle(";#Deltat [#tau_{S}];Efficiency [-]")
    self.histEfficiency.SetName("efficiency_per_time_difference")

    return self.histEfficiency

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