import ROOT
from lib.cuts.cut_evaluation import SingleCut, CutEvaluation
import numpy as np
from tqdm import tqdm

from dataclasses import dataclass

class CutOptimization:

  cut_values: list[float] = []
  efficiencies: list[float] = []
  purities: list[float] = []
  efficiencies_close: list[float] = []
  purities_close: list[float] = []

  def __init__(self, rdf: ROOT.RDataFrame, config: str, signal_filter: str = "(mctruth == 1 || mctruth == 0)"):
    self.signal_filter = signal_filter
    self.config = CutEvaluation(config).config

    self.rdf = rdf.Filter("mctruth != -1 && mcflag == 1")

    self.rdf_signal = self.rdf.Filter(f"{self.signal_filter} && mcflag == 1")
    self.rdf_background = self.rdf.Filter(f"!({self.signal_filter}) && mcflag == 1")

  def efficiency_purity_plot(self, prefix_cut: str = "", optimized_cut: SingleCut = None, num_points: int = 10, graph_path: str = ""):
    # Tutaj implementacja metody do tworzenia wykresu efektywności vs czystości
    
    if optimized_cut is None:
      print("Nothing to optimize, skipping plot.")
      return None

    # Przygotowanie zakresu optymalizacji
    if optimized_cut.optimization_range is not None:
      step = (optimized_cut.optimization_range[1] - optimized_cut.optimization_range[0]) / num_points
      start = optimized_cut.optimization_range[0]
      end = optimized_cut.optimization_range[1]
      self.cut_values = [start + i * step for i in range(num_points + 1)]
    else:
      print("No optimization range provided, skipping plot.")
      return None

    self.total_signal = self.rdf_signal.Count().GetValue()
    self.total_signal_close = self.rdf_signal.Filter("abs(deltaT) < 15").Count().GetValue()

    # Lista do przechowywania "obietnic" wyników (RResultPtr)
    booked_results = []
    all_handles = [] # Ta lista będzie potrzebna dla RunGraphs

    print(f"\nOptimizing cut: {optimized_cut.variable}")

    for value in self.cut_values:
        # Uwaga: Używamy formatowania f-string dla precyzji zmiennoprzecinkowej
        cut_string = optimized_cut.to_root_string(self.config).replace(str(optimized_cut.value), f"{value:.6f}")
        full_cut = f"{prefix_cut} && {cut_string}" if prefix_cut else cut_string

        # Ważne: NIE wywołujemy tutaj .GetValue()!
        res_signal = self.rdf_signal.Filter(full_cut).Count()
        res_bkg = self.rdf_background.Filter(full_cut).Count()
        res_signal_close = self.rdf_signal.Filter(full_cut + f" && abs(deltaT) < 15").Count()
        res_bkg_close = self.rdf_background.Filter(full_cut + f" && abs(deltaT) < 15").Count()
        
        booked_results.append((res_signal, res_bkg, res_signal_close, res_bkg_close))

        all_handles.append(res_signal)
        all_handles.append(res_bkg)
        all_handles.append(res_signal_close)
        all_handles.append(res_bkg_close)

    # 3. ETAP WYKONANIA (Execution)
    print(f"Running graphs for {len(all_handles)} nodes...")
    
    # Przekazujemy listę uchwytów do funkcji
    ROOT.RDF.RunGraphs(all_handles)

    for res_s, res_b, res_s_close, res_b_close in booked_results:
        signal_pass = res_s.GetValue()
        background_pass = res_b.GetValue()
        signal_pass_close = res_s_close.GetValue()
        background_pass_close = res_b_close.GetValue()

        efficiency = signal_pass / self.total_signal if self.total_signal > 0 else 0
        purity = signal_pass / (signal_pass + background_pass) if (signal_pass + background_pass) > 0 else 0

        efficiency_close = signal_pass_close / self.total_signal_close if self.total_signal_close > 0 else 0
        purity_close = signal_pass_close / (signal_pass_close + background_pass_close) if (signal_pass_close + background_pass_close) > 0 else 0

        self.efficiencies.append(efficiency)
        self.purities.append(purity)
        self.efficiencies_close.append(efficiency_close)
        self.purities_close.append(purity_close)

    # Przygotowanie wykresu
    self.mg, self.eff_graph, self.pur_graph = self._prepare_plot(close=False)

    img_name = f"optimization_plot_{optimized_cut.variable}_{num_points}.svg"
    img_name = f"{graph_path}/{img_name}" if graph_path else f"{img_name}"

    canvas = ROOT.TCanvas("c1", "Optimization", 800, 600)
    self.mg.Draw("ALP")

    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    legend.AddEntry(self.eff_graph, "Efficiency", "lp")
    legend.AddEntry(self.pur_graph, "Purity", "lp")
    legend.Draw()

    canvas.Update()
    canvas.SaveAs(img_name)

    # Close
    self.mg, self.eff_graph, self.pur_graph = self._prepare_plot(close=True)

    img_name = f"optimization_plot_{optimized_cut.variable}_{num_points}_close.svg"
    img_name = f"{graph_path}/{img_name}" if graph_path else f"{img_name}"

    canvas = ROOT.TCanvas("c2", "Optimization", 800, 600)
    self.mg.Draw("ALP")

    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    legend.AddEntry(self.eff_graph, "Efficiency", "lp")
    legend.AddEntry(self.pur_graph, "Purity", "lp")
    legend.Draw()

    canvas.Update()
    canvas.SaveAs(img_name)

  def _prepare_plot(self, close: bool):
    if close:
        efficiency_graph = ROOT.TGraph(len(self.cut_values), np.array(self.cut_values, dtype='d'), np.array(self.efficiencies_close, dtype='d'))
        purity_graph = ROOT.TGraph(len(self.cut_values), np.array(self.cut_values, dtype='d'), np.array(self.purities_close, dtype='d'))
    else:
        efficiency_graph = ROOT.TGraph(len(self.cut_values), np.array(self.cut_values, dtype='d'), np.array(self.efficiencies, dtype='d'))
        purity_graph = ROOT.TGraph(len(self.cut_values), np.array(self.cut_values, dtype='d'), np.array(self.purities, dtype='d'))

    # Stylizacja
    efficiency_graph.SetLineColor(ROOT.kBlue)
    efficiency_graph.SetMarkerColor(ROOT.kBlue)
    efficiency_graph.SetMarkerStyle(20)
    efficiency_graph.SetTitle("Efficiency")

    purity_graph.SetLineColor(ROOT.kRed)
    purity_graph.SetMarkerColor(ROOT.kRed)
    purity_graph.SetMarkerStyle(21)
    purity_graph.SetTitle("Purity")

    # Nakładanie za pomocą TMultiGraph
    mg = ROOT.TMultiGraph()
    mg.Add(efficiency_graph)
    mg.Add(purity_graph)
    mg.SetMinimum(0.0)
    mg.SetMaximum(1.1)
    mg.SetTitle(";Threshold value [a.u.];Rate [-]")

    return mg, efficiency_graph, purity_graph

