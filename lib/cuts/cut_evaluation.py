import json
import os
from typing import Optional, List, Dict, Union, Optional
from pydantic import BaseModel, Field, field_validator
from itertools import groupby

class SingleCut(BaseModel):
    variable: str
    value: float
    bias: Union[float, str] = None
    operator: str = Field(default="<")
    active: bool = True
    order: int = 0
    apply_only_in: Optional[str] = None

    def to_root_string(self, config: Optional['FullConfig'] = None) -> str:
        """Generuje string dla pojedynczego cięcia."""

        if self.bias is None:
          if self.operator == "abs_lt":
              core = f"abs({self.variable}) < {self.value}"
          else:
              core = f"{self.variable} {self.operator} {self.value}"
        else:
          if self.operator == "abs_lt":
              core = f"abs({self.variable} - {self.bias}) < {self.value}"
          else:
              core = f"({self.variable} - {self.bias}) {self.operator} {self.value}"

        core = f"({core})"

        if self.apply_only_in and config:
            fv_formula = config.get_region_formula(self.apply_only_in)
            # Logika: (W_REGIONIE && CIĘCIE) || !W_REGIONIE
            return f"(({fv_formula} && {core}) || !{fv_formula})"
        
        return core

class CutScenario(BaseModel):
    name: str
    description: Optional[str] = None
    cuts: List[SingleCut]

    @property
    def active_cuts(self) -> List[SingleCut]:
        # Sortujemy po polu 'order', aby mieć pewność co do kolejności
        return sorted([c for c in self.cuts if c.active], key=lambda x: x.order)

    @property
    def grouped_cuts(self) -> List[List[SingleCut]]:
        # Sortujemy, a potem grupujemy po wartości order
        sorted_active = sorted([c for c in self.cuts if c.active], key=lambda x: x.order)
        return [list(group) for key, group in groupby(sorted_active, lambda x: x.order)]

    def all_active_cuts(self, config: 'FullConfig') -> str:
        """Łączy wszystkie aktywne cięcia, uwzględniając ich FV."""
        active = sorted([c for c in self.cuts if c.active], key=lambda x: x.order)
        # Przekazujemy config do każdego cięcia
        result = " && ".join([f"({c.to_root_string(config)})" for c in active])
        return f"({result})" if result else ""

    def get_single_cut_by_order(self, order_index: int, config: 'FullConfig') -> str:
          """Zwraca cięcia o danym numerze 'order' z uwzględnieniem ich FV."""
          
          # Pobieramy wszystkie aktywne cięcia z tym samym numerem order
          cuts_in_order = [c for c in self.active_cuts if c.order == order_index]
          
          if not cuts_in_order:
              return ""
          
          cut_strings = []
          for c in cuts_in_order:
              # Kluczowy moment: wywołujemy to_root_string przekazując config
              raw_cut = c.to_root_string(config)
              cut_strings.append(f"({raw_cut})")
          
          # Łączymy (np. u i v z tym samym order) operatorem AND
          result = " && ".join(cut_strings)
          return f"({result})"

class FullConfig(BaseModel):
    scenarios: Dict[str, CutScenario]
    fiducial_volumes: Dict[str, List[SingleCut]]

    def get_region_formula(self, region_name: str) -> str:
        if region_name not in self.fiducial_volumes:
            return "1"
        
        region_cuts = self.fiducial_volumes[region_name]
        # Dla regionów FV wywołujemy to_root_string bez configu (żeby uniknąć rekurencji)
        active_region_cuts = [f"({c.to_root_string()})" for c in region_cuts if c.active]
        result = " && ".join(active_region_cuts) if active_region_cuts else "1"
        return f"({result})" if result else "1"

class CutEvaluation:
    def __init__(self, cut_file):
        if not os.path.exists(cut_file):
            raise FileNotFoundError(f"Cut file '{cut_file}' not found.")

        try:
          with open(cut_file, 'r') as file:
            data = json.load(file)
        except json.JSONDecodeError as e:
          raise ValueError(f"Error parsing cut file '{cut_file}': {e}")

        self.config = FullConfig(**data)

    def get_all_cuts(self, scenario_name: str) -> str:
        if scenario_name in self.config.scenarios:
            return self.config.scenarios[scenario_name].all_active_cuts(self.config)
        raise KeyError(f"Scenario {scenario_name} does not exist in JSON!")
  
    def get_single_cut(self, scenario_name: str, order_index: int) -> str:
        if scenario_name in self.config.scenarios:
            return self.config.scenarios[scenario_name].get_single_cut_by_order(order_index, self.config)
        raise KeyError(f"Scenario {scenario_name} does not exist in JSON!")