import json
import os
from typing import Optional, List, Dict
from pydantic import BaseModel, Field, field_validator

class SingleCut(BaseModel):
    variable: str
    value: float
    bias: float = 0.0
    operator: str = Field(default="<")
    active: bool = True
    order: int = 0

    # Przykład walidacji: upewnij się, że nie masz ujemnego Chi2
    @field_validator('value')
    def value_must_be_positive_for_chi2(cls, v, values):
        if 'Chi2' in values.get('variable', '') and v < 0:
            raise ValueError('Chi2 cut value must be positive!')
        return v

class CutScenario(BaseModel):
    name: str
    description: Optional[str] = None
    cuts: List[SingleCut]

    def to_root_string(self) -> str:
        active_cuts = [c for c in self.cuts if c.active]
        cut_strings = []
        for c in active_cuts:
            if c.operator == "abs_lt":
                cut_strings.append(f"abs({c.variable} - {c.bias}) < {c.value}")
            else:
                cut_strings.append(f"({c.variable} - {c.bias}) {c.operator} {c.value}")
        return " && ".join([f"({s})" for s in cut_strings])

class FullConfig(BaseModel):
    scenarios: Dict[str, CutScenario]

class CutEvaluation:
    def __init__(self, cut_file):
        if not os.path.exists(cut_file):
            raise FileNotFoundError(f"Cut file '{cut_file}' not found.")

        try:
          with open(cut_file, 'r') as file:
            data = json.load(file)
        except json.JSONDecodeError as e:
          raise ValueError(f"Error parsing cut file '{cut_file}': {e}")

        self.config = FullConfig(scenarios=data.get("scenarios", {}))

    def get_cut(self, scenario_name: str) -> str:
        if scenario_name in self.config.scenarios:
            return self.config.scenarios[scenario_name].to_root_string()
        raise KeyError(f"Scenario {scenario_name} does not exist in JSON!")