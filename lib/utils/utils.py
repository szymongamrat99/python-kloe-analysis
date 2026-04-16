import ROOT
import json
import os


def load_files_into_chain(config_path: str, chain: any) -> None:
    with open(config_path) as f:
            config = json.load(f)
        
    finished = config.get("finishedRootFiles", [])
    print(f"Loading files into chain from config: {finished}")

    for finished_path in finished:
        with open(finished_path) as f:
            for line in f:
                path = line.strip()

                if not path:
                    continue

                chain.Add(path)

 
def get_chann_dicts() -> dict:
    channName = dict(ROOT.KLOE.channName)
    channColor = dict(ROOT.KLOE.channColor)

    channName = {int(k): str(v) for k, v in channName.items()}
    channColor = {str(k): int(v) for k, v in channColor.items()}
    return channName, channColor
