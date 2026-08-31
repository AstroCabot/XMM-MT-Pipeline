import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHARED = json.loads((HERE.parent / "parameters.json").read_text())
CONFIG = json.loads((HERE / "parameters.json").read_text())
SPHEREX = SHARED["spherex"]
CONFIG["production_reference_r_au"] = SPHEREX["reference_distance_au"]
CONFIG["cross_section_cm2"]["H2"] = SHARED["h2_single_capture_cross_section_cm2"]
for name, values in SPHEREX["parent_priors"].items():
    rate = values["reference_rate_s-1"]
    beta = values["beta"]
    CONFIG["priors"][f"Qref_{name}"] = [rate["median"], rate["sigma_ln"]]
    CONFIG["priors"][f"beta_{name}"] = [beta["mean"], beta["sigma"]]
