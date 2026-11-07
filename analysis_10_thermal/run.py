import argparse
import json
from pathlib import Path

import numpy as np
from thermal import DAY_S, distance, simulate

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
U_KG = 1.660_539_066_60e-27
MASS_U = {"H2O": 18.01528, "CO2": 44.0095, "CO": 28.0101, "H2": 2.01588}
HEAT_CAPACITY_KB = {"H2O": 4.0, "CO2": 3.5, "CO": 3.5, "H2": 3.5}


def q3(values):
    lo, mid, hi = np.quantile(values, (0.16, 0.5, 0.84))
    return {"p16": float(lo), "median": float(mid), "p84": float(hi)}


def heated_depths(conductivity, model, orbit, xmm_days, temperatures):
    out = simulate(conductivity, model, *orbit, np.ceil(xmm_days) + 1, temperatures)
    return [
        np.array([np.interp(day, out["days"], column) for column in out["reached_m"].T])
        for day in (0.0, xmm_days)
    ]


def history_integral(orbit, r_xmm, start, end):
    days = np.linspace(start, end, 2001)
    ratio = distance(*orbit, days) / r_xmm
    betas = np.linspace(0.0, 15.0, 1501)
    integrals = [np.trapezoid(ratio**-beta, days * DAY_S) for beta in betas]
    return lambda beta: np.interp(beta, betas, integrals)


def budget(chain, orbit, xmm_days, parameters, molecules_per_m):
    with np.load(chain, allow_pickle=False) as data:
        names = [str(name) for name in data["physical_names"]]
        draws = dict(zip(names, np.asarray(data["physical_draws"], float).T))
    r_xmm = distance(*orbit, xmm_days)
    since_perihelion = history_integral(orbit, r_xmm, 0.0, xmm_days)
    combi_window = history_integral(orbit, r_xmm, *parameters["combi_window_days"])
    water = draws["Q_H2O"] * since_perihelion(draws["beta_H2O"])
    hydrogen = draws["Q_H2"] * since_perihelion(draws["beta_H2"])
    scale = distance(*orbit, parameters["iram_days"]) / r_xmm
    flux = {
        name: draws[f"Q_{name}"] * scale ** -draws[f"beta_{name}"] * MASS_U[name]
        for name in MASS_U
    }
    total = sum(flux.values())
    heavy = total - flux["H2"]
    enthalpy = {
        name: value * HEAT_CAPACITY_KB[name] / MASS_U[name]
        for name, value in flux.items()
    }
    mixed = sum(enthalpy.values()) / total
    unmixed = (mixed * total - enthalpy["H2"]) / heavy
    return {
        "H2O_molecules": q3(water),
        "H2_molecules": q3(hydrogen),
        "H2_over_H2O": q3(hydrogen / water),
        "H2O_layer_m": q3(water / molecules_per_m),
        "combi_window_H2O_kg": q3(
            draws["Q_H2O"] * combi_window(draws["beta_H2O"]) * MASS_U["H2O"] * U_KG
        ),
        "iram_H2_mass_fraction": q3(flux["H2"] / total),
        "iram_mixed_speed_ratio": q3((mixed / unmixed) ** 0.5),
    }


def run(task08_folder, steep_folder, task09_path, parameters_path, output_path):
    output_resolved = output_path.resolve()
    try:
        output_resolved.relative_to(ROOT.resolve())
    except ValueError as error:
        raise ValueError("output must remain inside analysis_repro") from error
    inputs = (
        task08_folder / "result.json",
        task08_folder / "chain.npz",
        steep_folder / "chain.npz",
        task09_path,
        parameters_path,
    )
    if any(
        output_resolved.exists() and path.exists() and output_resolved.samefile(path)
        for path in inputs
    ):
        raise ValueError("output must differ from every input")
    output_path.unlink(missing_ok=True)
    parameters = json.loads(parameters_path.read_text())
    model = parameters["model"]
    fit = json.loads((task08_folder / "result.json").read_text())
    area = json.loads(task09_path.read_text())["nucleus"]["surface_area_m2"]
    orbit = (fit["orbit"]["perihelion_q_au"], fit["orbit"]["eccentricity"])
    xmm_days = fit["xmm_observation_jd_tdb"] - fit["orbit"]["perihelion_jd_tdb"]
    water_mass = MASS_U["H2O"] * U_KG
    molecules_per_m = model["density_kg_m3"] * parameters["ice_fraction"] * area
    molecules_per_m /= water_mass
    budgets = {
        name: budget(folder / "chain.npz", orbit, xmm_days, parameters, molecules_per_m)
        for name, folder in (("adopted", task08_folder), ("steep", steep_folder))
    }
    ratios = parameters["h2_per_h2o"]
    lowest, highest = parameters["release_K"]
    temperatures = np.arange(lowest, highest + 0.5)
    depths, conduction = {}, {}
    for conductivity in parameters["conductivity_W_m_K"]:
        perihelion, xmm = heated_depths(
            conductivity, model, orbit, xmm_days, temperatures
        )
        name = f"{conductivity:g}_W_m_K"
        depths[name] = {
            f"{t:g}_K": round(float(np.interp(t, temperatures, xmm)), 3)
            for t in parameters["isotherms_K"]
        }
        layer = float((xmm - perihelion).max())
        heated = molecules_per_m * layer
        above = molecules_per_m * depths[name]["35_K"]
        conduction[name] = {
            "heated_m": round(layer, 3),
            "heated_H2O_molecules": heated,
            "max_H2_share": {
                history: max(ratios) * heated / water["H2_molecules"]["median"]
                for history, water in budgets.items()
            },
            "above_35_K_H2_share": {
                history: {
                    f"{ratio:g}": ratio * above / water["H2_molecules"]["median"]
                    for ratio in ratios
                }
                for history, water in budgets.items()
            },
        }
    days = np.linspace(0.0, xmm_days, 2001)
    absorbed = (1 - model["albedo"]) * model["solar_constant_W_m2"] * area / 4
    energy = np.trapezoid(absorbed / distance(*orbit, days) ** 2, days * DAY_S)
    result = {
        "orbit": fit["orbit"],
        "xmm_days": xmm_days,
        "surface_area_m2": area,
        "depth_reached_m": depths,
        "conduction": conduction,
        "sunlight_H2O_molecules": energy / parameters["latent_heat_J_kg"] / water_mass,
        "water": budgets,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--task08", type=Path, default=ROOT / "analysis_08_chemistry/outputs"
    )
    parser.add_argument(
        "--steep",
        type=Path,
        default=ROOT / "analysis_08_chemistry/outputs/water_frincke",
    )
    parser.add_argument(
        "--task09", type=Path, default=ROOT / "analysis_09_recoil/outputs/result.json"
    )
    parser.add_argument("--parameters", type=Path, default=HERE / "parameters.json")
    parser.add_argument("--output", type=Path, default=HERE / "outputs/result.json")
    args = parser.parse_args()
    run(args.task08, args.steep, args.task09, args.parameters, args.output)


if __name__ == "__main__":
    main()
