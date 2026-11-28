import argparse
import json

import numpy as np

import run
from config import CONFIG, SHARED
from inference import coordinates, unpack

FITS = {
    "parents_shallow": ("no_direct_h2", False, "beta"),
    "parents_steep": ("no_direct_h2/water_frincke", False, "beta_steep"),
    "radiolytic_shallow": ("", True, "beta"),
    "radiolytic_steep": ("water_frincke", True, "beta_steep"),
}
RADII = (0.75, 12.0)
SEED = 20261127
GROUPS = ("H2_direct", "H", "O", "OH", "parents", "C", "remaining")


def area_weights(edges, km_per_arcmin, order=24):
    radii = edges / km_per_arcmin
    bounds = np.r_[RADII[0], radii[(radii > RADII[0]) & (radii < RADII[1])], RADII[1]]
    node, weight = np.polynomial.legendre.leggauss(order)
    half = 0.5 * np.diff(bounds)
    radius = (0.5 * (bounds[1:] + bounds[:-1])[:, None] + half[:, None] * node).ravel()
    weight = (half[:, None] * weight).ravel() * radius
    rho = radius * km_per_arcmin
    inner, outer = edges[:-1, None], edges[1:, None]
    chord = np.sqrt(np.maximum(outer**2 - rho**2, 0))
    chord -= np.sqrt(np.maximum(inner**2 - rho**2, 0))
    return chord @ weight


def load(inputs, operator, folder, direct, key):
    beta = SHARED["spherex"]["parent_priors"]["H2O"][key]
    CONFIG["direct_h2"] = direct
    CONFIG["priors"]["beta_H2O"] = [beta["mean"], beta["sigma"]]
    posterior, _ = run.build_model(inputs, operator)
    edges = posterior.transport.edges
    area = area_weights(edges, inputs["geometry"]["km_per_arcmin"])
    with np.load(run.OUTPUT / folder / "chain.npz") as chain:
        if tuple(chain["coordinate_names"]) != coordinates(CONFIG):
            raise ValueError(f"{folder}: chain coordinates differ from the model")
        index = list(chain["best_labels"]).index("overall_map")
        best = chain["best_coordinates"][index]
        samples = chain["samples"].astype(float)
    return posterior, area, best, samples


def shares(theta, posterior, area):
    transport = posterior.transport
    _, rates, beta, speeds, _ = unpack(theta, posterior.config, transport.orbit)
    number, weighted = dict.fromkeys(GROUPS, 0.0), dict.fromkeys(GROUPS, 0.0)
    for parent, chains in posterior.chains.items():
        if rates[parent] == 0:
            continue
        for chain in chains:
            species = chain[0][-1][0]
            group = "parents" if species in ("H2O", "CO2", "CO") else species
            if group not in GROUPS:
                group = "H2_direct" if species == parent == "H2" else "remaining"
            density = transport.density(chain, beta[parent], speeds)
            column = rates[parent] * density @ area
            number[group] += column
            weighted[group] += column * transport.network[species]["sigma_cm2"]
    return {
        label: {group: value / sum(values.values()) for group, value in values.items()}
        for label, values in (("number", number), ("sigma_weighted", weighted))
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--draws", type=int, default=300)
    args = parser.parse_args()
    inputs = run.load_inputs()
    operator = run.load_psf(run.OUTPUT / "psf.npz", inputs, run.psf_dependencies())
    rng = np.random.default_rng(SEED)
    for name, (folder, *model) in FITS.items():
        posterior, area, best, samples = load(inputs, operator, folder, *model)
        selected = rng.choice(len(samples), args.draws, replace=False)
        draws = [shares(theta, posterior, area) for theta in samples[selected]]
        result = {
            "radius_arcmin": list(RADII),
            "seed": SEED,
            "draws": args.draws,
            "percentiles": [16, 50, 84],
            "map": shares(best, posterior, area),
            "draw_percentiles": {
                label: {
                    group: np.quantile(
                        [draw[label][group] for draw in draws], (0.16, 0.5, 0.84)
                    ).tolist()
                    for group in GROUPS
                }
                for label in ("number", "sigma_weighted")
            },
        }
        (run.OUTPUT / folder / "columns.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
        print(name, json.dumps(result["map"]["number"]), flush=True)


if __name__ == "__main__":
    main()
