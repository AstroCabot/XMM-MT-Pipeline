import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "outputs"
FITS = OUTPUT / "sensitivity_fits.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize(archive, baseline):
    wind = json.loads((ROOT / "analysis_06_solar_wind/outputs/result.json").read_text())
    prior = wind["heavy_ion_prior"]
    if prior["distribution"] != "uniform":
        raise ValueError("normalization check requires the current uniform flux prior")
    chains = {}
    for name, subdir in (("tan", ""), ("frincke", "water_frincke")):
        path = ROOT / "analysis_08_chemistry/outputs" / subdir / "chain.npz"
        with np.load(path) as chain:
            chains[name] = {
                "sha256": digest(path),
                "phi": chain["samples"][:, list(chain["coordinate_names"]).index("phi_heavy")],
                "q": chain["physical_draws"][:, list(chain["physical_names"]).index("Q_H2")].astype(float),
            }
    energy = baseline["mean_photon_energy_eV"]
    conversion = baseline["detector_K_ct_s_per_erg_cm2_s"]["mos2"]
    luminosity = baseline["regions"]["inner"]["luminosity_erg_s"]
    rows = []
    for fitted in archive["fits"]:
        scale = fitted["mean_photon_energy_eV"] * fitted["detector_K"]["mos2"] / (energy * conversion)
        sensitivity = {}
        for name, chain in chains.items():
            keep = (chain["phi"] / scale >= prior["lower"]) & (chain["phi"] / scale <= prior["upper"])
            quantiles = np.quantile(chain["q"][keep], [0.16, 0.50, 0.84])
            sensitivity[name] = {
                "retained_draws": int(keep.sum()),
                "total_draws": len(keep),
                "Q_H2_q16_median_q84": quantiles.tolist(),
                "median_change_percent": float(100 * (quantiles[1] / np.quantile(chain["q"], 0.5) - 1)),
            }
        rows.append({
            "variant": fitted["variant"],
            "chi2": fitted["chi2"],
            "dof": fitted["dof"],
            "delta_chi2": fitted["chi2"] - baseline["fit"]["chi2"],
            "mean_photon_energy_eV": fitted["mean_photon_energy_eV"],
            "energy_change_percent": 100 * (fitted["mean_photon_energy_eV"] / energy - 1),
            "luminosity_change_percent": 100 * (fitted["regions"]["inner"]["luminosity_erg_s"] / luminosity - 1),
            "MOS2_ECF_erg_cm2_per_count": 1 / fitted["detector_K"]["mos2"],
            "scale_energy_over_ECF": scale,
            "normalization_check": sensitivity,
        })
    checks = {path: digest(ROOT / path) == value if (ROOT / path).is_file() else None
              for path, value in archive["input_sha256"].items()}
    if any(value is False for value in checks.values()):
        raise RuntimeError("archived sensitivity fits differ from the current spectral inputs")
    return {
        "spectral_input_sha256_matches": checks,
        "chain_sha256": {name: chain["sha256"] for name, chain in chains.items()},
        "heavy_ion_prior": prior,
        "method": "At fixed extracted profile, a multiplicative change k in mean energy/ECF is absorbed by Phi -> Phi/k. The uniform prior gives constant weight inside its bounds and zero outside. The retained-sample check does not integrate newly accessible prior support or change detector exposure weights; it is not a new posterior fit.",
        "variants": rows,
    }


def refit(name, archive, baseline):
    import fit
    import model

    seed = next(row for row in archive["fits"] if row["variant"] == name)
    change, parameters = seed["change"], seed["parameters"]
    settings = json.loads((HERE / "parameters.json").read_text())["model"]
    settings["cx_collision_starts_keV_u"] = [parameters["cx_collision_keV_u"]] * 2
    configure0, link0, conversion0 = fit.configure, fit.link, fit.conversion
    captured = {}

    def configure(fitted, detector, area, config):
        configure0(fitted, detector, area, config)
        cx = model.components(fitted)[6]
        for key, value in change.items():
            if value == "free":
                model.free(cx.Hefrac, seed["Hefrac"], (0, 1))
            elif key == "Hefrac":
                model.fixed(cx.Hefrac, value)
            else:
                model.choice(getattr(cx, key), value)

    def link(groups, models):
        detectors, regions = link0(groups, models)
        first = model.components(models[0])
        for component, key in ((first[2], "lhb_norm"), (first[4], "halo_norm"), (first[5], "cxb_norm")):
            component.norm.values = parameters[key]
        first[6].temperature.values = parameters["cx_temperature_keV"]
        first[6].collnpar.values = parameters["cx_collision_keV_u"]
        for detector in ("mos1", "mos2"):
            model.components(detectors[detector])[1].factor.values = seed["crosscal"][detector]
        for region, fitted in regions.items():
            model.components(fitted)[6].norm.values = seed["regions"][region]["cx_norm_per_arcmin2"]
        if change.get("Hefrac") == "free":
            for fitted in models[1:]:
                model.components(fitted)[6].Hefrac.link = first[6].Hefrac
        return detectors, regions

    def conversion(xspec, models, regions, groups, band):
        result = conversion0(xspec, models, regions, groups, band)
        cx = model.components(models[0])[6]
        captured.update(Hefrac=float(cx.Hefrac.values[0]),
                        recombtype=int(cx.recombtype.values[0]),
                        acxmodel=int(cx.acxmodel.values[0]))
        return result

    with (OUTPUT / "spectra.tsv").open() as stream:
        groups = []
        for row in csv.DictReader(stream, delimiter="\t"):
            group = {key: float(row[key]) for key in ("exposure_s", "omega_arcmin2", "full_area_arcmin2", "response_throughput")}
            group.update(region=row["region"], det=row["det"])
            group.update({key: OUTPUT / "stack" / row[key] for key in ("source_pha", "background_pha", "rmf", "arf")})
            groups.append(group)
    fit.configure, fit.link, fit.conversion = configure, link, conversion
    try:
        result = fit.fit(groups, settings, baseline["distance_au"], os.environ["ATOMDB"], tuple(baseline["fit_band_keV"]))
    finally:
        fit.configure, fit.link, fit.conversion = configure0, link0, conversion0
    import importlib.metadata

    seed.update(
        chi2=result["statistic"], dof=result["dof"], parameters=result["parameters"],
        crosscal=result["crosscal"], mean_photon_energy_eV=result["mean_photon_energy_eV"],
        detector_K={g["det"]: result["region_detector_K"]["inner"][g["det"]] / g["response_throughput"] for g in groups if g["region"] == "inner"},
        regions={name: {key: row[key] for key in ("luminosity_erg_s", "cx_norm_per_arcmin2")} for name, row in result["regions"].items()},
        versions={package: importlib.metadata.version(package) for package in ("acx2", "pyatomdb", "numpy")},
        **captured,
    )
    FITS.write_text(json.dumps(archive, indent=2) + "\n")


def main():
    archive = json.loads(FITS.read_text())
    parser = argparse.ArgumentParser(description="Summarize saved donor/capture sensitivity fits; optionally refit one variant.")
    parser.add_argument("--refit", choices=[row["variant"] for row in archive["fits"]])
    args = parser.parse_args()
    baseline = json.loads((OUTPUT / "result.json").read_text())
    summary = summarize(archive, baseline)
    if args.refit:
        refit(args.refit, archive, baseline)
        summary = summarize(archive, baseline)
    (OUTPUT / "sensitivity.json").write_text(json.dumps(summary, indent=2) + "\n")
    for row in summary["variants"]:
        print(f"{row['variant']:24} E={row['mean_photon_energy_eV']:.3f} eV  dL={row['luminosity_change_percent']:+.3f}%")


if __name__ == "__main__":
    main()
