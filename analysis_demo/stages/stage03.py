import csv
import os
import sys

from _common import DATA, REPRO, read_json, restore_env, stage_dir, write_json

sys.path.insert(0, str(REPRO / "analysis_03_spectrum"))

import fit
import run as run3

restore_env()

import matplotlib.figure

SPEC = DATA / "spectrum"
FIT_BAND = (0.35, 1.1)
PROFILE_TARGETS = {"cx_norm_inner", "cx_norm_outer"}

_full_targets = fit.profile_targets
fit.profile_targets = lambda *args: {
    name: value
    for name, value in _full_targets(*args).items()
    if name in PROFILE_TARGETS
}
fit.profiled_fluxes = lambda *args, **kwargs: {}

_legend = matplotlib.figure.Figure.legend


def pn_legend(self, *args, **kwargs):
    handles = [
        handle
        for handle in kwargs.pop("handles", [])
        if handle.get_label().strip() not in ("", "data + model")
    ]
    kwargs["ncol"] = len(handles)
    return _legend(self, *args, handles=handles, **kwargs)


def load_groups():
    with (SPEC / "groups.tsv").open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    groups = []
    for row in rows:
        rmf = SPEC / (row["rmf"] + ".gz")
        groups.append(
            {
                "region": row["region"],
                "det": row["det"],
                "source_pha": SPEC / row["source_pha"],
                "background_pha": SPEC / row["background_pha"],
                "rmf": rmf if rmf.is_file() else SPEC / row["rmf"],
                "arf": SPEC / row["arf"],
                "n_frames": int(row["n_frames"]),
                "exposure_s": float(row["exposure_s"]),
                "omega_arcmin2": float(row["omega_arcmin2"]),
                "response_throughput": float(row["response_throughput"]),
                "full_area_arcmin2": float(row["full_area_arcmin2"]),
            }
        )
    return groups


def main():
    output = stage_dir(3)
    groups = load_groups()
    geometry = read_json(SPEC / "geometry.json")
    settings = read_json(REPRO / "analysis_03_spectrum/parameters.json")["model"]
    software = read_json(REPRO / "analysis_01_reduction/settings.json")
    atomdb = os.environ.get("ATOMDB", software["atomdb"])
    primary = fit.fit(
        groups,
        settings,
        geometry["distance_au"],
        atomdb,
        FIT_BAND,
        profile=True,
    )

    region_k = run3.on_axis_factors(primary["region_detector_K"], groups)
    k_pn = region_k["inner"]["pn"]
    ratio = geometry["archived_K_pn_over_K_mos2"]
    regions = primary["regions"]
    intervals = primary["profile_intervals_1sigma"]
    for name, region in regions.items():
        interval = intervals.get(f"cx_norm_{name}")
        if interval:
            center = region["cx_norm_per_arcmin2"]
            region["flux_1sigma_erg_cm2_s"] = [
                region["flux_erg_cm2_s"] * value / center for value in interval
            ]
            region["luminosity_1sigma_erg_s"] = [
                region["luminosity_erg_s"] * value / center for value in interval
            ]
        region["mos2_equivalent_surface_brightness_ct_s_arcmin2"] = (
            region["surface_brightness_erg_cm2_s_arcmin2"]
            * primary["region_detector_K"][name]["pn"]
            / ratio
        )
        if name != "inner":
            region["image_response_ratio_to_inner"] = (
                primary["region_detector_K"][name]["pn"]
                / primary["region_detector_K"]["inner"]["pn"]
            )
    result = {
        "distance_au": geometry["distance_au"],
        "heliocentric_distance_au": geometry["heliocentric_distance_au"],
        "fit_band_keV": list(FIT_BAND),
        "model": {
            "expression": fit.EXPRESSION,
            "abundances": "wilm",
            "cross_sections": "vern",
            "recombtype": 1,
            "acxmodel": 8,
            "collntype": 1,
            "free_parameters": primary["n_free_parameters"],
        },
        "fit": {
            "chi2": primary["statistic"],
            "dof": primary["dof"],
            "parameters": primary["parameters"],
            "crosscal": primary["crosscal"],
            "profile_intervals_1sigma": primary["profile_intervals_1sigma"],
        },
        "detector_K_ct_s_per_erg_cm2_s": {
            "pn": k_pn,
            "mos1": k_pn,
            "mos2": k_pn / ratio,
        },
        "response_throughput": {
            group["det"]: group["response_throughput"]
            for group in groups
            if group["region"] == "inner"
        },
        "mos2_ecf_erg_cm2_per_count": ratio / k_pn,
        "mos2_ecf_sigma_ln": geometry["mos2_ecf_sigma_ln"],
        "mean_photon_energy_eV": primary["mean_photon_energy_eV"],
        "regions": regions,
    }
    folded = []
    for group, arrays in zip(primary["groups"], primary["arrays"], strict=True):
        for index in range(len(arrays["energy_keV"])):
            folded.append(
                {
                    "region": group["region"],
                    "det": group["det"],
                    **{name: values[index] for name, values in arrays.items()},
                }
            )
    write_json(output / "result.json", result)
    blank = ("none", "", "")
    run3.figure.CAMERAS = {
        "pn": run3.figure.CAMERAS["pn"],
        "mos1": blank,
        "mos2": blank,
    }
    matplotlib.figure.Figure.legend = pn_legend
    run3.figure.draw(folded, output / "figure2.png", output / "figure2.pdf", FIT_BAND)
    print(
        f"stage 03: PN-only fit chi2/dof = {primary['statistic']:.1f}/{primary['dof']}"
        f", inner L_X = {regions['inner']['luminosity_erg_s']:.3e} erg/s"
        f", <E> = {primary['mean_photon_energy_eV']:.1f} eV"
    )


if __name__ == "__main__":
    main()
