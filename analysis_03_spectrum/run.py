import argparse
import csv
import json
import os
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RUNTIME = ROOT / "work" / "runtime" / "analysis_03_spectrum"
for name in ("tmp", "cache", "config", "matplotlib"):
    (RUNTIME / name).mkdir(parents=True, exist_ok=True)
os.environ.update(
    TMPDIR=str(RUNTIME / "tmp"),
    XDG_CACHE_HOME=str(RUNTIME / "cache"),
    XDG_CONFIG_HOME=str(RUNTIME / "config"),
    MPLCONFIGDIR=str(RUNTIME / "matplotlib"),
)

import numpy as np

import extract
import figure
import fit
import spectrum
from inputs import (
    FIT_BAND_KEV,
    OUTPUT,
    PARAMS,
    REGIONS,
    SHARED,
    SOFTWARE,
    WORK,
    load_inputs,
    reload_groups,
    validate_extraction_frames,
    validate_fresh_inputs,
)

FIT_PRODUCTS = ("result.json", "folded.tsv", "figure2.png", "figure2.pdf")
EXTRACTION_PRODUCTS = ("extraction.tsv", "spectra.tsv")
ABSOLUTE_EFFECTIVE_AREA_FRACTION = 0.10


def write_tsv(path, columns, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def ecf_sigma_ln(crosscal, interval):
    bounds = np.asarray(interval, float)
    if (
        crosscal <= 0
        or bounds.shape != (2,)
        or np.any(bounds <= 0)
        or not bounds[0] <= crosscal <= bounds[1]
    ):
        raise ValueError("invalid MOS2 profile interval")
    relative = np.max(np.abs(np.log(bounds / crosscal)))
    return float(np.hypot(np.log1p(ABSOLUTE_EFFECTIVE_AREA_FRACTION), relative))


def extract_all(frames, sources):
    jobs = [
        (
            frame,
            name,
            region,
            sources,
            PARAMS["channel_max"][frame["det"]],
            SHARED["pn_oot_fraction"],
            FIT_BAND_KEV,
            WORK,
        )
        for frame in frames
        for name, region in REGIONS.items()
    ]
    with ThreadPoolExecutor(max_workers=2) as mos, ThreadPoolExecutor(
        max_workers=1
    ) as pn:
        futures = [
            (pn if job[0]["det"] == "pn" else mos).submit(extract.one, *job)
            for job in jobs
        ]
        rows = [future.result() for future in futures]
    columns = (
        "frame",
        "det",
        "region",
        "exposure_s",
        "omega_arcmin2",
        "mean_exposure_s",
        "response_throughput",
        "qpb_fraction",
        "sp_fraction",
        "source_counts",
        "background_counts",
    )
    write_tsv(
        OUTPUT / "extraction.tsv",
        columns,
        [
            {
                name: f"{row[name]:.10g}" if isinstance(row[name], float) else row[name]
                for name in columns
            }
            for row in rows
        ],
    )
    return rows


def stacked_throughput(rows):
    weights = np.array([row["exposure_s"] * row["omega_arcmin2"] for row in rows])
    return float(
        np.average([row["response_throughput"] for row in rows], weights=weights)
    )


def stack_all(rows):
    output, groups = OUTPUT / "stack", []
    for region, (inner, outer) in REGIONS.items():
        full_area = np.pi * (outer**2 - inner**2) / 3600
        for detector in ("pn", "mos1", "mos2"):
            members = [
                row
                for row in rows
                if row["region"] == region and row["det"] == detector
            ]
            paths = spectrum.stack(
                members,
                output,
                f"{detector}-{region}",
                FIT_BAND_KEV,
                PARAMS["group_min_counts"],
                extract.run_tool,
            )
            source, background, rmf, arf, exposure, backscal, omega = paths
            groups.append(
                {
                    "region": region,
                    "det": detector,
                    "source_pha": source,
                    "background_pha": background,
                    "rmf": rmf,
                    "arf": arf,
                    "n_frames": len(members),
                    "exposure_s": exposure,
                    "backscal": backscal,
                    "omega_arcmin2": omega,
                    "response_throughput": stacked_throughput(members),
                    "full_area_arcmin2": full_area,
                    "fit_low_keV": FIT_BAND_KEV[0],
                    "fit_high_keV": FIT_BAND_KEV[1],
                    "group_min_counts": PARAMS["group_min_counts"],
                    "channel_max": PARAMS["channel_max"][detector],
                    "pn_oot_fraction": SHARED["pn_oot_fraction"],
                    "source_counts": sum(row["source_counts"] for row in members),
                    "background_counts": sum(
                        row["background_counts"] for row in members
                    ),
                }
            )
    columns = (
        "region",
        "det",
        "n_frames",
        "exposure_s",
        "omega_arcmin2",
        "response_throughput",
        "full_area_arcmin2",
        "source_counts",
        "background_counts",
        "fit_low_keV",
        "fit_high_keV",
        "group_min_counts",
        "channel_max",
        "pn_oot_fraction",
        "source_pha",
        "background_pha",
        "rmf",
        "arf",
    )
    write_tsv(
        OUTPUT / "spectra.tsv",
        columns,
        [
            {
                name: (
                    Path(group[name]).name
                    if name in ("source_pha", "background_pha", "rmf", "arf")
                    else group[name]
                )
                for name in columns
            }
            for group in groups
        ],
    )
    return groups


def on_axis_factors(region_factors, groups):
    throughput = {
        (group["region"], group["det"]): group["response_throughput"]
        for group in groups
    }
    return {
        region: {
            detector: value / throughput[region, detector]
            for detector, value in factors.items()
        }
        for region, factors in region_factors.items()
    }


def image_factors(region_factors, groups):
    reference = region_factors["inner"]["mos2"]
    output = {}
    for region, factors in region_factors.items():
        members = [group for group in groups if group["region"] == region]
        weights = {
            group["det"]: group["exposure_s"]
            * group["omega_arcmin2"]
            * group["response_throughput"]
            for group in members
        }
        numerator = sum(weights[detector] * factors[detector] for detector in weights)
        denominator = sum(
            weights[detector] * region_factors["inner"][detector]
            for detector in weights
        )
        output[region] = reference * numerator / denominator
    return output


def fit_all(groups, frames):
    distance = np.average(
        [frame["position"]["distance_au"] for frame in frames],
        weights=[frame["position"]["duration_s"] for frame in frames],
    )
    heliocentric = np.average(
        [frame["position"]["r_au"] for frame in frames],
        weights=[frame["position"]["duration_s"] for frame in frames],
    )
    primary = fit.fit(
        groups,
        PARAMS["model"],
        distance,
        SOFTWARE["atomdb"],
        FIT_BAND_KEV,
        profile=True,
    )
    if primary["n_free_parameters"] != 9:
        raise RuntimeError("primary fit does not have nine free parameters")
    pn = fit.fit(
        [group for group in groups if group["det"] == "pn"],
        PARAMS["model"],
        distance,
        SOFTWARE["atomdb"],
        FIT_BAND_KEV,
    )
    inner = fit.fit(
        [group for group in groups if group["region"] == "inner"],
        PARAMS["model"],
        distance,
        SOFTWARE["atomdb"],
        FIT_BAND_KEV,
    )
    region_k = on_axis_factors(primary["region_detector_K"], groups)
    detector_k = region_k["inner"]
    image_k = image_factors(region_k, groups)
    regions = primary["regions"]
    mos2_sigma = ecf_sigma_ln(
        primary["crosscal"]["mos2"],
        primary["profile_intervals_1sigma"]["crosscal_mos2"],
    )
    for name, region in regions.items():
        region["mos2_equivalent_surface_brightness_ct_s_arcmin2"] = (
            region["surface_brightness_erg_cm2_s_arcmin2"] * image_k[name]
        )
        if name != "inner":
            region["image_response_ratio_to_inner"] = image_k[name] / detector_k["mos2"]
    result = {
        "distance_au": float(distance),
        "heliocentric_distance_au": float(heliocentric),
        "fit_band_keV": list(FIT_BAND_KEV),
        "model": {
            "expression": fit.EXPRESSION,
            "abundances": "wilm",
            "cross_sections": "vern",
            "recombtype": 1,
            "acxmodel": 8,
            "collntype": 1,
            "free_parameters": 9,
        },
        "fit": {
            "chi2": primary["statistic"],
            "dof": primary["dof"],
            "parameters": primary["parameters"],
            "crosscal": primary["crosscal"],
            "profile_intervals_1sigma": primary["profile_intervals_1sigma"],
        },
        "detector_K_ct_s_per_erg_cm2_s": detector_k,
        "response_throughput": {
            group["det"]: group["response_throughput"]
            for group in groups
            if group["region"] == "inner"
        },
        "mos2_ecf_erg_cm2_per_count": 1 / detector_k["mos2"],
        "mos2_ecf_sigma_ln": mos2_sigma,
        "mean_photon_energy_eV": primary["mean_photon_energy_eV"],
        "regions": regions,
        "checks": {
            "pn_only": {
                "chi2": pn["statistic"],
                "dof": pn["dof"],
                "inner_luminosity_erg_s": pn["regions"]["inner"]["luminosity_erg_s"],
            },
            "inner_only": {
                "chi2": inner["statistic"],
                "dof": inner["dof"],
                "inner_luminosity_erg_s": inner["regions"]["inner"]["luminosity_erg_s"],
            },
        },
    }
    folded = []
    for group, arrays in zip(primary["groups"], primary["arrays"]):
        for index in range(len(arrays["energy_keV"])):
            folded.append(
                {
                    "region": group["region"],
                    "det": group["det"],
                    **{name: values[index] for name, values in arrays.items()},
                }
            )
    columns = tuple(folded[0])
    write_tsv(OUTPUT / "folded.tsv", columns, folded)
    (OUTPUT / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    figure.draw(folded, OUTPUT / "figure2.png", OUTPUT / "figure2.pdf", FIT_BAND_KEV)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("extract", "fit"))
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)
    for name in FIT_PRODUCTS:
        (OUTPUT / name).unlink(missing_ok=True)
    if args.stage == "extract":
        for name in EXTRACTION_PRODUCTS:
            (OUTPUT / name).unlink(missing_ok=True)
        stack = OUTPUT / "stack"
        if stack.exists():
            if stack.resolve().parent != OUTPUT.resolve():
                raise RuntimeError("invalid Task 03 stack path")
            shutil.rmtree(stack)
    frames, sources = load_inputs()
    if args.stage == "extract":
        missing = [
            tool
            for tool in (
                "evselect",
                "backscale",
                "arfgen",
                "rmfgen",
                "addrmf",
                "addarf",
            )
            if shutil.which(tool) is None
        ]
        if missing:
            raise RuntimeError(f"source SAS/HEASoft first: {', '.join(missing)}")
        if WORK.exists():
            shutil.rmtree(WORK)
        stack_all(extract_all(frames, sources))
        return
    validate_extraction_frames(frames)
    groups = reload_groups()
    validate_fresh_inputs(frames, groups)
    fit_all(groups, frames)


if __name__ == "__main__":
    main()
