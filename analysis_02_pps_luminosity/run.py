import csv
import hashlib
import json
import math

import numpy as np

from analysis_05_radial_profile.fit import (
    sample_shape,
    fit_pps,
    integrated_broken,
    radial_sum,
)
from .data import Dataset, HERE, NATIVE, SHARED, pattern_max

PARAMS = json.loads((HERE / "parameters.json").read_text())
OUTPUT = HERE / "outputs"
AU_CM = 1.495978707e13


def edges():
    lo, hi = SHARED["profile_fit_arcmin"]
    if lo != 0:
        raise ValueError("broken profiles must include the center")
    return (
        np.arange(round(hi / SHARED["profile_bin_width_arcmin"]) + 1)
        * SHARED["profile_bin_width_arcmin"]
    )


def profile_result(fit, sample):
    theta = fit["theta"]
    return {
        "model": "broken_power_law",
        "alpha_inner": float(theta[1]),
        "alpha_outer": float(theta[2]),
        "break_radius_arcmin": math.exp(theta[3]),
        "normalization_at_break_ct_s_arcmin2": math.exp(theta[0]),
        "sky_ct_s_arcmin2": math.exp(theta[4]),
        "soft_to_hard": fit["soft_to_hard"],
        "poisson_deviance": fit["stat"],
        "degrees_of_freedom": len(sample["soft"]) - 6,
    }


def apertures():
    distance = PARAMS["distance_au"]
    requests = [
        (
            f"{radius:g}_km",
            math.degrees(math.atan(radius / (distance * AU_CM / 1e5))) * 3600,
        )
        for radius in PARAMS["apertures_km"]
    ]
    return requests + [
        (f"{radius:g}_arcsec", radius) for radius in PARAMS["apertures_arcsec"]
    ]


def luminosities(fit, ecf):
    theta = fit["theta"]
    factor = 4 * math.pi * (PARAMS["distance_au"] * AU_CM) ** 2
    result = {}
    for name, radius in apertures():
        rate = integrated_broken(
            math.exp(theta[0]), theta[1], theta[2], math.exp(theta[3]), radius / 60
        )
        result[name] = {
            "radius_arcsec": radius,
            "model_rate_ct_s": rate,
            "energy_flux_erg_cm2_s": rate / ecf,
            "luminosity_erg_s": factor * rate / ecf,
        }
    return result


def calculate(samples, ecfs):
    bands, profiles = {}, []
    for band, sample in samples.items():
        fit = fit_pps(sample)
        low, high = band
        name = f"{low}-{high}"
        result = profile_result(fit, sample)
        result["energy_eV"] = list(band)
        result["pattern_max"] = pattern_max(band)
        result["exposure_weights"] = sample["exposure_weights"].tolist()
        if band in ecfs:
            result["ecf_counts_cm2_per_erg"] = ecfs[band]
            result["apertures"] = luminosities(fit, ecfs[band])
        bands[name] = result
        theta = fit["theta"]
        shape = sample_shape(sample, theta[1], theta[2], math.exp(theta[3]))
        comet = math.exp(theta[0]) * radial_sum(sample, shape)
        sky = math.exp(theta[4]) * sample["area"]
        particle = fit["soft_to_hard"] * fit["hard_expectation"]
        for i, (lo, hi) in enumerate(zip(edges()[:-1], edges()[1:])):
            profiles.append(
                [
                    name,
                    lo,
                    hi,
                    sample["soft"][i],
                    sample["hard"][i],
                    sample["area"][i],
                    comet[i],
                    sky[i],
                    particle[i],
                    fit["hard_expectation"][i],
                ]
            )
        print(
            f"PPS {name} eV: alpha={theta[1]:.4f}, {theta[2]:.4f}; break={math.exp(theta[3]):.4f} arcmin",
            flush=True,
        )
    sums = {}
    for low, high in PARAMS["sum_bands_eV"]:
        parts = [
            b
            for b in bands.values()
            if tuple(b["energy_eV"]) in NATIVE
            and "apertures" in b
            and low <= b["energy_eV"][0]
            and b["energy_eV"][1] <= high
        ]
        limits = sorted(b["energy_eV"] for b in parts)
        if (
            not limits
            or limits[0][0] != low
            or limits[-1][1] != high
            or any(a[1] != b[0] for a, b in zip(limits[:-1], limits[1:]))
        ):
            raise ValueError("summed bands must tile the requested photon band")
        sums[f"{low}-{high}"] = {
            name: {
                "radius_arcsec": radius,
                "energy_flux_erg_cm2_s": sum(
                    b["apertures"][name]["energy_flux_erg_cm2_s"] for b in parts
                ),
                "luminosity_erg_s": sum(
                    b["apertures"][name]["luminosity_erg_s"] for b in parts
                ),
            }
            for name, radius in apertures()
        }
    return {"bands": bands, "summed_bands": sums}, profiles


def write_outputs(result, profiles):
    OUTPUT.mkdir(exist_ok=True)
    (OUTPUT / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    with (OUTPUT / "profiles.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(
            [
                "band_eV",
                "r_lo_arcmin",
                "r_hi_arcmin",
                "events_count",
                "hard_count",
                "exposure_s_arcmin2",
                "comet_count",
                "sky_count",
                "particle_count",
                "hard_model_count",
            ]
        )
        writer.writerows(profiles)
    with (OUTPUT / "luminosities.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(
            [
                "band_eV",
                "aperture",
                "radius_arcsec",
                "energy_flux_erg_cm2_s",
                "luminosity_erg_s",
            ]
        )
        groups = {
            k: v["apertures"] for k, v in result["bands"].items() if "apertures" in v
        }
        groups.update(result["summed_bands"])
        for band, values in groups.items():
            for name, row in values.items():
                writer.writerow(
                    [
                        band,
                        name,
                        row["radius_arcsec"],
                        row["energy_flux_erg_cm2_s"],
                        row["luminosity_erg_s"],
                    ]
                )


def main():
    data = Dataset(edges())
    requests = [tuple(band) for band in PARAMS["bands_eV"]] + [
        tuple(SHARED["soft_band_ev"])
    ]
    samples = {
        band: data.band(*band, singles=pattern_max(band) == 0) for band in requests
    }
    calibration = json.loads((HERE / "ecfs.json").read_text())
    ecfs = {}
    for band in requests:
        row = calibration["bands"][f"{band[0]}-{band[1]}"]
        if row["pattern_max"] != pattern_max(band):
            raise ValueError("standard ECF must match event pattern selection")
        ecfs[band] = row["ecf_counts_cm2_per_erg"]
    result, profiles = calculate(samples, ecfs)
    for name, row in result["bands"].items():
        row["ecf_source"] = "fixed spectrum folded through standard response; ecfs.json"
    result.update(
        epoch_utc=data.epoch,
        distance_au=PARAMS["distance_au"],
        good_time_s=data.good_time,
        fit_range_arcmin=SHARED["profile_fit_arcmin"],
        selection="PPS CCD GTIs; RAWY>12; FLAG==0; PI in (low, high]; no field-source masks",
        calibration={
            "method": "direct standard-response ECFs in every band",
            "photon_index": 1.7,
            "nh_cm-2": 3e20,
            "response_conversion": calibration,
            "exposure_approximation": "reuse native PPS maps; sub-bands below 500 eV use EXPMAP1000",
        },
        inputs=[
            {
                "path": str(path),
                "sha256": hashlib.file_digest(path.open("rb"), "sha256").hexdigest(),
            }
            for path in data.inputs + [HERE / "ecfs.json"]
        ],
    )
    write_outputs(result, profiles)
    # Compact sufficient statistics let the demo rerun these fits without the event list.
    packed = {
        f"{lo}_{hi}_{key}": value
        for (lo, hi), sample in samples.items()
        for key, value in sample.items()
    }
    np.savez_compressed(OUTPUT / "samples.npz", **packed)


if __name__ == "__main__":
    main()
