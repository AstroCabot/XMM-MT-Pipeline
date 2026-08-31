import csv
import json
import math
from pathlib import Path

import numpy as np

from .data import PARAMS, aggregate_epic, epic_pixels, pps_profile
from .fit import (
    broken_shape,
    fit_epic,
    fit_pps,
    hard_expectation,
    integrated_broken,
    power_shape,
    radial_sum,
)

HERE = Path(__file__).resolve().parent
REPRO = HERE.parent
SPECTRUM = REPRO / "analysis_03_spectrum/outputs/result.json"
MORPHOLOGY = REPRO / "analysis_04_morphology/outputs/result.json"
OUTPUT = HERE / "outputs"
FINAL_PRODUCTS = (
    "profile.tsv",
    "annuli.tsv",
    "pps_profile.tsv",
    "apertures.tsv",
    "result.json",
)
AU_CM = 1.495978707e13
AU_KM = 1.495978707e8
SCIENCE_APERTURE = float(PARAMS["science_aperture_radius_arcsec"])
SPECTRUM_OUTER = tuple(map(float, PARAMS["spectrum_outer_annulus_arcsec"]))
MATCHED_APERTURES_KM = tuple(map(float, PARAMS["matched_apertures_km"]))
SCIENCE_NAME = f"{SCIENCE_APERTURE:g}_arcsec"
OUTER_NAME = f"outer_{SPECTRUM_OUTER[0]:g}_{SPECTRUM_OUTER[1]:g}"


def require_exposure_calibration(spectrum, morphology):
    conversion = spectrum["detector_K_ct_s_per_erg_cm2_s"]
    reference = float(conversion["mos2"])
    expected = {
        name: float(conversion[name]) / reference for name in ("pn", "mos1", "mos2")
    }
    actual = morphology.get("detector_exposure_weights", {})
    if set(actual) != set(expected) or not np.allclose(
        [actual[name] for name in expected], list(expected.values()), rtol=1e-12, atol=0
    ):
        raise RuntimeError("Task 04 exposure uses a different Task 03 calibration")


def subset(sample, first):
    n = len(sample["counts"] if "counts" in sample else sample["soft"])
    keep = sample["ring"] >= first
    output = {
        key: value[first:]
        for key, value in sample.items()
        if isinstance(value, np.ndarray) and value.shape == (n,)
    }
    for key in ("radius", "weight"):
        output[key] = sample[key][keep]
    output["ring"] = sample["ring"][keep] - first
    if np.any(output["area"] <= 0):
        raise ValueError("profile contains a zero-exposure ring")
    return output


def centers(edges):
    lo, hi = edges[:-1], edges[1:]
    return (2 / 3) * (hi**3 - lo**3) / (hi**2 - lo**2)


def write_table(name, columns, rows):
    with (OUTPUT / name).open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(columns)
        writer.writerows(rows)


def epic_rows(sample, edges, fit):
    theta = fit["theta"]
    shape = broken_shape(sample["radius"], theta[1], theta[2], math.exp(theta[3]))
    comet = math.exp(theta[0]) * radial_sum(sample, shape) / sample["area"]
    measured = (sample["counts"] - sample["fixed"]) / sample["area"]
    error = np.sqrt(sample["variance"]) / sample["area"]
    values = zip(
        edges[:-1],
        edges[1:],
        centers(edges),
        sample["counts"],
        sample["qpb"],
        sample["sp"],
        sample["oot"],
        sample["area"],
        sample["coverage"],
        measured,
        error,
        comet,
        strict=True,
    )
    return [list(row) + [math.exp(theta[4])] for row in values]


def pps_rows(sample, edges, fit):
    theta, ratio = fit["theta"], fit["soft_to_hard"]
    comet = (
        math.exp(theta[0])
        * radial_sum(sample, power_shape(sample["radius"], theta[1]))
        / sample["area"]
    )
    particle = ratio * fit["hard_expectation"]
    measured = (sample["soft"] - particle) / sample["area"]
    error = np.sqrt(sample["soft"] + ratio**2 * sample["hard"]) / sample["area"]
    return [
        list(row) + [math.exp(theta[2])]
        for row in zip(
            edges[:-1],
            edges[1:],
            centers(edges),
            sample["soft"],
            sample["hard"],
            sample["area"],
            measured,
            error,
            comet,
            strict=True,
        )
    ]


def aperture_rows(theta, distance, ecf):
    amplitude, radius_break = math.exp(theta[0]), math.exp(theta[3])
    requests = [
        (f"{radius:g}_km", math.degrees(math.atan(radius / (distance * AU_KM))) * 3600)
        for radius in MATCHED_APERTURES_KM
    ]
    requests.append((SCIENCE_NAME, SCIENCE_APERTURE))
    factor = 4 * math.pi * (distance * AU_CM) ** 2
    rows = []
    for name, arcsec in requests:
        rate = integrated_broken(
            amplitude, theta[1], theta[2], radius_break, arcsec / 60
        )
        flux = rate * ecf
        rows.append(
            [
                name,
                arcsec,
                distance * AU_KM * math.tan(math.radians(arcsec / 3600)),
                rate,
                flux,
                flux * factor,
            ]
        )
    return rows


def main():
    for name in FINAL_PRODUCTS:
        (OUTPUT / name).unlink(missing_ok=True)
    spectrum = json.loads(SPECTRUM.read_text())
    morphology = json.loads(MORPHOLOGY.read_text())
    require_exposure_calibration(spectrum, morphology)
    if not np.allclose(
        np.asarray(spectrum["fit_band_keV"]) * 1000,
        PARAMS["soft_band_ev"],
        rtol=0,
        atol=1e-9,
    ):
        raise RuntimeError("Task 03 fit band differs from the shared soft band")
    distance = float(spectrum["distance_au"])
    ecf = float(spectrum["mos2_ecf_erg_cm2_per_count"])
    width = float(PARAMS["profile_bin_width_arcmin"])
    bins = round(float(PARAMS["profile_fit_arcmin"][1]) / width)
    edges = np.arange(bins + 1) * width
    pixels = epic_pixels()
    epic = aggregate_epic(pixels, edges)
    first = int(
        round(PARAMS["profile_fit_arcmin"][0] / PARAMS["profile_bin_width_arcmin"])
    )
    headline_sample = subset(epic, first)
    headline = fit_epic(headline_sample)
    broken = fit_epic(epic, broken=True, background_interval=True)
    sector_slopes = [
        fit_epic(subset(aggregate_epic(pixels, edges, sector), first))["theta"][1]
        for sector in range(12)
    ]
    sector_error = math.sqrt(
        11 / 12 * np.sum((np.asarray(sector_slopes) - np.mean(sector_slopes)) ** 2)
    )

    pps = pps_profile(edges)
    pps_fit = fit_pps(subset(pps, first))
    pps_full = {**pps_fit}
    theta = pps_fit["theta"]
    vignetted = (
        math.exp(theta[0]) * radial_sum(pps, power_shape(pps["radius"], theta[1]))
        + math.exp(theta[2]) * pps["area"]
    )
    pps_full["hard_expectation"] = hard_expectation(
        pps["soft"], pps["hard"], vignetted, pps_fit["soft_to_hard"]
    )

    chemical_edges = np.asarray(PARAMS["chemistry_annulus_edges_arcsec"]) / 60
    chemical = aggregate_epic(pixels, chemical_edges)
    aperture = aperture_rows(broken["theta"], distance, ecf)
    outer = aggregate_epic(pixels, np.asarray(SPECTRUM_OUTER) / 60)
    profile_sb = float(
        (outer["counts"][0] - outer["fixed"][0]) / outer["area"][0]
        - math.exp(broken["theta"][4])
    )
    spectral_sb = float(
        spectrum["regions"]["outer"]["mos2_equivalent_surface_brightness_ct_s_arcmin2"]
    )

    OUTPUT.mkdir(exist_ok=True)
    profile_columns = (
        "r_lo_arcmin",
        "r_hi_arcmin",
        "r_arcmin",
        "events_count",
        "qpb_count",
        "sp_count",
        "oot_count",
        "exposure_s_arcmin2",
        "coverage",
        "net_ct_s_arcmin2",
        "error_ct_s_arcmin2",
        "comet_model_ct_s_arcmin2",
        "sky_ct_s_arcmin2",
    )
    write_table("profile.tsv", profile_columns, epic_rows(epic, edges, broken))
    write_table(
        "annuli.tsv",
        (
            "r_lo_arcmin",
            "r_hi_arcmin",
            "events_count",
            "qpb_count",
            "sp_count",
            "oot_count",
            "exposure_s_arcmin2",
        ),
        zip(
            chemical_edges[:-1],
            chemical_edges[1:],
            chemical["counts"],
            chemical["qpb"],
            chemical["sp"],
            chemical["oot"],
            chemical["area"],
            strict=True,
        ),
    )
    write_table(
        "pps_profile.tsv",
        (
            "r_lo_arcmin",
            "r_hi_arcmin",
            "r_arcmin",
            "soft_count",
            "hard_count",
            "exposure_s_arcmin2",
            "net_ct_s_arcmin2",
            "error_ct_s_arcmin2",
            "comet_model_ct_s_arcmin2",
            "sky_ct_s_arcmin2",
        ),
        pps_rows(pps, edges, pps_full),
    )
    write_table(
        "apertures.tsv",
        (
            "name",
            "radius_arcsec",
            "radius_km",
            "model_rate_ct_s",
            "energy_flux_erg_cm2_s",
            "luminosity_erg_s",
        ),
        aperture,
    )
    btheta, htheta = broken["theta"], headline["theta"]
    aperture_result = {
        row[0]: {
            "radius_arcsec": row[1],
            "radius_km": row[2],
            "model_rate_ct_s": row[3],
            "energy_flux_erg_cm2_s": row[4],
            "luminosity_erg_s": row[5],
        }
        for row in aperture
    }
    result = {
        "band_eV": PARAMS["soft_band_ev"],
        "epic": {
            "slope": float(htheta[1]),
            "normalization_ct_s_arcmin2": math.exp(htheta[0]),
            "sky_ct_s_arcmin2": math.exp(htheta[2]),
            "poisson_deviance": headline["stat"],
            "degrees_of_freedom": len(headline_sample["counts"]) - 3,
            "sector_standard_error": sector_error,
            "broken": {
                "normalization_at_break_ct_s_arcmin2": math.exp(btheta[0]),
                "alpha_inner": float(btheta[1]),
                "alpha_outer": float(btheta[2]),
                "break_radius_arcmin": math.exp(btheta[3]),
                "poisson_deviance": broken["stat"],
                "degrees_of_freedom": len(epic["counts"]) - 5,
            },
            "background_prior": {
                "mean_ct_s_arcmin2": math.exp(btheta[4]),
                "sigma_ct_s_arcmin2": broken["background_error"],
            },
        },
        "pps": {"slope": float(theta[1])},
        "apertures": aperture_result,
        "profile_to_spectrum_370_ratio": (
            aperture_result[SCIENCE_NAME]["luminosity_erg_s"]
            / spectrum["regions"]["inner"]["luminosity_erg_s"]
        ),
        OUTER_NAME: {
            "profile_ct_s_arcmin2": profile_sb,
            "spectrum_ct_s_arcmin2": spectral_sb,
            "ratio": profile_sb / spectral_sb,
        },
    }
    (OUTPUT / "result.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
