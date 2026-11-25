import csv
import json
import math
from pathlib import Path

import numpy as np

from .data import PARAMS, aggregate_epic, epic_pixels
from .fit import (
    sample_shape,
    fit_epic,
    integrated_broken,
    radial_sum,
)

HERE = Path(__file__).resolve().parent
REPRO = HERE.parent
SPECTRUM = REPRO / "analysis_03_spectrum/outputs/result.json"
PPS = REPRO / "analysis_02_pps_luminosity/outputs"
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
    shape = sample_shape(sample, theta[1], theta[2], math.exp(theta[3]))
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


def pps_rows():
    band = "-".join(map(str, PARAMS["soft_band_ev"]))
    result = json.loads((PPS / "result.json").read_text())
    if result["fit_range_arcmin"] != PARAMS["profile_fit_arcmin"]:
        raise RuntimeError("PPS profile has a different fit range")
    fit = result["bands"][band]
    with (PPS / "profiles.tsv").open() as stream:
        rows = [
            row
            for row in csv.DictReader(stream, delimiter="\t")
            if row["band_eV"] == band
        ]
    output = []
    for row in rows:
        lo, hi = float(row["r_lo_arcmin"]), float(row["r_hi_arcmin"])
        soft, hard = float(row["events_count"]), float(row["hard_count"])
        area = float(row["exposure_s_arcmin2"])
        output.append(
            [
                lo,
                hi,
                centers(np.array([lo, hi]))[0],
                soft,
                hard,
                area,
                (soft - float(row["particle_count"])) / area,
                math.sqrt(soft + fit["soft_to_hard"] ** 2 * hard) / area,
                float(row["comet_count"]) / area,
                fit["sky_ct_s_arcmin2"],
            ]
        )
    return fit, output


def fit_profile(pixels, edges):
    sample = aggregate_epic(pixels, edges)
    fit = fit_epic(sample, background_interval=True)
    sectors = np.array(
        [
            fit_epic(aggregate_epic(pixels, edges, sector))["theta"][1:3]
            for sector in range(12)
        ]
    )
    errors = np.sqrt(11 / 12 * np.sum((sectors - sectors.mean(axis=0)) ** 2, axis=0))
    theta = fit["theta"]
    broken = {
        "normalization_at_break_ct_s_arcmin2": math.exp(theta[0]),
        "alpha_inner": float(theta[1]),
        "alpha_outer": float(theta[2]),
        "break_radius_arcmin": math.exp(theta[3]),
        "poisson_deviance": fit["stat"],
        "degrees_of_freedom": len(sample["counts"]) - 5,
    }
    result = {
        "model": "broken_power_law",
        "fit_range_arcmin": [float(edges[0]), float(edges[-1])],
        "slope": float(theta[2]),
        "normalization_ct_s_arcmin2": math.exp(theta[0] + theta[2] * theta[3]),
        "sky_ct_s_arcmin2": math.exp(theta[4]),
        "poisson_deviance": fit["stat"],
        "degrees_of_freedom": len(sample["counts"]) - 5,
        "sector_standard_error": float(errors[1]),
        "sector_standard_errors": {
            "alpha_inner": float(errors[0]),
            "alpha_outer": float(errors[1]),
        },
        "broken": broken,
        "background_prior": {
            "mean_ct_s_arcmin2": math.exp(theta[4]),
            "sigma_ct_s_arcmin2": fit["background_error"],
        },
    }
    return sample, fit, result


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
    if PARAMS["profile_fit_arcmin"][0] != 0:
        raise ValueError("broken profile must include the center")
    epic, broken, epic_result = fit_profile(pixels, edges)
    pps_result, pps_values = pps_rows()

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
        pps_values,
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
        "epic": epic_result,
        "pps": pps_result,
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
    from .figure import make_figure

    make_figure(OUTPUT, epic_result, pps_result)
    print(
        f"EPIC alpha={epic_result['broken']['alpha_inner']:.4f}, "
        f"{epic_result['slope']:.4f} +/- {epic_result['sector_standard_error']:.4f}; "
        f"break={epic_result['broken']['break_radius_arcmin']:.4f} arcmin"
    )


if __name__ == "__main__":
    main()
