import math
import sys

import numpy as np

from _common import OUT, REPRO, read_json, stage_dir, write_json

sys.path.insert(0, str(REPRO))

import analysis_05_radial_profile.run as run5
from analysis_05_radial_profile import data as data5


def main():
    output = stage_dir(5)
    data5.IMAGE = OUT / "04" / "image.fits"
    run5.SPECTRUM = OUT / "03" / "result.json"
    run5.MORPHOLOGY = OUT / "04" / "result.json"
    run5.OUTPUT = output

    spectrum = read_json(run5.SPECTRUM)
    morphology = read_json(run5.MORPHOLOGY)
    run5.require_exposure_calibration(spectrum, morphology)
    if not np.allclose(
        np.asarray(spectrum["fit_band_keV"]) * 1000,
        run5.PARAMS["soft_band_ev"],
        rtol=0,
        atol=1e-9,
    ):
        raise RuntimeError("Task 03 fit band differs from the shared soft band")
    distance = float(spectrum["distance_au"])
    ecf = float(spectrum["mos2_ecf_erg_cm2_per_count"])
    width = float(run5.PARAMS["profile_bin_width_arcmin"])
    bins = round(float(run5.PARAMS["profile_fit_arcmin"][1]) / width)
    edges = np.arange(bins + 1) * width
    pixels = run5.epic_pixels()
    epic = run5.aggregate_epic(pixels, edges)
    first = int(
        round(
            run5.PARAMS["profile_fit_arcmin"][0]
            / run5.PARAMS["profile_bin_width_arcmin"]
        )
    )
    headline_sample = run5.subset(epic, first)
    headline = run5.fit_epic(headline_sample)
    broken = run5.fit_epic(epic, broken=True, background_interval=True)
    sector_slopes = [
        run5.fit_epic(run5.subset(run5.aggregate_epic(pixels, edges, sector), first))[
            "theta"
        ][1]
        for sector in range(12)
    ]
    sector_error = math.sqrt(
        11 / 12 * np.sum((np.asarray(sector_slopes) - np.mean(sector_slopes)) ** 2)
    )

    chemical_edges = np.asarray(run5.PARAMS["chemistry_annulus_edges_arcsec"]) / 60
    chemical = run5.aggregate_epic(pixels, chemical_edges)
    aperture = run5.aperture_rows(broken["theta"], distance, ecf)
    outer = run5.aggregate_epic(pixels, np.asarray(run5.SPECTRUM_OUTER) / 60)
    profile_sb = float(
        (outer["counts"][0] - outer["fixed"][0]) / outer["area"][0]
        - math.exp(broken["theta"][4])
    )
    spectral_sb = float(
        spectrum["regions"]["outer"]["mos2_equivalent_surface_brightness_ct_s_arcmin2"]
    )

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
    run5.write_table("profile.tsv", profile_columns, run5.epic_rows(epic, edges, broken))
    run5.write_table(
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
        "band_eV": run5.PARAMS["soft_band_ev"],
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
        "apertures": aperture_result,
        "profile_to_spectrum_370_ratio": (
            aperture_result[run5.SCIENCE_NAME]["luminosity_erg_s"]
            / spectrum["regions"]["inner"]["luminosity_erg_s"]
        ),
        run5.OUTER_NAME: {
            "profile_ct_s_arcmin2": profile_sb,
            "spectrum_ct_s_arcmin2": spectral_sb,
            "ratio": profile_sb / spectral_sb,
        },
    }
    write_json(output / "result.json", result)
    print(
        f"stage 05: slope {htheta[1]:.3f} +/- {sector_error:.3f} (sector),"
        f" L(370\") = {aperture_result[run5.SCIENCE_NAME]['luminosity_erg_s']:.3e} erg/s"
    )


if __name__ == "__main__":
    main()
