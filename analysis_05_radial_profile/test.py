import math

import numpy as np
from astropy.io import fits
from scipy.integrate import quad

from . import data as data_module
from .data import aggregate_epic, energy_masks, pps_maps
from .fit import (
    broken_shape,
    deviance,
    fit_epic,
    fit_pps,
    hard_expectation,
    integrated_broken,
    power_shape,
    radial_sum,
)
from .run import (
    MATCHED_APERTURES_KM,
    SCIENCE_APERTURE,
    SCIENCE_NAME,
    SPECTRUM_OUTER,
    FINAL_PRODUCTS,
    aperture_rows,
    require_exposure_calibration,
    subset,
)


def synthetic(alpha=1.12):
    radius = np.linspace(0.875, 11.875, 45)
    weight = np.full(45, 1e6)
    ring = np.arange(45)
    amplitude, sky = 0.012, 8e-5
    fixed = np.full(45, 300.0)
    counts = amplitude * weight * power_shape(radius, alpha) + sky * weight + fixed
    return {
        "counts": counts,
        "fixed": fixed,
        "area": weight,
        "radius": radius,
        "weight": weight,
        "ring": ring,
    }


def test_likelihood_and_fit():
    assert MATCHED_APERTURES_KM == (120000.0, 135000.0)
    sample = synthetic()
    fit = fit_epic(sample)
    assert deviance(sample["counts"], sample["counts"]) == 0
    assert abs(fit["theta"][1] - 1.12) < 1e-4
    assert abs(math.exp(fit["theta"][0]) / 0.012 - 1) < 1e-3
    assert abs(math.exp(fit["theta"][2]) / 8e-5 - 1) < 1e-3
    assert np.allclose(radial_sum(sample, np.ones(45)), sample["area"])


def test_profiled_background():
    soft = np.array([120.0, 20.0, 0.0])
    hard = np.array([80.0, 15.0, 0.0])
    vignetted = np.array([90.0, 12.0, 3.0])
    ratio = 0.2
    latent = hard_expectation(soft, hard, vignetted, ratio)
    residual = (
        ratio * (1 + ratio) * latent**2
        + (vignetted * (1 + ratio) - ratio * (soft + hard)) * latent
        - hard * vignetted
    )
    assert np.all(np.isfinite(latent)) and np.all(latent > 0)
    assert np.allclose(residual, 0, atol=1e-10)

    radius = np.linspace(0.875, 11.875, 45)
    area = np.full(45, 2e5)
    vignetted = 0.02 * area * power_shape(radius, 1.05) + 1e-4 * area
    latent = np.linspace(300, 900, 45)
    latent[0] = 0
    ratio = 0.12
    sample = {
        "soft": vignetted + ratio * latent,
        "hard": latent,
        "area": area,
        "radius": radius,
        "weight": area,
        "ring": np.arange(45),
    }
    fit = fit_pps(sample)
    assert abs(fit["theta"][1] - 1.05) < 0.01
    assert abs(fit["soft_to_hard"] - ratio) < 0.01


def test_broken_integral():
    amplitude, inner, outer, radius_break, radius = 0.02, 0.35, 1.1, 0.63, 6.2
    exact = integrated_broken(amplitude, inner, outer, radius_break, radius)
    numeric = (
        2
        * math.pi
        * quad(
            lambda r: amplitude
            * broken_shape(np.array([r]), inner, outer, radius_break)[0]
            * r,
            0,
            radius,
        )[0]
    )
    assert abs(exact / numeric - 1) < 2e-4


def test_shared_apertures():
    rows = aperture_rows(
        np.array([math.log(0.02), 0.3, 1.1, math.log(0.6)]), 1.8, 1e-11
    )
    assert rows[-1][:2] == [SCIENCE_NAME, SCIENCE_APERTURE]
    assert SCIENCE_APERTURE == 370 and SPECTRUM_OUTER == (400, 600)


def test_band_boundaries():
    values = np.array([349, 350, 1100, 1101, 2499, 2500, 8500, 8501])
    soft, hard = energy_masks(values)
    assert np.array_equal(soft, [False, True, True, False, False, False, False, False])
    assert np.array_equal(hard, [False, False, False, False, False, True, True, False])


def test_aggregation_and_subset():
    one = np.ones((2, 2))
    planes = {
        "EVENTS": 5 * one,
        "QPB": one,
        "SP": 2 * one,
        "OOT": 3 * one,
        "EXPOSURE": 10 * one,
        "radius": np.array([[0.1, 0.2], [0.6, 0.7]]),
        "sector": np.array([[0, 1], [0, 1]]),
        "pixel_area_arcmin2": 0.2,
    }
    sample = aggregate_epic(planes, np.array([0, 0.5, 1]))
    assert np.all(sample["counts"] == 10)
    assert np.all(sample["area"] == 4)
    assert np.allclose(sample["fixed"], 2 + 4 + 6 * 0.063)
    trimmed = subset(sample, 1)
    assert len(trimmed["counts"]) == 1 and np.all(trimmed["ring"] == 0)


def test_exposure_masking():
    planes = {
        "EVENTS": np.array([[5.0, 1000.0], [7.0, 11.0]]),
        "QPB": np.ones((2, 2)),
        "SP": 2 * np.ones((2, 2)),
        "OOT": 3 * np.ones((2, 2)),
        "EXPOSURE": np.array([[2.0, 0.0], [6.0, 4.0]]),
        "DISPLAY_MASK": np.zeros((2, 2)),
        "radius": np.array([[0.1, 0.2], [0.3, 0.4]]),
        "sector": np.zeros((2, 2), int),
        "pixel_area_arcmin2": 0.25,
    }
    sample = aggregate_epic(planes, np.array([0.0, 0.5]))
    assert sample["counts"][0] == 23
    assert sample["area"][0] == 3
    assert sample["coverage"][0] == 0.75
    assert np.isclose(sample["fixed"][0], 3 + 6 + 9 * 0.063)
    shape = power_shape(sample["radius"], 1.0)
    assert np.isclose(
        radial_sum(sample, shape), np.sum(sample["weight"] / sample["radius"])
    )


def test_exposure_calibration():
    spectrum = {"detector_K_ct_s_per_erg_cm2_s": {"pn": 4.0, "mos1": 2.0, "mos2": 1.0}}
    morphology = {"detector_exposure_weights": {"pn": 4.0, "mos1": 2.0, "mos2": 1.0}}
    require_exposure_calibration(spectrum, morphology)
    morphology["detector_exposure_weights"]["pn"] = 3.9
    try:
        require_exposure_calibration(spectrum, morphology)
    except RuntimeError as error:
        assert "different Task 03 calibration" in str(error)
    else:
        raise AssertionError("stale Task 04 calibration accepted")


def test_pps_map_grid():
    header = fits.Header()
    header["CTYPE1"], header["CTYPE2"] = "RA---TAN", "DEC--TAN"
    header["CUNIT1"], header["CUNIT2"] = "deg", "deg"
    header["CRPIX1"] = data_module.pps.MAP_CENTER[0] + 1
    header["CRPIX2"] = data_module.pps.MAP_CENTER[1] + 1
    header["CRVAL1"], header["CRVAL2"] = data_module.pps.REFERENCE
    header["CDELT1"] = -data_module.pps.MAP_PIXEL / 3600
    header["CDELT2"] = data_module.pps.MAP_PIXEL / 3600
    header["RADESYS"] = "FK5"

    def load(headers):
        rows = iter((np.ones((5, 6)), item.copy()) for item in headers)
        original = data_module.fits.getdata
        data_module.fits.getdata = lambda *args, **kwargs: next(rows)
        try:
            return pps_maps()
        finally:
            data_module.fits.getdata = original

    assert len(load([header] * 3)) == 3
    shifted = header.copy()
    shifted["CRPIX1"] += 1
    try:
        load([header, shifted, header])
    except RuntimeError as error:
        assert "PPS exposure" in str(error)
    else:
        raise AssertionError("misregistered PPS map accepted")
    transformed = header.copy()
    transformed["PC1_2"] = 0.01
    try:
        load([header, transformed, header])
    except RuntimeError as error:
        assert "PPS exposure" in str(error)
    else:
        raise AssertionError("PPS map with a different linear WCS accepted")
    try:
        load([transformed] * 3)
    except RuntimeError as error:
        assert "PPS exposure" in str(error)
    else:
        raise AssertionError("noncanonical PPS WCS accepted")


def test_output_contract():
    assert FINAL_PRODUCTS == (
        "profile.tsv",
        "annuli.tsv",
        "pps_profile.tsv",
        "apertures.tsv",
        "result.json",
    )


def main():
    test_likelihood_and_fit()
    test_profiled_background()
    test_broken_integral()
    test_shared_apertures()
    test_band_boundaries()
    test_aggregation_and_subset()
    test_exposure_masking()
    test_exposure_calibration()
    test_pps_map_grid()
    test_output_contract()


if __name__ == "__main__":
    main()
