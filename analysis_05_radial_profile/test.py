import math

import numpy as np
from scipy.integrate import quad

from .data import aggregate_epic
from .fit import (
    broken_shape,
    deviance,
    fit_epic,
    fit_pps,
    hard_expectation,
    integrated_broken,
    sample_shape,
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
)


def synthetic():
    edges = np.linspace(0, 12, 49)
    radius = np.linspace(0.005, 11.995, 2400)
    ring = np.searchsorted(edges, radius, side="right") - 1
    weight = 2 * np.pi * radius * 0.005 * 3e4 * np.exp(-radius / 30)
    area = np.bincount(ring, weight)
    sample = {
        "radius": radius,
        "weight": weight,
        "ring": ring,
        "area": area,
        "counts": np.zeros(48),
        "fixed": np.full(48, 300.0),
    }
    sample["counts"] = (
        0.02 * radial_sum(sample, broken_shape(radius, 0.4, 1.1, 0.74))
        + 8e-5 * area
        + sample["fixed"]
    )
    return sample


def test_likelihood_and_fit():
    assert MATCHED_APERTURES_KM == (120000.0, 135000.0)
    sample = synthetic()
    fit = fit_epic(sample)
    assert deviance(sample["counts"], sample["counts"]) == 0
    assert np.allclose(fit["theta"][1:3], [0.4, 1.1], atol=2e-4)
    assert abs(math.exp(fit["theta"][3]) / 0.74 - 1) < 2e-3
    assert abs(math.exp(fit["theta"][0]) / 0.02 - 1) < 2e-3
    assert abs(math.exp(fit["theta"][4]) / 8e-5 - 1) < 2e-3
    assert np.allclose(
        radial_sum(sample, np.ones(len(sample["radius"]))), sample["area"]
    )


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
    sample = synthetic()
    hard = 100 + 500 * (1 + np.sin(np.arange(48) * 0.7))
    hard[0] = 0
    sample["soft"] = sample["counts"] - sample["fixed"] + 0.12 * hard
    sample["hard"] = hard
    fit = fit_pps(sample)
    assert np.allclose(fit["theta"][1:3], [0.4, 1.1], atol=2e-3)
    assert abs(math.exp(fit["theta"][3]) / 0.74 - 1) < 3e-3
    assert abs(fit["soft_to_hard"] - 0.12) < 2e-3
    assert fit["stat"] < 1e-3


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


def test_center_without_core():
    # Well inside 2 arcsec, neither the intensity nor the aperture flux flattens.
    radius = np.array([0.1, 1.0, 2.0]) / 60
    assert np.allclose(broken_shape(radius, 0.5, 1.1, 0.7), (radius / 0.7) ** -0.5)
    for r in radius:
        flux = integrated_broken(0.02, 0.5, 1.1, 0.7, r)
        assert np.isclose(flux, 2 * math.pi * 0.02 * 0.7**0.5 * r**1.5 / 1.5)
    assert integrated_broken(1, 1.9, 2, 0.7, 0) == 0
    # The exact square-pixel mean of r^-1 is 2*asinh(1)/half_width.
    half_width = 2 / 60
    sample = {"radius": np.array([0.0, 0.2]), "pixel_area_arcmin2": (2 * half_width)**2}
    shape = sample_shape(sample, 1, 1.1, 0.7)
    assert np.isclose(shape[0], 2 * math.asinh(1) * 0.7 / half_width, rtol=1e-12)
    assert np.isclose(shape[1], 0.7 / 0.2)
    assert np.all(sample_shape(sample, 0, 0, 0.7) > 0)
    assert np.allclose(sample_shape(sample, 0, 0, 0.7), 1)


def test_shared_apertures():
    rows = aperture_rows(
        np.array([math.log(0.02), 0.3, 1.1, math.log(0.6)]), 1.8, 1e-11
    )
    assert rows[-1][:2] == [SCIENCE_NAME, SCIENCE_APERTURE]
    assert SCIENCE_APERTURE == 370 and SPECTRUM_OUTER == (400, 600)


def test_aggregation():
    one = np.ones((2, 2))
    planes = {
        "EVENTS": 5 * one,
        "QPB": one,
        "SP": 2 * one,
        "OOT": 3 * one,
        "EXPOSURE": 10 * one,
        "radius": np.array([[0.0, 0.2], [0.6, 0.7]]),
        "sector": np.array([[0, 1], [0, 1]]),
        "pixel_area_arcmin2": 0.2,
    }
    sample = aggregate_epic(planes, np.array([0, 0.5, 1]))
    assert np.all(sample["counts"] == 10)
    assert np.all(sample["area"] == 4)
    assert np.allclose(sample["fixed"], 2 + 4 + 6 * 0.063)
    assert sample["radius"][0] == 0
    assert np.all(np.isfinite(sample_shape(sample, 0.4, 1.1, 0.7)))


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
    shape = broken_shape(sample["radius"], 1.0, 1.0, 1.0)
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
    test_center_without_core()
    test_shared_apertures()
    test_aggregation()
    test_exposure_masking()
    test_exposure_calibration()
    test_output_contract()


if __name__ == "__main__":
    main()
