import math

import numpy as np

from .calibrate import channel_weights
from .data import Dataset, NATIVE, gti_mask, pattern_max, selected
from .run import AU_CM, PARAMS, edges, luminosities


def test_selection():
    pi = np.array([200, 201, 500, 501, 1000, 1001, 2000, 2001])
    masks = [selected(pi, lo, hi) for lo, hi in NATIVE]
    assert np.array_equal(np.sum(masks, axis=0), [0, 1, 1, 1, 1, 1, 1, 0])
    events = np.array(
        [(1, 0.0), (1, 1.0), (1, 2.0), (2, 1.0), (2, 3.0)],
        dtype=[("CCDNR", "i4"), ("TIME", "f8")],
    )
    mask = gti_mask(events, [np.array([[0.0, 1.0]]), np.array([[2.0, 3.0]])])
    assert np.array_equal(mask, [True, True, False, False, True])


def test_native_exposure_and_patterns():
    data = Dataset.__new__(Dataset)
    data.edges = np.array([0.0, 1.0, 2.0])
    data.pi = np.array([250, 300, 500, 501, 900, 1000, 1200])
    data.pattern = np.array([0, 1, 0, 0, 1, 0, 0])
    data.radius = np.full(7, 0.5)
    data.hard = np.array([3, 4])
    data.map_radius = np.array([[0.25, 0.75], [1.25, 1.75]])
    data.maps = [np.full((2, 2), value) for value in (10.0, 20.0, 30.0)]
    data.pixel_area = 0.5
    soft = data.band(200, 500, singles=True)
    assert soft["soft"].sum() == 2
    assert np.array_equal(soft["exposure_weights"], [1, 0, 0])
    assert np.allclose(soft["area"], 10)
    harder = data.band(500, 1000)
    assert harder["soft"].sum() == 3
    assert np.array_equal(harder["exposure_weights"], [0, 1, 0])
    assert np.allclose(harder["area"], 20)
    narrow = data.band(200, 350, singles=pattern_max((200, 350)) == 0)
    assert narrow["soft"].sum() == 1
    assert np.array_equal(narrow["exposure_weights"], [1, 0, 0])
    assert np.allclose(narrow["area"], 10)


def test_flux_integral():
    fit = {"theta": np.array([math.log(0.02), 0.0, 0.0, math.log(0.7), -10.0, -2.0])}
    result = luminosities(fit, 1e12)
    for row in result.values():
        expected = 0.02 * math.pi * (row["radius_arcsec"] / 60) ** 2
        assert np.isclose(row["model_rate_ct_s"], expected)
        assert np.isclose(row["energy_flux_erg_cm2_s"] / (expected / 1e12), 1)
        factor = 4 * math.pi * (PARAMS["distance_au"] * AU_CM) ** 2
        assert np.isclose(row["luminosity_erg_s"] / (factor * expected / 1e12), 1)
    assert edges()[0] == 0 and edges()[-1] == 12


def test_response_channels():
    channels = np.arange(101)
    narrow = channel_weights(channels, 200, 350)
    assert np.isclose(narrow[40], 0.8) and np.isclose(narrow[70], 0.2)
    assert np.isclose(narrow.sum(), 30)
    assert np.allclose(
        narrow + channel_weights(channels, 350, 500),
        channel_weights(channels, 200, 500),
    )


if __name__ == "__main__":
    test_selection()
    test_native_exposure_and_patterns()
    test_flux_integral()
    test_response_channels()
