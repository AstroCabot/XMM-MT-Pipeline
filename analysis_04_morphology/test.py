import json
import shutil

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

import figure
import run as morphology
import smooth

WORK = (morphology.REPRO / "work" / "morphology_tests").resolve()


def raises(error, function, *args):
    try:
        function(*args)
    except error:
        return
    raise AssertionError(f"{error.__name__} was not raised")


def frame(shape=(101, 101)):
    wcs = WCS(naxis=2)
    wcs.wcs.crval = [177.8, 1.4]
    wcs.wcs.crpix = [(shape[1] + 1) / 2, (shape[0] + 1) / 2]
    wcs.wcs.cdelt = [-4 / 3600, 4 / 3600]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    header = wcs.to_header()
    header["NAXIS1"], header["NAXIS2"] = shape[1], shape[0]
    header["CRPIX1"], header["CRPIX2"] = wcs.wcs.crpix
    return header


def test_combination():
    one = np.ones((3, 3))
    products = {}
    for detector in morphology.DETECTORS:
        products[detector] = {
            "events": one,
            "qpb": 0.1 * one,
            "sp": 0.2 * one,
            "exposure": 10 * one,
            "oot": np.zeros_like(one),
        }
    products["pn"]["oot"] = 2 * one
    conversion = {"pn": 2, "mos1": 3, "mos2": 4}
    events, qpb, sp, oot, exposure, net, variance = morphology.combine(
        products, conversion
    )
    scale = morphology.PARAMS["pn_oot_fraction"]
    assert np.all(events == 3)
    assert np.allclose(qpb, 0.3) and np.allclose(sp, 0.6)
    assert np.all(oot == 2) and np.all(exposure == 22.5)
    assert np.allclose(net, 3 - 0.3 - 0.6 - 2 * scale)
    assert np.allclose(variance, 3 + 0.3 + 0.6 + 2 * scale**2)


def test_stack_contract():
    root = WORK / "stack"
    root.mkdir()
    paths = []
    low, high = morphology.PARAMS["soft_band_ev"]
    for detector in morphology.DETECTORS:
        for kind in ("events", "qpb", "sp", "exposure"):
            header = frame((3, 3))
            header.update(
                DETECTOR=detector, BAND="soft", E_MIN=low, E_MAX=high, COMOVING=True
            )
            path = root / f"{detector}-{kind}-soft.fits"
            fits.PrimaryHDU(np.ones((3, 3)), header).writeto(
                path, output_verify="silentfix"
            )
            paths.append(path)
    header = frame((3, 3))
    header.update(DETECTOR="pn", BAND="soft", E_MIN=low, E_MAX=high, COMOVING=True)
    path = root / "pn-oot-soft.fits"
    fits.PrimaryHDU(np.ones((3, 3)), header).writeto(path, output_verify="silentfix")
    paths.append(path)

    original, morphology.STACK = morphology.STACK, root
    try:
        products, saved = morphology.load_stack()
        assert set(products) == set(morphology.DETECTORS)
        assert saved["DETECTOR"] == "EPIC"
        fits.setval(paths[0], "E_MIN", value=low + 1)
        raises(RuntimeError, morphology.load_stack)
        fits.setval(paths[0], "E_MIN", value=low)
        fits.setval(paths[0], "CRPIX1", value=2.5)
        raises(RuntimeError, morphology.load_stack)
        fits.setval(paths[0], "CRPIX1", value=2.0)
        for item in paths:
            fits.setval(item, "CDELT1", value=-8 / 3600)
            fits.setval(item, "CDELT2", value=8 / 3600)
        raises(RuntimeError, morphology.load_stack)
    finally:
        morphology.STACK = original


def test_mask_and_smoothing():
    shape = (61, 61)
    exposure = np.full(shape, 100.0)
    exposure[:25] = 5
    pn = np.full(shape, 100.0)
    pn[:, 30] = 0
    mask = morphology.display_mask(exposure, pn)
    assert not mask[:25].any() and not mask[:, 30].any() and mask[40, 10]
    assert morphology.display_mask(pn, exposure)[5, 10]


def test_adaptive_display_plane():
    root = WORK / "adaptive"
    background = root / "background/mos1-x/soft"
    background.mkdir(parents=True)
    template = background / f"mos1S001-fovimsky-{smooth.SOFT_TAG}.fits"
    fits.PrimaryHDU(np.zeros((2, 2)), frame((8, 8))).writeto(
        template, output_verify="silentfix"
    )
    original = smooth.REPRO, smooth.BACKGROUND, smooth.run

    def fake(command, work):
        if "withmask=yes" in command:
            path = next(
                value.split("=", 1)[1]
                for value in command
                if value.startswith("maskfile=")
            )
            saved = fits.getdata(path)
            assert saved.sum() == 32 and fits.getheader(path)["BITPIX"] == 32
        else:
            assert "withmask=no" in command
            assert not any(value.startswith("maskfile=") for value in command)
        assert (
            fits.getheader(work / f"comb-fovimsky-{smooth.SOFT_TAG}.fits")["BITPIX"]
            == -32
        )
        fits.PrimaryHDU(np.ones((smooth.SIZE // 2,) * 2)).writeto(
            work / f"comb-adaptimsky-{smooth.SOFT_TAG}.fits"
        )

    try:
        smooth.REPRO, smooth.BACKGROUND, smooth.run = root, root / "background", fake
        one = np.ones((8, 8))
        mask = np.indices((8, 8))[1] < 4
        result = smooth.adaptive(
            one, 0.1 * one, 0.2 * one, np.zeros_like(one), 10 * one, frame((8, 8))
        )
        assert np.all(np.isfinite(result))
        result = smooth.adaptive(
            one, 0.1 * one, 0.2 * one, np.zeros_like(one), 10 * one, frame((8, 8)), mask
        )
        assert np.all(np.isfinite(result[mask])) and np.all(np.isnan(result[~mask]))
    finally:
        smooth.REPRO, smooth.BACKGROUND, smooth.run = original


def test_surface_brightness():
    shape = (81, 101)  # non-square: catches a transposed crop
    counts, seconds, pixel_arcsec = 200, 1000.0, 4.0
    per_deg2 = counts / seconds / (pixel_arcsec / 3600) ** 2
    x0, y0 = (smooth.SIZE - shape[1]) // 2, (smooth.SIZE - shape[0]) // 2
    blocked = np.zeros((smooth.SIZE // smooth.BINFACTOR,) * 2)
    blocked[
        y0 // smooth.BINFACTOR : (y0 + shape[0]) // smooth.BINFACTOR,
        x0 // smooth.BINFACTOR : (x0 + shape[1]) // smooth.BINFACTOR,
    ] = per_deg2
    field = smooth.surface_brightness(blocked, shape, x0, y0)
    assert field.shape == shape
    interior = field[2:-2, 2:-2]
    assert np.allclose(interior, counts / seconds / (pixel_arcsec / 60) ** 2)
    assert not np.allclose(interior, per_deg2)  # the 3600 must be applied


def test_slice_and_products():
    header = frame()
    dx, dy = morphology.coordinates((101, 101), header)
    radius = np.hypot(dx, dy)
    exposure = np.full(radius.shape, 100.0)
    brightness = 1 / (radius + 0.5)
    pixel_area = (4 / 60) ** 2
    net = brightness * exposure * pixel_area
    variance = np.ones(radius.shape)
    mask = np.ones(radius.shape, bool)
    profile = morphology.slice_rows(
        net, variance, exposure, brightness, mask, dx, dy, pixel_area
    )
    table = np.asarray(profile)
    assert np.allclose(table[:, 2], table[::-1, 2], rtol=0.02, equal_nan=True)
    assert np.allclose(table[:, 5], table[:, 6], rtol=0.02, equal_nan=True)

    angle = np.deg2rad(morphology.SUN_PA_DEG)
    along = dx * np.sin(angle) + dy * np.cos(angle)
    across = dx * np.cos(angle) - dy * np.sin(angle)
    asymmetric = np.exp(-((along - 2) ** 2 + across**2) / 0.08)
    directed = np.asarray(
        morphology.slice_rows(
            asymmetric * exposure * pixel_area,
            variance,
            exposure,
            asymmetric,
            mask,
            dx,
            dy,
            pixel_area,
        )
    )
    assert directed[np.nanargmax(directed[:, 5]), 0] > 1.5

    morphology.OUTPUT = WORK
    arrays = tuple(
        (name, np.ones(radius.shape), "count")
        for name in (
            "EVENTS",
            "QPB",
            "SP",
            "OOT",
            "EXPOSURE",
            "DISPLAY_MASK",
            "NET",
            "VARIANCE",
            "SMOOTH",
        )
    )
    morphology.write_products(
        arrays, header, profile, {"pn": 1.0, "mos1": 2.0, "mos2": 4.0}, 6.0
    )
    with fits.open(WORK / "image.fits") as hdul:
        assert [h.name for h in hdul[1:]] == [a[0] for a in arrays]
    result = json.loads((WORK / "result.json").read_text())
    assert result["detector_exposure_weights"] == {"pn": 0.25, "mos1": 0.5, "mos2": 1.0}
    assert result["contour_seed_radius_arcmin"] == 0.75
    figure.draw(
        brightness,
        mask,
        dx,
        dy,
        profile,
        morphology.SUN_PA_DEG,
        morphology.MOTION_PA_DEG,
        morphology.PARAMS["science_aperture_radius_arcsec"] / 60,
        WORK,
    )
    assert (WORK / "figure1.pdf").stat().st_size > 0
    assert (WORK / "figure1.png").stat().st_size > 0


def test_display_profile_keeps_excluded_bins():
    header = frame((401, 401))
    dx, dy = morphology.coordinates((401, 401), header)
    angle = np.deg2rad(morphology.SUN_PA_DEG)
    along = dx * np.sin(angle) + dy * np.cos(angle)
    exposure = np.full(along.shape, 100.0)
    brightness = 1 / (np.hypot(dx, dy) + 0.5)
    mask = ~((along >= 1.0) & (along < 1.25))
    masked = morphology.slice_rows(
        brightness * exposure,
        np.ones_like(brightness),
        exposure,
        np.where(mask, brightness, np.nan),
        mask,
        dx,
        dy,
        1.0,
        max_radius=3,
    )
    display = morphology.slice_rows(
        brightness * exposure,
        np.ones_like(brightness),
        exposure,
        brightness,
        np.ones_like(mask),
        dx,
        dy,
        1.0,
        max_radius=3,
    )
    masked, display = np.asarray(masked), np.asarray(display)
    gap = np.argmin(abs(display[:, 0] - 1.125))
    assert np.isnan(masked[gap, 2])
    assert np.isfinite(display[gap, 2])
    assert np.isfinite(display[gap, 6])


def test_reflected_curve_uses_raw_profile():
    import matplotlib.pyplot as plt

    distance = np.array([-0.75, -0.25, 0.25, 0.75])
    raw = np.array([1.0, 2.0, 3.0, 4.0]) * 1e-3
    smooth = 10 * raw
    coverage = np.array([1.0, 0.4, 0.5, 1.0])
    profile = np.column_stack(
        [distance, distance, raw, 0.1 * raw, coverage, smooth, smooth[::-1]]
    )
    fig, ax = plt.subplots()
    figure._profile_panel(ax, profile, 0, 1, 2, 1)
    reflected = next(line for line in ax.lines if line.get_label() == "Reflected")
    assert np.allclose(reflected.get_xdata(), distance)
    assert np.allclose(
        reflected.get_ydata(),
        np.where(coverage >= figure.MIN_PROFILE_COVERAGE, raw * 1e3, np.nan),
        equal_nan=True,
    )
    plt.close(fig)


def test_source_extent():
    yy, xx = np.indices((121, 121))
    dx, dy = (xx - 60) / 10, (yy - 60) / 10
    radius = np.hypot(dx, dy)
    detected = radius < 3
    detected[:, 60] = False
    detected[5:8, 5:8] = True
    extent = morphology.source_extent(detected, dx, dy)
    assert 2.8 < extent < 3.1


def test_rotated_mask():
    header = frame()
    dx, dy = morphology.coordinates((101, 101), header)
    mask = dx < 0
    values, support, _, _ = figure.rotated(
        np.where(mask, 1.0, np.nan), mask, dx, dy, morphology.SUN_PA_DEG, 3
    )
    assert np.allclose(values[support], 1)
    assert np.all(np.isnan(values[~support]))


def main():
    assert smooth.SETTINGS["python_lib"] == "/home/scabot/miniconda3/envs/xmm/lib"
    assert smooth.SOFT_BAND == tuple(morphology.PARAMS["soft_band_ev"])
    assert smooth.SOFT_TAG == "-".join(map(str, smooth.SOFT_BAND))
    assert WORK.parent == (morphology.REPRO / "work").resolve()
    if WORK.exists():
        shutil.rmtree(WORK)
    try:
        WORK.mkdir(parents=True)
        test_combination()
        test_stack_contract()
        test_mask_and_smoothing()
        test_adaptive_display_plane()
        test_surface_brightness()
        test_slice_and_products()
        test_display_profile_keeps_excluded_bins()
        test_reflected_curve_uses_raw_profile()
        test_source_extent()
        test_rotated_mask()
    finally:
        if WORK.exists():
            shutil.rmtree(WORK)


if __name__ == "__main__":
    main()
