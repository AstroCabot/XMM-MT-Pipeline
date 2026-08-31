import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

import lightcurve
import run as analysis
import solar_wind
from figure import (
    LINES,
    MAP_LONGITUDE_POINTS,
    MAP_RADIAL_POINTS,
    PLOT_LIMITS_UTC,
    _hourly,
    _lightcurve_edges,
    nearest_source_epoch,
    parker_field,
)


def raises(error, function, *args):
    try:
        function(*args)
    except error:
        return
    raise AssertionError(f"{function.__name__} did not raise {error.__name__}")


def fixed_track(longitude, radius):
    return {
        "mjd": np.linspace(-10, 20, 301),
        "lon": np.full(301, longitude),
        "r": np.full(301, radius),
    }


def test_wrap_and_mapping():
    assert solar_wind.wrap180(181) == -179
    assert solar_wind.wrap180(-181) == 179
    series = {
        "mjd": np.array([0.0, 1.0]),
        "density": np.array([5.0, 5.0]),
        "speed": np.array([400.0, 400.0]),
        "segment": np.array([0, 0]),
    }
    mapped = solar_wind.map_samples(
        series, fixed_track(0, 1), fixed_track(solar_wind.OMEGA_DEG_DAY, 1)
    )
    assert np.allclose(mapped["arrival"], series["mjd"] + 1)
    assert np.allclose(mapped["flux"], 2e8)
    assert "delay" not in mapped

    speed = solar_wind.AU_KM / 86400
    radial = {**series, "speed": np.array([speed, speed])}
    mapped = solar_wind.map_samples(radial, fixed_track(0, 1), fixed_track(0, 2))
    assert np.allclose(mapped["arrival"], radial["mjd"] + 1)
    assert np.allclose(mapped["flux"], 5 * speed * 1e5 / 4)


def test_segment_smoothing():
    hours = np.r_[np.arange(5), np.arange(8, 13)] / 24
    mapped = {
        "arrival": hours,
        "flux": np.r_[np.ones(5), np.full(5, 3.0)],
        "segment": np.r_[np.zeros(5, int), np.ones(5, int)],
    }
    grid = np.arange(0, 12 * 6 + 1) / 144
    smooth = solar_wind.smooth_segments(mapped, grid)
    assert solar_wind.SMOOTH_MINUTES // solar_wind.GRID_MINUTES + 1 == 13
    assert np.all(np.isnan(smooth[(grid > 4 / 24) & (grid < 8 / 24)]))
    assert np.allclose(smooth[(grid >= 1 / 24) & (grid <= 3 / 24)], 1)
    assert np.allclose(smooth[(grid >= 9 / 24) & (grid <= 11 / 24)], 3)


def test_f107_proxy():
    values = solar_wind.f107()
    daily = values["daily"]
    assert set(daily) == set(solar_wind.F107_DATES)
    for day in daily.values():
        assert np.isclose(
            day["f107p_sfu"], 0.5 * (day["f107_sfu"] + day["f107a_sfu"]), atol=1e-4
        )
    assert np.isclose(
        values["f107p_sfu"],
        np.mean([day["f107p_sfu"] for day in daily.values()]),
        atol=1e-4,
    )
    assert np.isclose(daily["2025-12-03"]["f107p_sfu"], 173.27037)
    assert np.isclose(daily["2025-12-04"]["f107p_sfu"], 178.15062)
    assert np.isclose(values["f107p_sfu"], 175.7105)


def test_declared_series_regression():
    with tempfile.TemporaryDirectory(dir=solar_wind.HERE) as directory:
        result = solar_wind.products(Path(directory))
        assert {path.name for path in Path(directory).iterdir()} == {"mapped_flux.tsv"}
    stats = result["stats"]
    expected = {
        "stereo_a": (8.340256284e7, 3.563975394),
        "solar_orbiter": (1.404654165e8, 4.018791085),
        "ace": (3.260094742e7, 3.071594967),
        "ace_scaled": (5.542161062e7, 3.071594967),
        "wind": (6.019930670e7, 2.585809248),
    }
    for name, (median, dynamic) in expected.items():
        assert np.isclose(stats[name]["median"], median, rtol=1e-9)
        assert np.isclose(stats[name]["dynamic_range"], dynamic, rtol=1e-9)
    assert solar_wind.EVAL_PAD_DAYS == 3
    limits = Time(PLOT_LIMITS_UTC, format="isot", scale="utc").mjd
    for name, values in result["curves"].items():
        for limit in limits:
            index = np.argmin(abs(result["grid"] - limit))
            if name == "wind" and not np.isfinite(values[index]):
                finite = result["grid"][np.isfinite(values)]
                assert np.min(abs(finite - limit)) < 6 / 24
            else:
                assert np.isfinite(values[index])
    quality = result["ace_quality"]
    assert quality["source_samples"] == 22
    assert quality["quality"] == "preliminary browse"
    assert quality["analysis_scale"] == 1.7
    assert quality["source_interval_utc"] == [
        "2025-11-24T02:00:00.000",
        "2025-11-26T09:00:00.000",
    ]
    assert np.isclose(quality["density_ace_over_wind"], 0.577502056276482)
    assert np.isclose(quality["speed_ace_over_wind"], 0.9886269037979841)


def test_parker_support():
    assert {name: style[3] for name, style in LINES.items()} == {
        "stereo_a": 5,
        "solar_orbiter": 4,
        "wind": 3,
        "ace_scaled": 2,
        "ace": 2,
    }
    source = np.linspace(-30, 30, 121)
    mapped = {
        "source": source,
        "lon": np.zeros(source.size),
        "r": np.ones(source.size),
        "speed": np.full(source.size, 400.0),
        "source_flux": np.ones(source.size),
    }
    longitude, radius, field = parker_field(mapped, 0)
    assert longitude.size == MAP_LONGITUDE_POINTS
    assert radius.size == MAP_RADIAL_POINTS
    assert np.all(np.isfinite(field))
    assert PLOT_LIMITS_UTC == ("2025-12-01T23:00:00", "2025-12-06T01:00:00")
    raises(RuntimeError, parker_field, mapped, 31.0)

    hourly = {
        **mapped,
        "source": np.arange(24) / 24,
        "lon": np.zeros(24),
        "r": np.ones(24),
        "speed": np.full(24, 400.0),
        "source_flux": np.ones(24),
    }
    requested = (23 + 8 / 60) / 24
    assert nearest_source_epoch(hourly, requested) == 23 / 24
    raises(RuntimeError, nearest_source_epoch, hourly, 2.0)

    time, values = _hourly(
        np.array([0.0, 1 / 48, 1 / 24]), np.array([1.0, 3.0, np.nan])
    )
    assert np.allclose(time, [1 / 48, 1.5 / 24])
    assert np.allclose(values[0], 2) and np.isnan(values[1])
    assert np.allclose(_lightcurve_edges([1.0, 2.0, 4.0]), [0.5, 1.5, 3.0, 5.0])


def test_ion_fraction():
    ion = analysis.ION_FRACTION
    assert set(ion) == {
        "units",
        "contemporaneous",
        "schwadron_cravens",
        "ace_swics",
        "solar_orbiter_his",
        "adopted_range",
    }
    assert ion["units"] == "dimensionless" and ion["contemporaneous"] is False
    thresholds = ion["schwadron_cravens"]["principal_transition_energy_gt_eV"]
    sc00 = (
        (459, 0.318, 0.085),
        (348, 0.210, 0.440),
        (640, 0.006, 0),
        (518, 0.058, 0.011),
        (74, 0.065, 0.127),
        (871, 0.070, 0),
        (709, 0.200, 0.030),
        (103, 0.730, 0.970),
        (201, 0.084, 0.102),
        (174, 0.004, 0.005),
        (368, 0.098, 0.029),
        (328, 0.052, 0.044),
        (266, 0.041, 0.028),
        (225, 0.017, 0.007),
        (187, 0.009, 0.003),
        (356, 0.021, 0.024),
        (312, 0.049, 0.045),
        (255, 0.057, 0.022),
        (208, 0, 0.002),
        (565, 0, 0.001),
        (505, 0.005, 0.008),
        (447, 0.016, 0.027),
        (379, 0.019, 0.023),
        (328, 0.006, 0.005),
        (281, 0.002, 0.001),
        (310, 0.002, 0.005),
        (285, 0.007, 0.017),
        (237, 0.023, 0.025),
        (197, 0.031, 0.025),
        (176, 0.041, 0.015),
        (97, 0.034, 0.005),
        (80, 0.007, 0.001),
    )
    for limit in (200, 50):
        values = [
            sum(row[column] for row in sc00 if row[0] > limit) / abundance
            for column, abundance in ((1, 1780), (2, 1550))
        ]
        assert thresholds[str(limit)] == {
            name: float(f"{value:.4g}")
            for name, value in zip(("slow", "fast"), values, strict=True)
        }
    expected = {
        "ace_swics": (
            "ACE/SWICS 2.0 AC_H3_SW2 v04 two-hour",
            ["2025-01-01T00:46:42", "2025-03-16T22:40:10"],
            {"all": 888, "retained": 370},
            [1.386e-4, 2.491e-4, 4.395e-4],
        ),
        "solar_orbiter_his": (
            "Solar Orbiter SWA/HIS L3 V04 ten-minute",
            ["2023-04-01T00:00:02", "2023-04-27T14:10:15"],
            {"all": 2802, "retained": 2441},
            [4.694e-4, 6.053e-4, 9.664e-4],
        ),
    }
    for name, (dataset, window, records, values) in expected.items():
        assert ion[name]["dataset"] == dataset
        assert ion[name]["window_utc"] == window and ion[name]["records"] == records
        assert np.array_equal(ion[name]["f_16_50_84"], values)
        assert np.all(np.diff(values) > 0)
    assert ion["ace_swics"]["o_h_scatter_fwhm_dex"] == 0.501
    assert ion["solar_orbiter_his"]["includes_sub_200_eV_charge_states"] is True
    assert ion["adopted_range"] == [3e-4, 1e-3]


def test_lightcurve_equation():
    value, error = lightcurve.rate(100, 50, 10, 4, 2000, 3000)
    area = np.pi * (370 / 60) ** 2
    net = 100 - 0.063 * 10, 50 - 0.063 * 4
    variance = (100 + 0.063**2 * 10) / 2000**2 + (50 + 0.063**2 * 4) / 3000**2
    assert np.isclose(value, area * (net[0] / 2000 - net[1] / 3000))
    assert np.isclose(error, area * np.sqrt(variance))


def test_geometry():
    assert np.isclose(lightcurve.separation(0, 0, 1, 0), 3600)
    assert np.isclose(lightcurve.separation(12, -3, 12, -3), 0)
    assert abs(lightcurve.disk(20).sum() - np.pi * 400) / (np.pi * 400) < 0.02
    ring = lightcurve.annulus(12, 20)
    assert set(np.unique(ring)) <= {0.0, 1.0}
    assert abs(ring.sum() - np.pi * (400 - 144)) / (np.pi * 256) < 0.05


def test_input_contracts():
    with tempfile.TemporaryDirectory(dir=lightcurve.HERE) as directory:
        reduction = Path(directory)
        frame_dir = reduction / "outputs/03_events"
        source_dir = reduction / "outputs/04_sources"
        frame_dir.mkdir(parents=True)
        source_dir.mkdir(parents=True)
        frames = [{"frame": name, "det": "pn"} for name in lightcurve.EXPECTED]
        masked = [
            {"frame": name, "det": "pn", "event": name} for name in lightcurve.EXPECTED
        ]
        solar_wind.write_tsv(frame_dir / "frames.tsv", tuple(frames[0]), frames)
        solar_wind.write_tsv(source_dir / "events.tsv", tuple(masked[0]), masked)
        sources = [
            {
                "ra": index,
                "dec": index / 2,
                "radius_arcsec": lightcurve.PARAMS["source_mask_radius_arcsec"],
            }
            for index in range(15)
        ]
        solar_wind.write_tsv(source_dir / "sources.tsv", tuple(sources[0]), sources)
        old = lightcurve.REDUCTION
        lightcurve.REDUCTION = reduction
        try:
            assert len(lightcurve.pn_frames()) == 8 and len(lightcurve.sources()) == 15
            solar_wind.write_tsv(
                frame_dir / "frames.tsv", tuple(frames[0]), [frames[0], *frames[:-1]]
            )
            raises(RuntimeError, lightcurve.pn_frames)
            solar_wind.write_tsv(
                source_dir / "sources.tsv",
                tuple(sources[0]),
                [{**sources[0], "ra": np.nan}, *sources[1:]],
            )
            raises(RuntimeError, lightcurve.sources)
        finally:
            lightcurve.REDUCTION = old


def test_event_selection_and_exposure():
    with tempfile.TemporaryDirectory(dir=lightcurve.HERE) as directory:
        directory = Path(directory)
        header = fits.Header(
            {
                "EXTNAME": "EVENTS",
                "REFXCRVL": 0.0,
                "REFYCRVL": 0.0,
                "REFXCRPX": 10000.0,
                "REFYCRPX": 10000.0,
                "REFXCDLT": -0.05 / 3600,
                "REFYCDLT": 0.05 / 3600,
                "REFXCTYP": "RA---TAN",
                "REFYCTYP": "DEC--TAN",
                "MJDREF": lightcurve.XMM_MJDREF,
                "TIMESYS": "TT",
            }
        )
        columns = [
            fits.Column(
                name="X", format="J", array=[10000, 20000, 40000, 10000, 10000]
            ),
            fits.Column(name="Y", format="J", array=np.full(5, 10000)),
            fits.Column(name="TIME", format="D", array=np.zeros(5)),
            fits.Column(name="PI", format="J", array=[500, 500, 500, 1100, 500]),
            fits.Column(name="PATTERN", format="B", array=[0, 0, 0, 0, 1]),
            fits.Column(name="FLAG", format="J", array=np.zeros(5, int)),
        ]
        event = directory / "event.fits"
        fits.HDUList(
            [
                fits.PrimaryHDU(),
                fits.BinTableHDU.from_columns(columns, header=header, name="EVENTS"),
            ]
        ).writeto(event)
        epoch = lightcurve.event_mjd(np.array([0.0]), header)[0]
        eph = {
            "mjd": np.array([epoch - 1, epoch + 1]),
            "ra_deg": np.zeros(2),
            "dec_deg": np.zeros(2),
        }
        # The PATTERN==1 event is rejected: PN below 1.1 keV keeps singles only.
        assert lightcurve.counts(event, eph) == (1, 1)
        assert lightcurve.valid_fits(event)
        assert not lightcurve.valid_fits(directory / "missing.fits")

        wcs = WCS(naxis=2)
        wcs.wcs.crval = [0, 0]
        wcs.wcs.crpix = [401, 401]
        wcs.wcs.cdelt = [-4 / 3600, 4 / 3600]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        exposure = directory / "exposure.fits"
        fits.PrimaryHDU(np.full((801, 801), 10.0, dtype="f4"), wcs.to_header()).writeto(
            exposure
        )
        source, background = lightcurve.region_exposure(
            exposure, np.array([epoch]), np.array([1.0]), eph
        )
        source_area = np.pi * (370 / 60) ** 2
        background_area = np.pi * (720**2 - 480**2) / 3600
        assert abs(source / 10 - source_area) / source_area < 0.01
        assert abs(background / 10 - background_area) / background_area < 0.01

        empty = directory / "empty.fits"
        fits.PrimaryHDU(np.zeros((801, 801), dtype="f4"), wcs.to_header()).writeto(
            empty
        )
        raises(
            RuntimeError,
            lightcurve.region_exposure,
            empty,
            np.array([epoch]),
            np.array([1.0]),
            eph,
        )


def test_gti_exposure_cache_and_contract():
    with tempfile.TemporaryDirectory(dir=lightcurve.HERE) as directory:
        directory = Path(directory)
        event = directory / "event.fits"
        event_header = fits.Header(
            {"DSTYP2": "TIME", "DSREF2": ":GTI00003", "ONTIME": 10.0}
        )
        columns = [
            fits.Column(name=key, format="D", array=values)
            for key, values in (("START", [0.0, 20.0]), ("STOP", [5.0, 25.0]))
        ]
        broad = [
            fits.Column(name=key, format="D", array=values)
            for key, values in (("START", [0.0]), ("STOP", [30.0]))
        ]
        fits.HDUList(
            [
                fits.PrimaryHDU(),
                fits.BinTableHDU.from_columns([], header=event_header, name="EVENTS"),
                fits.BinTableHDU.from_columns(broad, name="STDGTI01"),
                fits.BinTableHDU.from_columns(columns, name="GTI00003"),
            ]
        ).writeto(event)
        times, weights = lightcurve.gti_samples(
            {"event": event, "ontime": 10, "frame": "unit"}
        )
        assert np.allclose(times, [2.5, 22.5]) and np.allclose(weights, [5, 5])

        reduction = directory / "reduction"
        source_dir = reduction / "outputs/04_sources"
        source_dir.mkdir(parents=True)
        solar_wind.write_tsv(
            source_dir / "sources.tsv",
            ("ra", "dec", "radius_arcsec"),
            [
                {
                    "ra": 0 if index == 0 else index * 10,
                    "dec": 0 if index == 0 else index,
                    "radius_arcsec": lightcurve.PARAMS["source_mask_radius_arcsec"],
                }
                for index in range(15)
            ],
        )
        work = directory / "work"
        calls = []
        wcs = WCS(naxis=2)
        wcs.wcs.crval = [0, 0]
        wcs.wcs.crpix = [11, 11]
        wcs.wcs.cdelt = [-4 / 3600, 4 / 3600]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        header = wcs.to_header()
        header["RADECSYS"] = "FK5"
        header.pop("RADESYS", None)

        def fake_sas(tool, cwd, env, **parameters):
            calls.append((tool, parameters))
            if tool == "evselect":
                fits.PrimaryHDU(np.ones((21, 21), dtype="i4"), header).writeto(
                    parameters["imageset"]
                )
            else:
                fits.PrimaryHDU(np.full((21, 21), 100.0, dtype="f4"), header).writeto(
                    parameters["expimageset"]
                )

        original = lightcurve.WORK, lightcurve.REDUCTION, lightcurve.sas
        lightcurve.WORK, lightcurve.REDUCTION, lightcurve.sas = (
            work,
            reduction,
            fake_sas,
        )
        row = {"frame": "unit", "event": event}
        try:
            path = lightcurve.exposure_map(row, {})
            data = fits.getdata(path)
            assert data[10, 10] == 0 and data[0, 0] > 0
            saved = fits.getheader(path)
            assert (
                not saved["VIGNET"]
                and saved["SRCMASK"]
                and (saved["PI_MIN"], saved["PI_MAX"], saved["PATTERN"])
                == (200, 1000, 0)
            )
            assert "RADECSYS" not in saved and saved["RADESYS"] == "FK5"
            assert (
                calls[0][0] == "evselect" and "PATTERN==0" in calls[0][1]["expression"]
            )
            assert calls[1][0] == "eexpmap" and calls[1][1]["withvignetting"] is False
            with fits.open(path, mode="update") as hdul:
                hdul[0].header["VIGNET"] = True
            lightcurve.exposure_map(row, {})
            assert len(calls) == 4 and not fits.getheader(path)["VIGNET"]
        finally:
            lightcurve.WORK, lightcurve.REDUCTION, lightcurve.sas = original

        wind = {
            "stats": {"wind": {"median": 6e7, "dynamic_range": 2.5}},
            "ace_quality": {"source_samples": 2},
            "f107": {"f107p_sfu": 175.7105},
        }
        result = analysis.summarize(
            wind,
            [{"rate_ct_s": "1.0"}, {"rate_ct_s": "2.0"}],
            lightcurve.BAND_EV,
            lightcurve.PATTERN,
        )
        assert result["lightcurve"]["dynamic_range"] == 2
        assert result["heavy_ion_prior"] == analysis.SHARED["heavy_ion_prior"]
        assert result["ion_fraction"] == {
            "adopted_range": analysis.ION_FRACTION["adopted_range"]
        }
        assert result["f107"]["f107p_sfu"] == 175.7105
        raises(
            RuntimeError,
            analysis.summarize,
            wind,
            [{"rate_ct_s": "-1"}],
            lightcurve.BAND_EV,
            lightcurve.PATTERN,
        )


def test_disposable_frame_fits():
    with tempfile.TemporaryDirectory(dir=lightcurve.HERE) as directory:
        directory = Path(directory)
        work, output = directory / "work", directory / "outputs"
        event = directory / "event.fits"
        header = fits.Header({"MJDREF": lightcurve.XMM_MJDREF, "TIMESYS": "TT"})
        fits.HDUList(
            [
                fits.PrimaryHDU(),
                fits.BinTableHDU.from_columns([], header=header, name="EVENTS"),
            ]
        ).writeto(event)
        frame = {"frame": "pn-unit", "event": event, "oot": event}

        def fake_exposure(row, env):
            paths = lightcurve.frame_fits(row["frame"])
            for path in paths:
                path.parent.mkdir(parents=True, exist_ok=True)
                fits.PrimaryHDU(np.ones((2, 2), dtype="f4")).writeto(path)
            return paths[0]

        original = (
            lightcurve.WORK,
            lightcurve.OUTPUT,
            lightcurve.ephemeris,
            lightcurve.sas_environment,
            lightcurve.gti_samples,
            lightcurve.counts,
            lightcurve.exposure_map,
            lightcurve.region_exposure,
        )
        lightcurve.WORK, lightcurve.OUTPUT = work, output
        lightcurve.ephemeris = lambda: {}
        lightcurve.sas_environment = lambda: {}
        lightcurve.gti_samples = lambda row: (np.array([0.0]), np.array([1.0]))
        lightcurve.counts = lambda path, eph: (100, 50)
        lightcurve.exposure_map = fake_exposure
        lightcurve.region_exposure = lambda *args: (2000, 3000)
        try:
            rows = lightcurve.build([frame])
            assert len(rows) == 1 and (output / "lightcurve.tsv").is_file()
            assert not any(path.exists() for path in lightcurve.frame_fits("pn-unit"))

            def fail(*args):
                raise ValueError("synthetic failure")

            lightcurve.region_exposure = fail
            raises(ValueError, lightcurve.build, [frame])
            assert not any(path.exists() for path in lightcurve.frame_fits("pn-unit"))

            outside = directory / "outside.fits"
            outside.write_bytes(b"unchanged")
            raises(RuntimeError, lightcurve.remove_frame_fits, "../../outside")
            assert outside.read_bytes() == b"unchanged"
        finally:
            (
                lightcurve.WORK,
                lightcurve.OUTPUT,
                lightcurve.ephemeris,
                lightcurve.sas_environment,
                lightcurve.gti_samples,
                lightcurve.counts,
                lightcurve.exposure_map,
                lightcurve.region_exposure,
            ) = original


if __name__ == "__main__":
    assert lightcurve.SETTINGS["python_lib"] == "/home/scabot/miniconda3/envs/xmm/lib"
    test_wrap_and_mapping()
    test_segment_smoothing()
    test_f107_proxy()
    test_declared_series_regression()
    test_parker_support()
    test_ion_fraction()
    test_lightcurve_equation()
    test_geometry()
    test_input_contracts()
    test_event_selection_and_exposure()
    test_gti_exposure_cache_and_contract()
    test_disposable_frame_fits()
