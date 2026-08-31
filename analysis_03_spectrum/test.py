import contextlib
import json
import math
import os
import shutil
import sys
import tempfile
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS

import extract
import fit
import inputs
import model
import run as pipeline
import spectrum

HERE = Path(__file__).resolve().parent

# Response mixing: 3% in rate and 5% in normalized shape.
RATE_TOLERANCE = 0.03
SHAPE_TOLERANCE = 0.05

rows = [{"frame": "a"}, {"frame": "b"}]
assert set(inputs.exact_index(rows, lambda row: row["frame"], {"a", "b"}, "frame")) == {
    "a",
    "b",
}
for bad in (rows[:1], [rows[0], rows[0]]):
    try:
        inputs.exact_index(bad, lambda row: row["frame"], {"a", "b"}, "frame")
    except RuntimeError:
        pass
    else:
        raise AssertionError("incomplete Task 01 table accepted")


def write_pha(path, counts, exposure=100, backscal=10, rate=False, error=None, start=0):
    counts = np.asarray(counts, float)
    values = counts / exposure if rate else counts
    columns = [
        fits.Column(
            name="CHANNEL", format="I", array=np.arange(start, start + len(counts))
        ),
        fits.Column(name="RATE" if rate else "COUNTS", format="E", array=values),
    ]
    if error is not None:
        values = (
            np.asarray(error, float) / exposure if rate else np.asarray(error, float)
        )
        columns.append(fits.Column(name="STAT_ERR", format="E", array=values))
    hdu = fits.BinTableHDU.from_columns(columns, name="SPECTRUM")
    hdu.header.update(
        EXPOSURE=exposure,
        BACKSCAL=backscal,
        SPECDELT=5.0,
        TELESCOP="XMM",
        INSTRUME="EMOS2",
    )
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)


def write_response(path, channels):
    channels = np.asarray(channels, int)
    columns = [
        fits.Column(name="CHANNEL", format="I", array=channels),
        fits.Column(name="E_MIN", format="E", unit="keV", array=channels * 0.005),
        fits.Column(name="E_MAX", format="E", unit="keV", array=(channels + 1) * 0.005),
    ]
    fits.HDUList(
        [fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns, name="EBOUNDS")]
    ).writeto(path)


def write_event(path):
    columns = [
        fits.Column(name="X", format="E", array=[0.0]),
        fits.Column(name="Y", format="E", array=[0.0]),
    ]
    hdu = fits.BinTableHDU.from_columns(columns, name="EVENTS")
    hdu.header.update(
        REFXCRVL=10.0,
        REFYCRVL=20.0,
        REFXCRPX=1.0,
        REFYCRPX=1.0,
        REFXCDLT=-1 / 3600,
        REFYCDLT=1 / 3600,
        REFXCTYP="RA---TAN",
        REFYCTYP="DEC--TAN",
        MJDREF=50814.0,
        TIMESYS="TT",
        ONTIME=10.0,
    )

    def gti(spans, name):
        spans = np.asarray(spans, float)
        columns = [
            fits.Column(name="START", format="D", array=spans[:, 0]),
            fits.Column(name="STOP", format="D", array=spans[:, 1]),
        ]
        return fits.BinTableHDU.from_columns(columns, name=name)

    fits.HDUList(
        [
            fits.PrimaryHDU(),
            hdu,
            gti([[0, 30]], "STDGTI"),
            gti([[5, 10], [20, 25]], "GTI00003"),
        ]
    ).writeto(path)


class Parameter:
    def __init__(self):
        self._values = [1.0, 0.01, 0.0, 0.0, 100.0, 100.0]
        self.frozen, self.link, self.index = False, "", 1

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, value):
        values = list(value) if isinstance(value, (list, tuple)) else [value]
        self._values[: len(values)] = values
        if len(values) > 1:
            self.frozen = values[1] < 0


def component(names):
    item = SimpleNamespace(parameterNames=names)
    for name in names:
        setattr(item, name, Parameter())
    return item


class Model:
    def __init__(self):
        definitions = [
            ("constant", ["factor"]),
            ("constant_2", ["factor"]),
            ("apec", ["kT", "Abundanc", "Redshift", "norm"]),
            ("TBabs", ["nH"]),
            ("apec_5", ["kT", "Abundanc", "Redshift", "norm"]),
            ("powerlaw", ["PhoIndex", "norm"]),
            (
                "vacx2",
                [
                    "temperature",
                    "collnpar",
                    "collntype",
                    "acxmodel",
                    "recombtype",
                    "Hefrac",
                    "H",
                    "He",
                    "C",
                    "N",
                    "O",
                    "Ne",
                    "Mg",
                    "Al",
                    "Si",
                    "S",
                    "Ar",
                    "Ca",
                    "Fe",
                    "Ni",
                    "vbroad",
                    "tbroad",
                    "Redshift",
                    "norm",
                ],
            ),
        ]
        self.componentNames = [name for name, _ in definitions]
        parameters = []
        for name, names in definitions:
            item = component(names)
            setattr(self, name, item)
            parameters.extend(getattr(item, par) for par in names)
        self.parameters = parameters
        self.nParameters = len(parameters)
        for index, parameter in enumerate(parameters, 1):
            parameter.index = index

    def __call__(self, index):
        return self.parameters[index - 1]


def test_geometry():
    assert extract.circle_overlap(10, 2, 20) == 0
    assert math.isclose(extract.circle_overlap(10, 2, 0), math.pi * 4)
    half = extract.circle_overlap(10, 10, 10)
    expected = 2 * 100 * math.acos(0.5) - 50 * math.sqrt(3)
    assert math.isclose(half, expected)
    position = {"ra_deg": 10.0, "dec_deg": 0.0}
    source = {"ra": 10.0, "dec": 0.0, "radius": 2.0}
    assert math.isclose(
        extract.solid_angle((0, 10), position, [source]), math.pi * (100 - 4) / 3600
    )
    assert math.isclose(
        extract.solid_angle((5, 10), position, [source]), math.pi * 75 / 3600
    )


def test_sas_environment(work):
    captured = {}

    def invoke(*args, **kwargs):
        captured.update(kwargs["env"])
        return SimpleNamespace(returncode=0)

    old_run, old_work = extract.subprocess.run, extract.WORK
    extract.subprocess.run, extract.WORK = invoke, work / "sas"
    try:
        extract.run_tool("fake", work)
    finally:
        extract.subprocess.run, extract.WORK = old_run, old_work
    assert all(captured[name] == value for name, value in extract.SAS_ENV.items())


def test_center_and_selection(work):
    event = work / "event.fits"
    write_event(event)
    endpoints = np.array([0.0, 30.0])
    mjd = np.asarray(
        (
            Time(50814, format="mjd", scale="tt") + TimeDelta(endpoints, format="sec")
        ).utc.mjd
    )
    eph = {
        "mjd": mjd,
        "ra_deg": np.array([10.0, 40.0]),
        "dec_deg": np.array([0.0, 3.0]),
        "r_au": np.array([2.0, 4.0]),
        "delta_au": np.array([1.0, 2.0]),
    }
    frame = {"frame": "x", "event": event}
    center = extract.center(frame, eph)
    assert math.isclose(center["duration_s"], 10)
    assert math.isclose(center["ra_deg"], 25, abs_tol=1e-8)
    assert math.isclose(center["dec_deg"], 1.5, abs_tol=1e-8)
    position = {"ra_deg": 10.0, "dec_deg": 20.0}
    sources = [
        {"ra": 10.0, "dec": 20.0, "radius": 45.0},
        {"ra": 30.0, "dec": 20.0, "radius": 45.0},
    ]
    pn = extract.selection(event, "pn", (0, 370), position, sources, "T T T T")
    mos = extract.selection(
        event, "mos2", (400, 600), position, sources, "T F T F T T T"
    )
    assert "PATTERN==0" in pn and "PATTERN<=12" in mos
    assert pn.count("!((X,Y) IN circle") == 1
    assert "annulus" in mos and "FLAG==0" in mos
    assert "CCDNR!=2" in mos and "CCDNR!=4" in mos and "CCDNR" not in pn
    try:
        extract.selection(event, "pn", (0, 370), position, sources, "T F T T")
    except RuntimeError as error:
        assert "invalid chip selection" in str(error)
    else:
        raise AssertionError("disabled PN quadrant accepted")


def test_map_fraction(work):
    wcs = WCS(naxis=2)
    wcs.wcs.crval, wcs.wcs.crpix = [10.0, 20.0], [11.0, 11.0]
    wcs.wcs.cdelt, wcs.wcs.ctype = [-1 / 3600, 1 / 3600], ["RA---TAN", "DEC--TAN"]
    path = work / "map.fits"
    header = wcs.to_header()
    header["RADECSYS"] = header.pop("RADESYS")
    fits.PrimaryHDU(np.ones((21, 21)), header).writeto(path)
    fraction = extract.map_fraction(path, (0, 4.5), {"ra_deg": 10.0, "dec_deg": 20.0})
    y, x = np.indices((21, 21))
    expected = np.count_nonzero((x - 10) ** 2 + (y - 10) ** 2 < 4.5**2) / 441
    assert math.isclose(fraction, expected)
    data, use, pixel_area = extract.map_region(
        path, (0, 4.5), {"ra_deg": 10.0, "dec_deg": 20.0}
    )
    omega = use.sum() * pixel_area
    assert np.all(data[use] == 1)
    assert math.isclose(
        extract.mean_exposure(path, (0, 4.5), {"ra_deg": 10.0, "dec_deg": 20.0}, omega),
        1,
    )


def test_response_arguments(work):
    event, oot = work / "event-args.fits", work / "oot-args.fits"
    write_event(event)
    write_event(oot)
    qpb, sp = work / "qpb-args.pi", work / "sp-args.pi"
    write_pha(qpb, np.ones(16), rate=True, error=np.ones(16))
    write_pha(sp, np.ones(16), error=np.zeros(16))
    wcs = WCS(naxis=2)
    wcs.wcs.crval, wcs.wcs.crpix = [10.0, 20.0], [2.0, 2.0]
    wcs.wcs.cdelt, wcs.wcs.ctype = [-1 / 3600, 1 / 3600], ["RA---TAN", "DEC--TAN"]
    qpb_map, sp_map = work / "qpb-args.fits", work / "sp-args.fits"
    exposure_map = work / "exposure-args.fits"
    fits.PrimaryHDU(np.ones((3, 3)), wcs.to_header()).writeto(qpb_map)
    fits.PrimaryHDU(np.ones((3, 3)), wcs.to_header()).writeto(sp_map)
    fits.PrimaryHDU(np.full((3, 3), 100.0), wcs.to_header()).writeto(exposure_map)
    calls = []

    def tool(name, cwd, **parameters):
        calls.append((name, parameters))
        path = Path(cwd)
        if name == "evselect":
            write_pha(path / parameters["spectrumset"], np.ones(16))
        elif name == "arfgen":
            (path / parameters["arfset"]).write_bytes(b"arf")
        elif name == "rmfgen":
            (path / parameters["rmfset"]).write_bytes(b"rmf")

    old_tool, old_link = extract.run_tool, extract.link
    extract.run_tool = tool
    extract.link = lambda path, target: shutil.copy(target, path)
    try:
        frame = {
            "frame": "pn-test",
            "det": "pn",
            "chips": "T T T T",
            "event": event,
            "oot": oot,
            "position": {"ra_deg": 10.0, "dec_deg": 20.0},
            "qpb_pha": qpb,
            "sp_pha": sp,
            "qpb_map": qpb_map,
            "sp_map": sp_map,
            "exposure_map": exposure_map,
        }
        region = (0, math.sqrt(1 / math.pi))
        row = extract.one(
            frame, "inner", region, [], 15, 0.063, (0.025, 0.05), work / "response"
        )
    finally:
        extract.run_tool, extract.link = old_tool, old_link
    arf = next(parameters for name, parameters in calls if name == "arfgen")
    assert arf["detxbins"] == 400 and arf["detybins"] == 400
    assert arf["detmaptype"] == "flat" and arf["filterdss"] == "yes"
    rmf = next(parameters for name, parameters in calls if name == "rmfgen")
    assert rmf["detmaptype"] == "flat" and rmf["format"] == "var"
    assert math.isclose(row["response_throughput"], 1, rel_tol=1e-9)
    assert row["source_counts"] == 6
    with fits.open(exposure_map, mode="update") as hdul:
        hdul[0].data[1, 1] = 50
    half = extract.mean_exposure(exposure_map, region, frame["position"], 1 / 3600)
    assert math.isclose(half / 100, 0.5, rel_tol=1e-9)


def test_background_and_grouping(work):
    unused = work / "unused-error-negative.pi"
    write_pha(unused, [0, 1], error=[-123, 1])
    assert np.array_equal(spectrum.read(unused).variance, [0, 1])
    unused = work / "unused-error-positive.pi"
    write_pha(unused, [0, 1], error=[1e30, 1])
    assert np.array_equal(spectrum.read(unused).variance, [0, 1])
    invalid = work / "invalid-error.pi"
    write_pha(invalid, [0, 1], error=[0, -1])
    try:
        spectrum.read(invalid)
    except ValueError:
        pass
    else:
        raise AssertionError("negative in-band PHA error accepted")
    source, qpb, sp, oot = (
        work / name for name in ("source.pi", "qpb.pi", "sp.pi", "oot.pi")
    )
    write_pha(source, [100, 200, 300, 400], exposure=100, backscal=12)
    write_pha(qpb, [10, 20, 30, 40], exposure=50, rate=True, error=[1, 2, 3, 4])
    write_pha(sp, [2, 4, 6, 8], exposure=100, error=[0, 0, 0, 0])
    write_pha(oot, [5, 10, 15, 20], exposure=100)
    output = work / "background.pi"
    spectrum.background(output, source, qpb, sp, (0.2, 0.3), oot, 0.063)
    actual = spectrum.read(output, poisson=False)
    expected = (
        0.4 * np.array([10, 20, 30, 40])
        + 0.3 * np.array([2, 4, 6, 8])
        + 0.063 * np.array([5, 10, 15, 20])
    )
    variance = 0.4**2 * np.array([1, 2, 3, 4]) ** 2 + 0.063**2 * np.array(
        [5, 10, 15, 20]
    )
    assert np.allclose(actual.counts, expected)
    assert np.allclose(actual.variance, variance, rtol=1e-6)
    assert actual.exposure == 100 and actual.backscal == 12
    assert spectrum.band_counts(source, 0, 5) == 300
    with fits.open(output) as hdul:
        assert hdul["SPECTRUM"].header["AREASCAL"] == 1
        assert not hdul["SPECTRUM"].header["POISSERR"]

    grouped = work / "group.pi"
    response = work / "group.rmf"
    write_pha(grouped, [0, 0, 10, 10, 10, 10, 10, 10], start=68)
    write_response(response, np.arange(68, 76))
    groups, quality = spectrum.grouping(grouped, response, (0.35, 0.37), 25)
    assert np.array_equal(
        quality == 0, [False, False, True, True, True, True, True, False]
    )
    starts = np.flatnonzero((groups == 1) & (quality == 0))
    assert len(starts) == 1 and groups[starts[0] + 1] == -1
    selected = work / "selected.pi"
    spectrum.write(selected, grouped, spectrum.read(grouped), groups, quality)
    assert spectrum.grouped_range(selected, 350, 374) == (3, 3)
    assert spectrum.grouped_range(selected, 349.999999999, 374.000000001) == (3, 3)
    try:
        spectrum.grouped_range(selected, 350, 369)
    except ValueError:
        pass
    else:
        raise AssertionError("grouped image-band edge accepted")


def test_extraction_stamp(work):
    output = work / "reload"
    (output / "stack").mkdir(parents=True)
    rows = []
    for region, (inner, outer) in pipeline.REGIONS.items():
        for detector in ("pn", "mos1", "mos2"):
            stem = f"{detector}-{region}"
            rows.append(
                {
                    "region": region,
                    "det": detector,
                    "n_frames": 1,
                    "exposure_s": 1,
                    "omega_arcmin2": 1,
                    "response_throughput": 1,
                    "full_area_arcmin2": np.pi * (outer**2 - inner**2) / 3600,
                    "source_counts": 1,
                    "background_counts": 1,
                    "fit_low_keV": pipeline.FIT_BAND_KEV[0],
                    "fit_high_keV": pipeline.FIT_BAND_KEV[1],
                    "group_min_counts": pipeline.PARAMS["group_min_counts"],
                    "channel_max": pipeline.PARAMS["channel_max"][detector],
                    "pn_oot_fraction": pipeline.SHARED["pn_oot_fraction"],
                    "source_pha": f"{stem}-source.pi",
                    "background_pha": f"{stem}-background.pi",
                    "rmf": f"{stem}.rmf",
                    "arf": f"{stem}.arf",
                }
            )
    columns = tuple(rows[0])
    pipeline.write_tsv(output / "spectra.tsv", columns, rows)
    original = pipeline.OUTPUT, inputs.OUTPUT
    pipeline.OUTPUT = inputs.OUTPUT = output
    try:
        assert inputs.frame_detector("mos2-s002-p06") == "mos2"
        frame = {"frame": "pn-p01", "det": "pn"}
        pipeline.write_tsv(
            output / "extraction.tsv",
            ("frame", "det", "region"),
            [{**frame, "region": region} for region in pipeline.REGIONS],
        )
        inputs.validate_extraction_frames([frame])
        assert len(inputs.reload_groups()) == 6
        for bad in (rows[:-1], [rows[0], *rows[:-1]]):
            pipeline.write_tsv(output / "spectra.tsv", columns, bad)
            try:
                inputs.reload_groups()
            except RuntimeError:
                pass
            else:
                raise AssertionError(
                    "incomplete or duplicated spectrum groups accepted"
                )
        pipeline.write_tsv(
            output / "spectra.tsv",
            columns,
            [{**rows[0], "rmf": "wrong.rmf"}, *rows[1:]],
        )
        try:
            inputs.reload_groups()
        except RuntimeError as error:
            assert "filenames" in str(error)
        else:
            raise AssertionError("noncanonical stacked filename accepted")
        row = rows[0]
        changes = {
            "fit_high_keV": row["fit_high_keV"] + 0.01,
            "full_area_arcmin2": row["full_area_arcmin2"] + 1,
            "group_min_counts": row["group_min_counts"] + 1,
            "channel_max": row["channel_max"] + 1,
            "pn_oot_fraction": row["pn_oot_fraction"] + 0.001,
        }
        for key, value in changes.items():
            pipeline.write_tsv(
                output / "spectra.tsv", columns, [{**row, key: value}, *rows[1:]]
            )
            try:
                inputs.reload_groups()
            except RuntimeError as error:
                assert "extraction parameters changed" in str(error)
            else:
                raise AssertionError(f"stale extraction parameter accepted: {key}")
        pipeline.write_tsv(output / "spectra.tsv", columns, rows)
        pipeline.write_tsv(
            output / "extraction.tsv",
            ("frame", "det", "region"),
            [{**frame, "region": "inner"}],
        )
        try:
            inputs.validate_extraction_frames([frame])
        except RuntimeError as error:
            assert "extraction frames changed" in str(error)
        else:
            raise AssertionError("stale extraction frames accepted")
    finally:
        pipeline.OUTPUT, inputs.OUTPUT = original

    old, stamp = work / "old", work / "stamp"
    old.write_text("old")
    stamp.write_text("stamp")
    moment = stamp.stat().st_mtime_ns
    os.utime(old, ns=(moment - 2_000_000_000,) * 2)
    inputs.require_not_newer([old], stamp)
    os.utime(old, ns=(moment + 2_000_000_000,) * 2)
    try:
        inputs.require_not_newer([old], stamp)
    except RuntimeError as error:
        assert "input changed" in str(error)
    else:
        raise AssertionError("newer extraction input accepted")


def test_stack(work):
    rows = []
    for index, (exposure, omega) in enumerate(((100.0, 10.0), (200.0, 20.0))):
        source, background = work / f"s{index}.pi", work / f"b{index}.pi"
        rmf, arf = work / f"r{index}.rmf", work / f"a{index}.arf"
        counts = np.full(80, 10 + index)
        write_pha(source, counts, exposure=exposure, backscal=10 + index)
        write_pha(
            background,
            counts / 10,
            exposure=exposure,
            backscal=10 + index,
            error=np.sqrt(counts) / 10,
        )
        write_response(rmf, np.arange(80))
        shutil.copy(rmf, arf)
        rows.append(
            {
                "source_pha": source,
                "background_pha": background,
                "rmf": rmf,
                "arf": arf,
                "exposure_s": exposure,
                "omega_arcmin2": omega,
            }
        )
    captured = []

    def tool(name, cwd, **parameters):
        lines = (
            (Path(cwd) / parameters["list"].removeprefix("@")).read_text().splitlines()
        )
        captured.append([float(line.split()[-1]) for line in lines])
        source = (Path(cwd) / lines[0].split()[0]).resolve()
        destination = (
            Path(cwd) / parameters["rmffile" if name == "addrmf" else "out_ARF"]
        )
        shutil.copy(source, destination)

    result = spectrum.stack(rows, work / "stack", "mos2-inner", (0.35, 0.40), 25, tool)
    source, background, _, _, exposure, backscal, omega = result
    assert captured == [[0.2, 0.8], [0.2, 0.8]]
    assert exposure == 300 and math.isclose(backscal, 32 / 3)
    assert math.isclose(omega, 50 / 3)
    assert np.allclose(spectrum.read(source).counts, 21)
    assert np.allclose(spectrum.read(background, poisson=False).counts, 2.1)
    with fits.open(source) as hdul:
        header = hdul["SPECTRUM"].header
        assert header["BACKFILE"] == "mos2-inner-background.pi"
        assert {"GROUPING", "QUALITY"} <= set(hdul["SPECTRUM"].data.names)

    weighted = [
        {"exposure_s": 100.0, "omega_arcmin2": 10.0, "response_throughput": 0.5},
        {"exposure_s": 200.0, "omega_arcmin2": 20.0, "response_throughput": 0.25},
    ]
    assert math.isclose(pipeline.stacked_throughput(weighted), 0.3)
    factors = {"inner": {"pn": 2.0, "mos2": 1.0}, "outer": {"pn": 1.2, "mos2": 0.2}}
    groups = [
        {
            "region": "inner",
            "det": "pn",
            "exposure_s": 10.0,
            "omega_arcmin2": 1.0,
            "response_throughput": 0.5,
        },
        {
            "region": "inner",
            "det": "mos2",
            "exposure_s": 10.0,
            "omega_arcmin2": 1.0,
            "response_throughput": 0.25,
        },
        {
            "region": "outer",
            "det": "pn",
            "exposure_s": 10.0,
            "omega_arcmin2": 1.0,
            "response_throughput": 0.1,
        },
        {
            "region": "outer",
            "det": "mos2",
            "exposure_s": 10.0,
            "omega_arcmin2": 1.0,
            "response_throughput": 0.3,
        },
    ]
    corrected = pipeline.on_axis_factors(factors, groups)
    assert all(
        math.isclose(corrected[region][detector], expected)
        for region, detector, expected in (
            ("inner", "pn", 4),
            ("inner", "mos2", 4),
            ("outer", "pn", 12),
            ("outer", "mos2", 2 / 3),
        )
    )
    image = pipeline.image_factors(corrected, groups)
    assert math.isclose(image["inner"], 4) and math.isclose(image["outer"], 3.5)


def response_metrics(stacked, direct):
    stacked, direct = np.asarray(stacked, float), np.asarray(direct, float)
    epsilon = stacked.sum() / direct.sum() - 1
    shape_l1 = np.abs(stacked / stacked.sum() - direct / direct.sum()).sum()
    return epsilon, shape_l1


def weighted_response(spectra, weights):
    return np.average(
        np.asarray(spectra, float), axis=0, weights=np.asarray(weights, float)
    )


def test_response_arithmetic():
    spectra = np.array([[4.0, 0.0, 0.0], [0.0, 1.0, 3.0]])
    weights = np.array([1.0, 4.0])
    expected = (spectra[0] + 4 * spectra[1]) / 5
    assert np.allclose(weighted_response(spectra, weights), expected)
    epsilon, shape_l1 = response_metrics(
        weighted_response(spectra, weights[::-1]), expected
    )
    assert abs(epsilon) > RATE_TOLERANCE or shape_l1 > SHAPE_TOLERANCE


def response_inputs():
    frames = []
    for row in inputs.read_tsv(inputs.OUTPUT / "extraction.tsv"):
        directory = inputs.WORK / "frames" / f"{row['frame']}-{row['region']}"
        frames.append(
            {
                "frame": row["frame"],
                "region": row["region"],
                "det": row["det"],
                "weight": float(row["exposure_s"]) * float(row["omega_arcmin2"]),
                "source_pha": directory / "source.pi",
                "rmf": directory / "source.rmf",
                "arf": directory / "source.arf",
            }
        )
    groups = inputs.reload_groups()
    for item in frames + groups:
        for name in ("source_pha", "rmf", "arf"):
            if not Path(item[name]).is_file():
                raise FileNotFoundError(item[name])
    return frames, groups


def check_response_grid(item):
    with fits.open(item["source_pha"]) as pha, fits.open(item["rmf"]) as rmf:
        channels = np.asarray(pha["SPECTRUM"].data["CHANNEL"])
        bounds = rmf["EBOUNDS"].data
        response_channels = np.asarray(bounds["CHANNEL"])
        low = np.asarray(bounds["E_MIN"], float)
        high = np.asarray(bounds["E_MAX"], float)
    if not np.array_equal(channels, response_channels):
        raise RuntimeError(f"PHA/RMF channels differ for {item['source_pha']}")
    use = (high > inputs.FIT_BAND_KEV[0]) & (low <= inputs.FIT_BAND_KEV[1])
    indices = np.flatnonzero(use)
    if (
        not len(indices)
        or np.any(np.diff(indices) != 1)
        or low[indices[0]] > inputs.FIT_BAND_KEV[0]
        or high[indices[-1]] < inputs.FIT_BAND_KEV[1]
        or np.any(low[indices[1:]] > high[indices[:-1]] + 1e-6)
    ):
        raise RuntimeError(f"response does not cover fit band: {item['source_pha']}")


def folded_response(xspec, item, settings, shape, bins=None):
    check_response_grid(item)
    previous = Path.cwd()
    xspec.AllData.clear()
    xspec.AllModels.clear()
    os.chdir(Path(item["source_pha"]).parent)
    try:
        data = xspec.Spectrum(Path(item["source_pha"]).name)
        data.response = Path(item["rmf"]).name
        data.response.arf = Path(item["arf"]).name
        xspec.AllData.ignore("bad")
        xspec.AllData.ignore(
            f"**-{inputs.FIT_BAND_KEV[0]:g} {inputs.FIT_BAND_KEV[1]:g}-**"
        )
        fitted = xspec.Model(model.EXPRESSION)
        model.configure(fitted, item["det"], 1, settings)
        _, camera, lhb, _, halo, cxb, cx = model.components(fitted)
        model.fixed(camera.factor, 1)
        for component in (lhb, halo, cxb):
            model.fixed(component.norm, 0)
        model.fixed(cx.temperature, shape[0])
        model.fixed(cx.collnpar, shape[1])
        model.fixed(cx.norm, 1)
        xspec.Plot("data")
        values = np.asarray(xspec.Plot.model(1), float)
        energy = np.asarray(xspec.Plot.x(1), float)
        half_width = np.asarray(xspec.Plot.xErr(1), float)
        low, high = energy - half_width, energy + half_width
        if not len(values) or np.any(low[1:] > high[:-1] + 1e-6):
            raise RuntimeError(f"invalid folded bins: {item['source_pha']}")
        if bins is None:
            bins = np.column_stack((low, high))
        else:
            low, high = np.maximum(low, bins[0, 0]), np.minimum(high, bins[-1, 1])
            use = high > low
            values, low, high = values[use], low[use], high[use]
        if np.any(bins[:, 1] <= bins[:, 0]) or np.any(
            bins[1:, 0] > bins[:-1, 1] + 1e-6
        ):
            raise RuntimeError("invalid common response bins")
        output = np.zeros(len(bins))
        coverage = np.zeros(len(bins))
        for value, left, right in zip(values, low, high):
            overlap = np.maximum(
                0, np.minimum(right, bins[:, 1]) - np.maximum(left, bins[:, 0])
            )
            output += value * overlap
            coverage += overlap
        native = np.sum(values * (high - low))
        if (
            not np.isfinite(output).all()
            or output.sum() <= 0
            or not math.isclose(output.sum(), native, rel_tol=1e-6)
            or not np.allclose(coverage, bins[:, 1] - bins[:, 0], rtol=0, atol=1e-6)
        ):
            raise RuntimeError(
                f"response does not cover common bins: {item['source_pha']}"
            )
        return output, bins
    finally:
        xspec.AllData.clear()
        xspec.AllModels.clear()
        os.chdir(previous)


def validate_response_mixture(xspec, frames, groups, settings, shape):
    report = []
    order = {"inner": 0, "outer": 1, "pn": 0, "mos1": 1, "mos2": 2}
    for group in sorted(
        groups, key=lambda row: (order[row["region"]], order[row["det"]])
    ):
        members = [
            row
            for row in frames
            if row["region"] == group["region"] and row["det"] == group["det"]
        ]
        weights = np.asarray([row["weight"] for row in members])
        if len(members) != group["n_frames"] or np.any(weights <= 0):
            raise RuntimeError(
                f"invalid response members: {group['det']}/{group['region']}"
            )
        stacked, bins = folded_response(xspec, group, settings, shape)
        direct = weighted_response(
            [folded_response(xspec, row, settings, shape, bins)[0] for row in members],
            weights,
        )
        epsilon, shape_l1 = response_metrics(stacked, direct)
        if abs(epsilon) > RATE_TOLERANCE or shape_l1 > SHAPE_TOLERANCE:
            raise RuntimeError(
                f"stacked response differs: {group['det']}/{group['region']}"
            )
        report.append(
            (
                group["region"],
                group["det"],
                len(members),
                direct.sum(),
                stacked.sum(),
                epsilon,
                shape_l1,
            )
        )
    return report


def assert_relative(actual, expected, tolerance=0.01):
    if not math.isclose(
        float(actual),
        float(expected),
        rel_tol=tolerance,
        abs_tol=tolerance * max(abs(float(expected)), 1e-30),
    ):
        raise AssertionError(f"{actual} does not reproduce {expected}")


def validate_fit_replay(groups, settings, result):
    saved = result["fit"]
    replay = fit.fit(
        groups,
        settings,
        result["distance_au"],
        inputs.SOFTWARE["atomdb"],
        tuple(result["fit_band_keV"]),
    )
    if replay["n_free_parameters"] != 9 or replay["dof"] != saved["dof"]:
        raise AssertionError("nine-parameter fit contract changed")
    if abs(replay["statistic"] - saved["chi2"]) > 0.1:
        raise AssertionError("saved fit statistic did not replay")
    for name, expected in saved["parameters"].items():
        assert_relative(replay["parameters"][name], expected)
    for detector in ("mos1", "mos2"):
        assert_relative(replay["crosscal"][detector], saved["crosscal"][detector])
    for region in ("inner", "outer"):
        for name in ("cx_norm_per_arcmin2", "flux_erg_cm2_s", "luminosity_erg_s"):
            assert_relative(
                replay["regions"][region][name], result["regions"][region][name]
            )
    factors = pipeline.on_axis_factors(replay["region_detector_K"], groups)["inner"]
    for detector, expected in result["detector_K_ct_s_per_erg_cm2_s"].items():
        assert_relative(factors[detector], expected)
    assert_relative(replay["mean_photon_energy_eV"], result["mean_photon_energy_eV"])
    assert_relative(1 / factors["mos2"], result["mos2_ecf_erg_cm2_per_count"])

    displaced = json.loads(json.dumps(settings))
    displaced["cx_temperature_start_keV"] = 0.75
    displaced["cx_collision_starts_keV_u"] = [80.0, 60.0]
    alternative = fit.fit(
        groups,
        displaced,
        result["distance_au"],
        inputs.SOFTWARE["atomdb"],
        tuple(result["fit_band_keV"]),
    )
    if alternative["statistic"] < saved["chi2"] - 0.1:
        raise AssertionError("displaced seed found a lower fit")


def test_products():
    os.environ["ATOMDB"] = inputs.SOFTWARE["atomdb"]
    import xspec

    __import__("acx2_xspec")
    xspec.Xset.chatter = xspec.Xset.logChatter = 0
    xspec.Xset.abund, xspec.Xset.xsect = "wilm", "vern"
    xspec.Plot.device, xspec.Plot.xAxis = "/null", "keV"
    xspec.Plot.background, xspec.Plot.add = False, False
    result = json.loads((inputs.OUTPUT / "result.json").read_text())
    interval = result["fit"]["profile_intervals_1sigma"]["crosscal_mos2"]
    expected_sigma = pipeline.ecf_sigma_ln(result["fit"]["crosscal"]["mos2"], interval)
    assert math.isclose(result["mos2_ecf_sigma_ln"], expected_sigma, rel_tol=1e-12)
    assert "products" not in result
    settings = inputs.PARAMS["model"]
    frames, groups = response_inputs()
    parameters = result["fit"]["parameters"]
    shape = (parameters["cx_temperature_keV"], parameters["cx_collision_keV_u"])
    validate_response_mixture(xspec, frames, groups, settings, shape)
    validate_fit_replay(groups, settings, result)


def test_model_contract():
    settings = json.loads((HERE / "parameters.json").read_text())["model"]
    groups = [
        {"region": region, "det": detector}
        for region in ("inner", "outer")
        for detector in model.DETECTORS
    ]
    models = [Model() for _ in groups]
    for group, fitted in zip(groups, models):
        model.configure(fitted, group["det"], 100, settings)
    _, regions = model.link(groups, models)
    assert len(model.free_parameters(models)) == 9
    cx = model.components(models[0])[-1]
    assert cx.recombtype.values[0] == 1 and cx.acxmodel.values[0] == 8
    assert cx.collntype.values[0] == 1
    assert cx.recombtype.frozen and cx.acxmodel.frozen and cx.collntype.frozen
    assert cx.Hefrac.values[0] == settings["cx_donor_helium_fraction"]
    assert model.components(models[0])[0].factor.values[0] == 100
    assert all(
        getattr(cx, name).values[0] == 1
        for name in (
            "H",
            "He",
            "C",
            "N",
            "O",
            "Ne",
            "Mg",
            "Al",
            "Si",
            "S",
            "Ar",
            "Ca",
            "Fe",
            "Ni",
        )
    )
    assert set(regions) == {"inner", "outer"}
    assert math.isclose(fit.mean_energy_eV(5.024e-10, 1), 313.57, rel_tol=1e-3)
    expected = math.hypot(math.log1p(0.10), max(-math.log(0.99), math.log(1.02)))
    assert math.isclose(pipeline.ecf_sigma_ln(1, [0.99, 1.02]), expected)
    for crosscal, interval in ((0, [0.9, 1.1]), (1, [1.1, 0.9]), (1, [0, 1.1])):
        try:
            pipeline.ecf_sigma_ln(crosscal, interval)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid ECF interval accepted")

    pn_groups = [{"region": region, "det": "pn"} for region in ("inner", "outer")]
    pn_models = [Model() for _ in pn_groups]
    for group, fitted in zip(pn_groups, pn_models):
        model.configure(fitted, "pn", 100, settings)
    detectors, regions = model.link(pn_groups, pn_models)
    targets = fit.profile_targets(pn_models, detectors, regions)
    assert not any(name.startswith("crosscal_") for name in targets)
    assert {"cx_norm_inner", "cx_norm_outer"} <= set(targets)

    inner_groups = [
        {"region": "inner", "det": detector} for detector in model.DETECTORS
    ]
    inner_models = [Model() for _ in inner_groups]
    for group, fitted in zip(inner_groups, inner_models):
        model.configure(fitted, group["det"], 100, settings)
    detectors, regions = model.link(inner_groups, inner_models)
    targets = fit.profile_targets(inner_models, detectors, regions)
    assert {"crosscal_mos1", "crosscal_mos2", "cx_norm_inner"} <= set(targets)
    assert "cx_norm_outer" not in targets


def test_profile_failure():
    fitted = Model()
    fitted.startParIndex = 1
    target = {"cx": (fitted, fitted.vacx2.temperature)}
    fake = SimpleNamespace(
        AllModels=lambda index, name: fitted,
        AllData=SimpleNamespace(nGroups=1),
        Fit=SimpleNamespace(error=lambda command: (_ for _ in ()).throw(ValueError())),
    )
    original = fit.minimise
    fit.minimise = lambda xspec: 0
    try:
        try:
            fit.profile_errors(fake, target)
        except RuntimeError as error:
            assert "profile errors failed" in str(error)
        else:
            raise AssertionError("failed profile errors accepted")
    finally:
        fit.minimise = original


def test_conversion_regions():
    groups = [
        {
            "region": region,
            "det": detector,
            "omega_arcmin2": area,
            "source_pha": Path("unused"),
        }
        for (region, detector), area in zip(
            (("inner", "pn"), ("inner", "mos2"), ("outer", "pn"), ("outer", "mos2")),
            (2.0, 3.0, 4.0, 5.0),
        )
    ]
    cameras = (1.0, 0.8, 1.0, 0.8)
    raw_k = (5.0, 4.0, 3.0, 2.0)
    intrinsic = [
        (group["omega_arcmin2"] * camera * 2, group["omega_arcmin2"] * camera)
        for group, camera in zip(groups, cameras)
    ]
    rates = [
        group["omega_arcmin2"] * camera * factor * 2
        for group, camera, factor in zip(groups, cameras, raw_k)
    ]
    models = [
        SimpleNamespace(
            componentNames=["area", "camera"],
            area=None,
            camera=SimpleNamespace(factor=SimpleNamespace(values=[camera])),
        )
        for camera in cameras
    ]
    calls = []

    def all_data(index):
        return SimpleNamespace(
            notice=lambda value: calls.append((index, "notice", value)),
            ignore=lambda value: calls.append((index, "ignore", value)),
        )

    all_data.notice = lambda value: calls.append(("all", "notice", value))
    all_data.ignore = lambda value: calls.append(("all", "ignore", value))
    xspec = SimpleNamespace(AllData=all_data)
    mismatch = [False]
    old_isolated, old_fluxes, old_plot = fit.isolated, fit.fluxes, fit.plot_arrays
    old_range = fit.spectrum.grouped_range
    fit.isolated = lambda *args: contextlib.nullcontext()
    fit.fluxes = lambda xspec, count, band: intrinsic
    fit.plot_arrays = lambda xspec, count: [
        {
            "energy_keV": np.arange(2 if mismatch[0] else 1),
            "model_ct_s_keV": np.array([rate]),
            "half_width_keV": np.array([0.5]),
        }
        for rate in rates
    ]
    fit.spectrum.grouped_range = lambda *args: (2, 2)
    try:
        factors, mean, intrinsic = fit.conversion(xspec, models, {}, groups, (0.4, 0.9))
        assert calls[-2:] == [("all", "notice", "all"), ("all", "ignore", "bad")]
        assert [call for call in calls if call[1] == "ignore"][:-1] == [
            (index, "ignore", "**-1 3-**") for index in range(1, 5)
        ]
        mismatch[0] = True
        try:
            fit.conversion(xspec, models, {}, groups, (0.4, 0.9))
        except RuntimeError:
            pass
        else:
            raise AssertionError("XSPEC image-band mismatch accepted")
        assert calls[-2:] == [("all", "notice", "all"), ("all", "ignore", "bad")]
    finally:
        fit.isolated, fit.fluxes, fit.plot_arrays = old_isolated, old_fluxes, old_plot
        fit.spectrum.grouped_range = old_range
    assert factors == {
        "inner": {"pn": 5.0, "mos2": 3.2},
        "outer": {"pn": 3.0, "mos2": 1.6},
    }
    assert math.isclose(mean, fit.mean_energy_eV(*intrinsic[0]))
    assert fit.band_text(pipeline.FIT_BAND_KEV) == "0.35 1.1"
    assert pipeline.FIT_BAND_KEV == tuple(
        value / 1000 for value in pipeline.SHARED["soft_band_ev"]
    )
    assert pipeline.REGIONS == {
        "inner": (0, float(pipeline.SHARED["science_aperture_radius_arcsec"])),
        "outer": tuple(map(float, pipeline.SHARED["spectrum_outer_annulus_arcsec"])),
    }


def main():
    with tempfile.TemporaryDirectory(dir=HERE) as temporary:
        work = Path(temporary)
        test_geometry()
        test_sas_environment(work)
        test_center_and_selection(work)
        test_map_fraction(work)
        test_response_arguments(work)
        test_background_and_grouping(work)
        test_extraction_stamp(work)
        test_stack(work)
        test_response_arithmetic()
        test_model_contract()
        test_profile_failure()
        test_conversion_regions()


if __name__ == "__main__":
    if sys.argv[1:] == ["products"]:
        warnings.filterwarnings("ignore", category=ResourceWarning)
        test_products()
    elif not sys.argv[1:]:
        main()
    else:
        raise SystemExit("usage: test.py [products]")
