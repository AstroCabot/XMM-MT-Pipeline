import re
import shutil
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

import background
import common
import comove
import events
import figure
import flare
import pha
import prune
import reprocess
import soft_protons
import sources

HERE = Path(__file__).resolve().parent
REPRO = HERE.parent.resolve()
WORK = (REPRO / "work" / "reduction_tests").resolve()
FLARE_NEGATIVE = (slice(120, 160), slice(220, 260))
assert common.CFG["python_lib"] == "/home/scabot/miniconda3/envs/xmm/lib"


def raises(error, function, *args, **kwargs):
    try:
        function(*args, **kwargs)
    except error:
        return
    raise AssertionError(f"{function.__name__} did not raise {error.__name__}")


def clear_work():
    expected = (REPRO / "work").resolve()
    assert WORK.parent == expected and WORK.name == "reduction_tests"
    if WORK.exists():
        shutil.rmtree(WORK)


def write_events(
    path,
    x,
    y,
    pi=None,
    pattern=None,
    flag=None,
    ccd=None,
    time=None,
    detx=None,
    dety=None,
):
    x, y = np.asarray(x, float), np.asarray(y, float)
    size = len(x)
    pi = np.full(size, 500) if pi is None else np.asarray(pi)
    pattern = np.zeros(size, int) if pattern is None else np.asarray(pattern)
    flag = np.zeros(size, int) if flag is None else np.asarray(flag)
    ccd = np.ones(size, int) if ccd is None else np.asarray(ccd)
    time = np.arange(size, dtype=float) if time is None else np.asarray(time, float)
    detx = x.copy() if detx is None else np.asarray(detx, float)
    dety = y.copy() if dety is None else np.asarray(dety, float)
    columns = [
        fits.Column(name="X", format="D", array=x),
        fits.Column(name="Y", format="D", array=y),
        fits.Column(name="DETX", format="D", array=detx),
        fits.Column(name="DETY", format="D", array=dety),
        fits.Column(name="TIME", format="D", array=time),
        fits.Column(name="PI", format="J", array=pi),
        fits.Column(name="PATTERN", format="J", array=pattern),
        fits.Column(name="FLAG", format="J", array=flag),
        fits.Column(name="CCDNR", format="J", array=ccd),
        fits.Column(name="RAWY", format="J", array=np.full(size, 20)),
    ]
    hdu = fits.BinTableHDU.from_columns(columns, name="EVENTS")
    hdu.header.update(
        MJDREF=60000.0,
        TIMESYS="UTC",
        FILTER="Thin1",
        INSTRUME="EMOS1",
        TSTART=float(time.min()),
        TSTOP=float(time.max()),
        REFXCRVL=100.0,
        REFYCRVL=20.0,
        REFXCRPX=0.0,
        REFYCRPX=0.0,
        REFXCDLT=-0.05 / 3600,
        REFYCDLT=0.05 / 3600,
        REFXCTYP="RA---TAN",
        REFYCTYP="DEC--TAN",
    )
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)


def write_canonical(path, states=None):
    states = states or {}
    header = fits.Header({f"ANOMFL{ccd}": states.get(ccd, "U") for ccd in range(2, 8)})
    fits.PrimaryHDU(header=header).writeto(path)


def write_gti(path, spans):
    spans = np.asarray(spans, float)
    columns = [
        fits.Column(name="START", format="D", array=spans[:, 0]),
        fits.Column(name="STOP", format="D", array=spans[:, 1]),
    ]
    fits.HDUList(
        [fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns, name="STDGTI")]
    ).writeto(path)


def gti_hdu(spans, name):
    spans = np.asarray(spans, float)
    columns = [
        fits.Column(name="START", format="D", array=spans[:, 0]),
        fits.Column(name="STOP", format="D", array=spans[:, 1]),
    ]
    return fits.BinTableHDU.from_columns(columns, name=name)


def write_pha(
    path, channel, values, exposure=100.0, backscal=1.0, field="COUNTS", errors=None
):
    columns = [
        fits.Column(name="CHANNEL", format="I", array=channel),
        fits.Column(name=field, format="E", array=np.asarray(values, "f4")),
    ]
    if errors is not None:
        columns.append(
            fits.Column(name="STAT_ERR", format="E", array=np.asarray(errors, "f4"))
        )
    hdu = fits.BinTableHDU.from_columns(columns, name="SPECTRUM")
    hdu.header.update(EXPOSURE=exposure, BACKSCAL=backscal, SPECDELT=5.0)
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)


def write_flare_ccf(path):
    instruments, bands, filters, maps = [], [], [], []
    for instrument in ("EMOS1", "EMOS2", "EPN"):
        for band_index, band in enumerate(("B4", "B5", "B6"), 1):
            for filter_index, filter_name in enumerate(("Thin1", "Medium", "Thick"), 1):
                instruments.append(instrument)
                bands.append(band)
                filters.append(filter_name)
                flare_map = np.full(
                    (780, 780), band_index * 10 + filter_index, dtype=np.float32
                )
                flare_map[FLARE_NEGATIVE] *= -1  # Real FLARE maps are signed.
                maps.append(flare_map)
    columns = [
        fits.Column(name="INSTRUME", format="5A", array=instruments),
        fits.Column(name="BAND", format="2A", array=bands),
        fits.Column(name="FILTER", format="6A", array=filters),
        fits.Column(name="FLAREMAP", format="608400E", dim="(780,780)", array=maps),
    ]
    fits.HDUList(
        [fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns, name="FLAREMAP")]
    ).writeto(path)


def image_wcs(pixel_arcsec=4.0, size=101):
    wcs = WCS(naxis=2)
    wcs.wcs.crval = [100.0, 20.0]
    wcs.wcs.crpix = [(size + 1) / 2] * 2
    wcs.wcs.cdelt = [-pixel_arcsec / 3600, pixel_arcsec / 3600]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return wcs


def test_output_boundary():
    original = common.OUTPUT
    try:
        common.OUTPUT = (WORK / "outputs").resolve()
        assert common.setup_output() == common.OUTPUT
        raises(RuntimeError, common.directory, "../../../escape")
        common.OUTPUT = (REPRO.parent / "escape").resolve()
        raises(RuntimeError, common.setup_output)
    finally:
        common.OUTPUT = original
    table = WORK / "table.tsv"
    common.write_tsv(table, ("a", "b"), [{"a": 1, "b": 2}])
    assert common.read_tsv(table) == [{"a": "1", "b": "2"}]


def test_gti():
    event = WORK / "gti-events.fits"
    write_events(event, [0, 1, 2], [0, 1, 2], time=[100, 150, 200])
    all_pass, valid, outside = (
        WORK / name for name in ("all-pass.fits", "valid.fits", "outside.fits")
    )
    write_gti(all_pass, [[0, 200000]])
    write_gti(valid, [[120, 180]])
    write_gti(outside, [[300, 400]])
    assert not reprocess.real_gti(all_pass, event)
    assert reprocess.real_gti(valid, event)
    assert not reprocess.real_gti(outside, event)

    with fits.open(event, mode="update") as hdul:
        hdul["EVENTS"].header["ONTIME"] = 20.0
        hdul.append(gti_hdu([[0, 50]], "STDGTI"))
        hdul.append(gti_hdu([[5, 10], [20, 25], [30, 40]], "GTI00003"))
    frame = {"frame": "synthetic", "event": event, "ontime": "20.0"}
    times, weights = comove.samples(frame)
    assert np.isclose(weights.sum(), 20)
    assert len(times) == 4 and np.all(weights == 5)


def test_event_support_and_regions():
    canonical = WORK / "mos1-canonical.fits"
    write_canonical(canonical, {2: "B", 3: "O", 4: "I", 6: "O"})
    path = WORK / "selection.fits"
    write_events(
        path,
        np.arange(8),
        np.arange(8),
        pi=[500, 500, 500, 500, 3000, 3000, 3000, 9000],
        pattern=[0, 13, 0, 0, 12, 12, 13, 0],
        flag=[0, 0, 1, 0, 0, 0, 0, 0],
        ccd=[1, 1, 1, 2, 2, 3, 1, 1],
    )
    data = fits.getdata(path, "EVENTS")
    assert np.flatnonzero(
        events.band_mask(data, "mos1", "soft", events.anomaly_drops(canonical, "soft"))
    ).tolist() == [0]
    assert np.flatnonzero(
        events.band_mask(data, "mos1", "hard", events.anomaly_drops(canonical, "hard"))
    ).tolist() == [4]
    assert events.anomaly_drops(canonical, "soft") == {2, 3, 4, 6}
    assert events.anomaly_drops(canonical, "hard") == {3, 6}
    assert events.FLAG_EXPRESSION == "FLAG==0"
    audit = events.selection_rows(
        [
            {
                "exposure": "mos1-test",
                "det": "mos1",
                "frame": "test-p01",
                "pointing": 1,
                "ontime": "420.0",
                "soft": "800",
            }
        ],
        set(),
    )
    measured = next(row for row in audit if row["exposure"] == "mos1-test")
    assert (
        measured["decision"] == "real espfilt GTI; maximum frame ONTIME 420 s < 500 s"
    )
    frame = {"det": "mos1", "exposure": "mos1-s001"}
    assert background.chip_string(frame, "hard", canonical) == "T T F T T F T"

    axis = np.linspace(-2000, 2000, 30)
    x, y = np.meshgrid(axis, axis)
    x, y = x.ravel(), y.ravel()
    angle = 0.3
    detx = np.cos(angle) * x - np.sin(angle) * y + 30
    dety = np.sin(angle) * x + np.cos(angle) * y - 20
    region_event = WORK / "region-events.fits"
    write_events(region_event, x, y, detx=detx, dety=dety)
    source_xy = np.array([500.0, -750.0])
    ra, dec = sources.event_wcs(region_event).wcs_pix2world(*source_xy, 1)
    source = {"ra": str(ra), "dec": str(dec), "radius_arcsec": "45"}
    region = sources.detector_region(region_event, [source])
    match = re.search(r"circle\(([-0-9.]+),([-0-9.]+),([-0-9.]+)\)", region)
    expected = np.array(
        [
            np.cos(angle) * source_xy[0] - np.sin(angle) * source_xy[1] + 30,
            np.sin(angle) * source_xy[0] + np.cos(angle) * source_xy[1] - 20,
        ]
    )
    assert match and np.allclose([float(match[1]), float(match[2])], expected, atol=0.1)
    assert float(match[3]) == 900.0
    assert np.hypot(7131.2 - 97, 18131.8 + 172) > 17152 + 900
    assert sources.FOV["EMOS2"] == (-61, -228, 17085)
    assert sources.intersects_fov("EPN", -15474.4, -929.3, 900)
    assert not sources.intersects_fov("EPN", 15385.7, -598.4, 900)
    far_ra, far_dec = sources.event_wcs(region_event).wcs_pix2world(1e6, 1e6, 1)
    far = {"ra": str(far_ra), "dec": str(far_dec), "radius_arcsec": "45"}
    assert sources.detector_region(region_event, [far]) == ""
    sky_region = sources.source_expression(region_event, [source])
    assert "circle(500.000,-750.000,900.000)" in sky_region
    distorted = WORK / "distorted-events.fits"
    write_events(distorted, x, y, detx=2 * x, dety=y)
    raises(RuntimeError, sources.detector_region, distorted, [source])

    chip_event = WORK / "chip-events.fits"
    write_events(
        chip_event, [0.5, 0.5], [0.5, 0.5], pi=[500, 500], ccd=[1, 5], time=[10, 20]
    )
    out_wcs, shape = comove.output_wcs(100, 20)
    eph = {
        "mjd": np.array([59999.99, 60000.01]),
        "ra_deg": np.array([99.99, 100.01]),
        "dec_deg": np.array([20.0, 20.0]),
    }
    image = comove.event_image(
        chip_event, "mos1", "soft", "T T T T F T T", out_wcs, shape, eph
    )
    assert image.sum() == 1


def test_pha():
    assert soft_protons.safe_index(1) == 1.001
    assert soft_protons.safe_index(0.8) == 0.8
    response = WORK / "diagonal.rsp"
    bounds = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="CHANNEL", format="J", array=[1, 2]),
            fits.Column(name="E_MIN", format="E", array=[0.0, 0.006]),
            fits.Column(name="E_MAX", format="E", array=[0.006, 0.010]),
        ],
        name="EBOUNDS",
    )
    fits.HDUList([fits.PrimaryHDU(), bounds]).writeto(response)
    template = pha.Spectrum(np.array([0, 1]), np.zeros(2), np.zeros(2), 1.0, 1.0, 5.0)
    rebinned = pha.response_counts(np.array([6.0, 4.0]), response, template)
    assert np.allclose(rebinned, [5.0, 5.0]) and rebinned.sum() == 10
    channel = np.array([70, 100, 150, 220, 600])
    counts = np.array([10, 20, 30, 40, 50.0])
    source, qpb = (WORK / name for name in ("source.pi", "qpb.pi"))
    write_pha(source, channel, counts, exposure=100, backscal=0.8)
    write_pha(qpb, channel, [1, 2, 3, 4, 5], exposure=100, backscal=0.8)
    read = pha.read(source)
    assert np.array_equal(read.counts, counts)
    assert np.array_equal(read.variance, counts)
    unused = WORK / "unused-error-negative.pi"
    write_pha(unused, [0, 1], [0, 1], errors=[-123, 1])
    assert np.array_equal(pha.read(unused).variance, [0, 1])
    unused = WORK / "unused-error-positive.pi"
    write_pha(unused, [0, 1], [0, 1], errors=[1e30, 1])
    assert np.array_equal(pha.read(unused).variance, [0, 1])
    invalid = WORK / "invalid-error.pi"
    write_pha(invalid, [0, 1], [0, 1], errors=[0, -1])
    raises(ValueError, pha.read, invalid)
    net = pha.net(source, qpb)
    expected = counts - np.array([1, 2, 3, 4, 5])
    variance = counts + np.array([1, 2, 3, 4, 5])
    assert np.allclose(net.counts, expected) and np.allclose(net.variance, variance)
    bad_norm = WORK / "bad-norm.pi"
    write_pha(bad_norm, channel, counts, exposure=101, backscal=0.8)
    raises(ValueError, pha.net, source, bad_norm)
    other = net._replace(
        counts=2 * net.counts, variance=2 * net.variance, exposure=200, backscal=0.5
    )
    stacked = pha.stack([net, other])
    assert np.allclose(stacked.counts, 3 * expected)
    assert np.allclose(stacked.variance, 3 * variance)
    assert stacked.exposure == 300 and np.isclose(stacked.backscal, 0.6)
    assert np.array_equal(pha.grouping([10, 10, 10, 10]), [1, -1, -1, -1])
    assert np.array_equal(pha.grouping([30, 30, 30]), [1, 1, 1])
    assert np.isclose(pha.band_counts(read, 350, 1100), 100)
    assert np.isclose(pha.area_arcmin2(read), 0.8 * (0.05 / 60) ** 2)
    scale, norm = soft_protons.frame_norm(read, {"hard_net_rate": 1.0, "norm": 2.0})
    assert np.isclose(scale, 0.5) and np.isclose(norm, 1.0)

    output = WORK / "stacked.pi"
    pha.write(output, source, stacked, rate=True, group_counts=stacked.counts)
    reread = pha.read(output)
    assert np.allclose(reread.counts, stacked.counts, rtol=2e-7)
    assert np.allclose(reread.variance, stacked.variance, rtol=2e-7)
    with fits.open(output) as hdul:
        assert len(hdul) == 2
        assert "RATE" in hdul["SPECTRUM"].data.dtype.names
        assert "GROUPING" in hdul["SPECTRUM"].data.dtype.names
        assert hdul["SPECTRUM"].header["EXPOSURE"] == 300

    rate = WORK / "rate.pi"
    write_pha(
        rate,
        channel,
        counts / 100,
        exposure=100,
        field="RATE",
        errors=np.sqrt(counts) / 100,
    )
    rate_read = pha.read(rate)
    assert np.allclose(rate_read.counts, counts)
    assert np.allclose(rate_read.variance, counts)
    raises(ValueError, pha.stack, [])
    raises(ValueError, pha.stack, [read, read._replace(channel=channel + 1)])
    raises(ValueError, pha.stack, [read, read._replace(channel_ev=10)])


def test_hard_flare_correction():
    calibration = WORK / "flare.ccf"
    write_flare_ccf(calibration)
    expected = sum(
        weight * sum(band_index * 10 + filter_index for filter_index in (1, 2, 3))
        for band_index, (_, weight) in enumerate(flare.HARD_WEIGHTS, 1)
    )
    signed = np.full((780, 780), float(expected))
    signed[FLARE_NEGATIVE] = -expected
    template = flare.hard_template("mos1", calibration)
    assert template.shape == (780, 780) and np.allclose(template, signed)

    detector = WORK / "hard-detector.fits"
    support = WORK / "hard-support.fits"
    original = np.zeros((780, 780), np.float32)
    original[20:30, 40:60] = 2
    mask = np.zeros((780, 780), np.float32)
    mask[100:500, 200:600] = 1
    fits.PrimaryHDU(original).writeto(detector)
    fits.PrimaryHDU(mask).writeto(support)
    target = float(original.sum(dtype=np.float64))
    flare.replace_hard_map(detector, "mos1", "T T F T T F T", support, calibration)
    corrected, header = fits.getdata(detector, header=True)
    assert np.isclose(corrected.sum(dtype=np.float64), target, rtol=2e-7)
    assert np.array_equal(corrected != 0, mask != 0)
    assert header["SPCORR"] == "CCF FLARE B4-B6" and header["SPFILT"] == "ALL"
    assert np.isclose(header["SPB4W"], 1499 / 1999)
    # Signed template pixels must survive scaling and chip masking.
    kept = signed * (mask != 0)
    assert np.allclose(
        corrected, kept * (target / kept.sum(dtype=np.float64)), rtol=2e-6, atol=1e-9
    )
    negative = np.zeros((780, 780), bool)
    negative[FLARE_NEGATIVE] = True
    assert np.array_equal(corrected < 0, negative)

    unmasked = WORK / "hard-unmasked.fits"
    fits.PrimaryHDU(original).writeto(unmasked)
    flare.replace_hard_map(unmasked, "pn", "T T T T", calibration=calibration)
    assert np.count_nonzero(fits.getdata(unmasked)) == 780 * 780
    assert np.array_equal(fits.getdata(unmasked) < 0, negative)


def test_comove_maps():
    expected = Time(comove.XMM_MJDREF, format="mjd", scale="tt").utc.mjd
    assert np.isclose(comove.mjd(np.array([0.0]), fits.Header())[0], expected)
    wrap = {
        "mjd": np.array([0.0, 1.0]),
        "ra_deg": np.array([359.9, 0.1]),
        "dec_deg": np.array([10.0, 12.0]),
    }
    ra, dec = comove.position(np.array([0.5]), wrap)
    assert min(abs(ra[0]), abs(ra[0] - 360)) < 1e-10 and dec[0] == 11
    out_wcs, shape = comove.output_wcs(100, 20)
    assert shape == (751, 751)
    assert np.allclose(out_wcs.world_to_pixel_values(100, 20), [375, 375])
    deposited = comove.deposit(
        (10, 10), np.array([2.2, 7.4]), np.array([3.1, 6.8]), np.array([4, 9])
    )
    assert np.isclose(deposited.sum(), 13)

    sky_wcs = image_wcs(2, 21)
    data = np.zeros((21, 21))
    data[10, 10] = 100
    path = WORK / "sky-map.fits"
    fits.PrimaryHDU(data, sky_wcs.to_header()).writeto(path)
    times = np.arange(5, 100, 10.0)
    weights = np.full(len(times), 10.0)
    header = fits.Header({"MJDREF": 60000.0, "TIMESYS": "UTC"})
    eph = {
        "mjd": np.array([59999.99, 60000.01]),
        "ra_deg": np.array([99.99, 100.01]),
        "dec_deg": np.array([20.0, 20.0]),
    }
    counts = comove.map_image(path, "qpb", times, weights, header, out_wcs, shape, eph)
    exposure = comove.map_image(
        path, "exposure", times, weights, header, out_wcs, shape, eph
    )
    assert np.isclose(counts.sum(), 100, atol=1e-8)
    assert np.isclose(exposure.sum(), 25, atol=1e-8)
    written = WORK / "comove-map.fits"
    comove.write_image(written, counts, out_wcs, "mos1", "soft", "qpb")
    with fits.open(written) as hdul:
        saved = hdul[0].header
        assert saved["COMOVING"] and saved["BUNIT"] == "count"
        assert (saved["DETECTOR"], saved["BAND"], saved["E_MIN"], saved["E_MAX"]) == (
            "mos1",
            "soft",
            *common.CFG["bands"]["soft"]["pi"],
        )
        assert np.isclose(np.asarray(hdul[0].data, float).sum(), 100, atol=1e-5)


def test_mask_equality():
    wcs = image_wcs()
    soft, hard = WORK / "soft.fits", WORK / "hard.fits"
    fits.PrimaryHDU(np.ones((101, 101)), wcs.to_header()).writeto(soft)
    fits.PrimaryHDU(np.full((101, 101), 2.0), wcs.to_header()).writeto(hard)
    ra0, dec0 = wcs.pixel_to_world_values(50, 50)
    ra1, dec1 = wcs.pixel_to_world_values(70, 35)
    catalog = [
        {"ra": str(ra0), "dec": str(dec0), "radius_arcsec": "12"},
        {"ra": str(ra1), "dec": str(dec1), "radius_arcsec": "8"},
    ]
    background.mask_sky_map(soft, catalog)
    background.mask_sky_map(hard, catalog)
    soft_data, hard_data = fits.getdata(soft), fits.getdata(hard)
    mask = soft_data == 0
    assert np.array_equal(mask, hard_data == 0)
    assert mask.any() and (~mask).any()
    assert np.all(soft_data[~mask] == 1) and np.all(hard_data[~mask] == 2)


def test_background_stamp():
    region = WORK / "region.txt"
    region.write_text("&&!((DETX,DETY) IN circle(1,2,3))\n")
    frame = {"soft": "10", "hard": "20"}
    event = {"region": str(region)}
    stamp = background.inputs_stamp(frame, event, "soft", "T F T", 350, 1100, 12)
    assert stamp == (
        "soft\n350-1100\n12\nT F T\n10\n20\n" "&&!((DETX,DETY) IN circle(1,2,3))\n"
    )
    region.write_text("&&!((DETX,DETY) IN circle(1,2,4))\n")
    assert (
        background.inputs_stamp(frame, event, "soft", "T F T", 350, 1100, 12) != stamp
    )
    assert (
        background.inputs_stamp(frame, event, "hard", "T F T", 2500, 8500, 12) != stamp
    )
    assert (
        background.inputs_stamp(frame, event, "soft", "F F T", 350, 1100, 12) != stamp
    )

    oot = WORK / "frame-oot.fits"
    output = WORK / "oot-detector.fits"
    output.touch()
    event["oot"] = str(oot)
    assert (
        background.inputs_stamp(frame, event, "soft", "T F T", 350, 1100, 12) != stamp
    )
    calls = []
    old_sas = background.sas
    background.sas = lambda tool, cwd, **kw: calls.append((tool, cwd, kw))
    try:
        background.make_oot_detector(WORK, oot, region, output, 0, 350, 1100)
    finally:
        background.sas = old_sas
    tool, cwd, parameters = calls[0]
    assert tool == "evselect" and cwd == WORK
    assert parameters["table"] == "frame-oot.fits:EVENTS"
    assert "PATTERN<=0" in parameters["expression"]
    assert "PI in [350:1100]" in parameters["expression"]
    assert "CCDNR == 1" in parameters["expression"]
    assert "CCDNR == 12" in parameters["expression"]
    assert "BOX(-2196,-1110,16060,15510,0)" in parameters["expression"]
    assert "circle(-2200,-1110,17980)" in parameters["expression"]
    assert "circle(1,2,4)" in parameters["expression"]
    assert parameters["imageset"] == output.name and not output.exists()
    assert parameters["xcolumn"] == "DETX" and parameters["ycolumn"] == "DETY"
    assert parameters["ximagesize"] == parameters["yimagesize"] == 780
    assert parameters["ximagemin"] == parameters["yimagemin"] == -19499
    assert parameters["ximagemax"] == parameters["yimagemax"] == 19500
    assert parameters["imagedatatype"] == "Int32"


def test_reduction_figure():
    path = WORK / "figure-events.fits"
    write_events(
        path,
        [0, 20, 40, 60],
        [0, 20, 40, 60],
        pi=[500, 500, 3000, 500],
        pattern=[0, 4, 0, 0],
        flag=[0, 0, 0, 1],
    )
    wcs, shape = figure.event_bounds([path], (100, 20))
    soft = figure.event_image([path], "soft", wcs, shape, True)
    hard = figure.event_image([path], "hard", wcs, shape, True)
    assert soft.sum() == 1 and hard.sum() == 1
    raw = figure.event_image([path], "soft", wcs, shape, False)
    assert raw.sum() == 3
    exposure = WORK / "figure-exposure.fits"
    fits.PrimaryHDU(np.ones(shape), wcs.to_header()).writeto(exposure)
    sampled = figure.exposure_image([exposure], wcs, shape)
    assert np.isclose(sampled.max(), 1) and sampled.sum() > 0

    edge = WORK / "figure-edge.fits"
    x = np.zeros(40001)
    x[39999] = 800
    write_events(edge, x, np.zeros_like(x))
    edge_wcs, edge_shape = figure.event_bounds([edge], (100, 20))
    assert figure.event_image([edge], "soft", edge_wcs, edge_shape, False).sum() == len(
        x
    )


def test_prune():
    root = WORK / "prune" / "outputs"
    raises(RuntimeError, prune.prune, REPRO.parent / "outside-analysis-repro")
    for stage in prune.GATES:
        (root / stage).mkdir(parents=True, exist_ok=True)
        (root / stage / ".done").write_text("complete\n")
    for relative in prune.MANIFESTS:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("value\n")

    retained = root / "02_reprocess" / "events" / "retained.fits"
    retained.parent.mkdir(parents=True)
    retained.touch()
    common.write_tsv(
        root / "02_reprocess" / "exposures.tsv", ("path",), [{"path": retained}]
    )
    for relative in ("03_events/frames.tsv", "04_sources/events.tsv"):
        common.write_tsv(root / relative, ("event",), [{"event": retained}])

    work = root / "05_background" / "pn-test" / "soft"
    work.mkdir(parents=True)
    names = {
        "source_pha": "pn-fovt.pi",
        "qpb_pha": "pn-bkg.pi",
        "rmf": "pn.rmf",
        "arf": "pn.arf",
        "events_image": "pn-fovimsky-350-1100.fits",
        "exposure_map": "pn-expimsky-350-1100.fits",
        "qpb_map": "pn-qpb-sky-350-1100.fits",
    }
    paths = {key: work / name for key, name in names.items()}
    paths.update(
        fovimdet=work / "pn-fovimdet-350-1100.fits",
        fovtootsub=work / "pn-fovtootsub.pi",
    )
    for path in paths.values():
        path.touch()
    background = {
        "frame": "pn-test",
        "det": "pn",
        "band": "soft",
        "dir": work,
        **{key: paths[key] for key in names},
    }
    common.write_tsv(
        root / "05_background" / "background.tsv", tuple(background), [background]
    )

    sp_pha = root / "06_soft_protons" / "pn-test-sp.pi"
    sp_map = root / "06_soft_protons" / "pn-test-sp.fits"
    sp_pha.touch()
    sp_map.touch()
    common.write_tsv(
        root / "06_soft_protons" / "soft_protons.tsv",
        ("sp_pha", "sp_map"),
        [{"sp_pha": sp_pha, "sp_map": sp_map}],
    )
    summary = root / "01_init" / "odf" / "testSUM.SAS"
    summary.parent.mkdir(parents=True)
    summary.touch()
    odf_copy = summary.parent / "testODF.FIT"
    odf_copy.touch()
    (root / "01_init" / "ccf.cif").touch()
    (root / "01_init" / ".done").touch()
    attitude = root / "02_reprocess" / "atthk.fits"
    attitude.touch()
    resume = [
        root / "02_reprocess" / "emproc" / "test_ImagingEvts.ds",
        root / "02_reprocess" / "epproc" / "test_ImagingEvts.ds",
        root / "02_reprocess" / "epproc_oot" / "test_ImagingEvts.ds",
        root / "02_reprocess" / "filter" / "test" / "test-gti.fits",
        root / "02_reprocess" / "filter" / "test" / "test-allevc.fits",
        root / "02_reprocess" / "filter" / "test" / "test-allevcoot.fits",
        root / "02_reprocess" / "filter" / "test" / "test-corevc.fits",
    ]
    for path in resume:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    regions = root / "04_sources" / "sources.reg"
    regions.touch()

    final = root / "07_comove" / "stack" / "final.fits"
    final.parent.mkdir(parents=True)
    final.touch()
    figure_product = root / "08_figure" / "figure.pdf"
    figure_product.touch()
    scratch = work / "unreferenced.fits"
    scratch.touch()
    external = WORK / "external.fits"
    external.touch()
    link = work / "external-link.fits"
    link.symlink_to(external)
    external_dir = WORK / "external-directory"
    external_dir.mkdir()
    external_child = external_dir / "child.fits"
    external_child.touch()
    directory_link = work / "external-directory"
    directory_link.symlink_to(external_dir, target_is_directory=True)
    pfiles = root / "07_comove" / ".pfiles"
    pfiles.mkdir()
    (pfiles / "task.par").touch()
    bootstrap = root / "tmp" / "cache"
    bootstrap.mkdir(parents=True)
    (bootstrap / "fontlist.json").touch()
    logs = root / "logs"
    logs.mkdir()
    success = logs / "success.log"
    success.write_text("tool started\ntool ended:    2026-08-30T01:00:00.000\n")
    failure = logs / "failure.log"
    failure.write_text("** tool: error (broken), stopped\n")

    outside = WORK / "outside.tmp"
    outside.touch()
    raises(RuntimeError, prune.remove, outside, root)
    raises(RuntimeError, prune.remove, directory_link / "child.fits", root)
    removed = prune.prune(root)
    assert removed >= 4
    assert all(
        path.exists()
        for path in (
            *paths.values(),
            retained,
            sp_pha,
            sp_map,
            summary,
            odf_copy,
            attitude,
            final,
            figure_product,
            regions,
            *resume,
        )
    )
    assert not scratch.exists() and not link.exists() and external.exists()
    assert not directory_link.exists() and external_child.exists()
    assert (
        not success.exists()
        and failure.exists()
        and not pfiles.exists()
        and not (root / "tmp").exists()
    )

    resume[0].unlink()
    sentinel = root / "02_reprocess" / "sentinel.tmp"
    sentinel.touch()
    raises(FileNotFoundError, prune.prune, root)
    assert sentinel.exists()

    missing_gate = WORK / "missing-gate" / "outputs"
    missing = missing_gate / "01_init" / "scratch"
    missing.parent.mkdir(parents=True)
    missing.touch()
    raises(RuntimeError, prune.prune, missing_gate)
    assert missing.exists()


def main():
    clear_work()
    WORK.mkdir(parents=True)
    try:
        test_output_boundary()
        test_gti()
        test_event_support_and_regions()
        test_pha()
        test_hard_flare_correction()
        test_comove_maps()
        test_mask_equality()
        test_background_stamp()
        test_reduction_figure()
        test_prune()
    finally:
        clear_work()


if __name__ == "__main__":
    main()
