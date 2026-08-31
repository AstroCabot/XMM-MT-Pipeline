import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time

from analysis_02_pps_luminosity import run, selection

assert run.SETTINGS["python_lib"] == "/home/scabot/miniconda3/envs/xmm/lib"
assert run.FINAL_PRODUCTS == (
    "result.json",
    "response.arf",
    "source-350-500.pi",
    "background-350-500.pi",
    "source-500-1000.pi",
    "background-500-1000.pi",
    "source-1000-1100.pi",
    "background-1000-1100.pi",
)
assert run.SOFT_BAND == tuple(run.PARAMS["soft_band_ev"])
assert (run.BANDS[0][0], run.BANDS[-1][1]) == run.SOFT_BAND
energy = np.array([349, 350, 499, 500, 999, 1000, 1100, 1101])
subbands = [
    selection.band_mask(energy, lo, hi, include_high=i == len(run.BANDS) - 1)
    for i, (lo, hi) in enumerate(run.BANDS)
]
assert np.array_equal(np.sum(subbands, axis=0), [0, 1, 1, 1, 1, 1, 1, 0])
assert not subbands[0][3] and subbands[1][3]
assert subbands[-1][-2]
assert selection.gti_expression().count("CCDNR==") == 12
times, weights = selection.time_samples([np.array([[0.0, 70.0], [100.0, 130.0]])], 60)
assert np.isclose(weights.sum(), 100)
assert np.all((times >= 0) & (times <= 130))
header = fits.Header({"MJDREF": 50814.0, "TIMESYS": "TT"})
assert np.all(np.diff(selection.event_mjd([0, 1], header)) > 0)

ephemeris = np.array([[1.0, 10.0, 2.0], [2.0, 11.0, 3.0]])
ra, dec = selection.comet_position(np.array([1.25, 1.75]), ephemeris)
assert np.allclose(ra, [10.25, 10.75])
assert np.allclose(dec, [2.25, 2.75])

event_header = fits.Header({"MJDREF": 60000.0, "TIMESYS": "UTC"})
sample_mjd = selection.event_mjd([0, 30, 90, 120], event_header)
shift = 20 / 3600 / np.cos(np.deg2rad(2))
moving_ephemeris = np.column_stack(
    (sample_mjd, [10, 10, 10 + shift, 10 + shift], np.full(4, 2.0))
)
coverage = selection.exposure_coverage(
    [np.array([[0.0, 120.0]])],
    event_header,
    [(10.0, 2.0)],
    moving_ephemeris,
    (600, 620),
    radius_arcsec=4,
)
assert np.isclose(coverage.min(), 0.5)
assert np.any((coverage > 0) & (coverage < 1)) and coverage[0, 0] == 1

source = (10.0, 2.0)
east, north = selection.source_offsets(source, np.array([10.0]), np.array([2.0]))
assert np.allclose([east[0], north[0]], 0)
east, north = selection.source_offsets(
    (177.882617, 1.381687), np.array([178.105677774]), np.array([1.308813404])
)
assert np.allclose([east[0], north[0]], [-802.7974583, 262.3449456])
keep = selection.source_keep(
    [selection.EVENT_CENTER[0], selection.EVENT_CENTER[0] + 1000],
    [selection.EVENT_CENTER[1], selection.EVENT_CENTER[1]],
    np.array([1.0, 1.0]),
    [source],
    ephemeris,
    run.PARAMS["source_mask_radius_arcsec"],
)
assert np.array_equal(keep, [False, True])

yy, xx = np.indices((400, 400), dtype=float)
radius = np.hypot(xx - 199.5, yy - 199.5)
source_region = radius <= 100
background_region = (radius >= 150) & (radius < 180)
uniform = np.ones(radius.shape)
measured = selection.correction(uniform, radius, source_region, background_region)
expected = 1 / (1 - 100 / (150 + 180))
assert abs(measured / expected - 1) < 3e-3
source_area, background_area = 90.0, 180.0
source_time, background_time = 30.0, 20.0
alpha = (source_area / source_time) / (background_area / background_time)
assert np.isclose(
    alpha / background_time, (source_area / background_area) / source_time
)

with tempfile.TemporaryDirectory(dir=run.HERE) as name:
    path = Path(name) / "flux.fits"
    table = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="ENERGY", format="E", array=np.array([0.4, 0.5, 0.6])),
            fits.Column(name="ENERGY_BIN", format="E", array=np.full(3, 0.1)),
            fits.Column(name="FLUX", format="E", array=np.full(3, 2.0)),
        ]
    )
    fits.HDUList([fits.PrimaryHDU(), table]).writeto(path)
    assert np.isclose(run.integrate_flux(path, 0.4, 0.6), 0.4)

    pha = Path(name) / "spectrum.pi"
    spectrum = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="CHANNEL", format="J", array=np.arange(4096)),
            fits.Column(name="COUNTS", format="J", array=np.zeros(4096, int)),
        ],
        name="SPECTRUM",
    )
    fits.HDUList([fits.PrimaryHDU(), spectrum]).writeto(pha)
    counts = np.zeros(4096, int)
    counts[:3] = [2, 3, 4]
    run._set_pha(pha, counts, 1.5)
    with fits.open(pha) as hdul:
        assert np.array_equal(hdul["SPECTRUM"].data["COUNTS"], counts)
        assert hdul["SPECTRUM"].header["BACKSCAL"] == 1.5
        assert hdul["SPECTRUM"].header["BACKFILE"] == "none"
    assert len(fits.getdata(run.RMF, "EBOUNDS")) == 4096

with tempfile.TemporaryDirectory(dir=run.HERE) as name:
    folder = Path(name)
    event = folder / "event.fits"
    primary = fits.PrimaryHDU()
    primary.header["MTFLAG"] = True
    events = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="TIME", format="D", array=np.array([0.0])),
        ],
        name="EVENTS",
    )
    fits.HDUList([primary, events]).writeto(event)
    sources_path, ephemeris_path = folder / "sources.tsv", folder / "ephemeris.tsv"
    sources_path.write_text(
        "ra\tdec\n" + "".join(f"{index}\t{index + 1}\n" for index in range(15))
    )
    ephemeris_path.write_text("mjd\tra_deg\tdec_deg\n1\t2\t3\n2\t3\t4\n")
    sources, ephemeris, _ = run.motion_inputs(event, sources_path, ephemeris_path)
    assert len(sources) == 15 and ephemeris.shape == (2, 3)
    sources_path.write_text("ra\tdec\n1\t2\n")
    try:
        run.motion_inputs(event, sources_path, ephemeris_path)
    except ValueError:
        pass
    else:
        raise AssertionError("truncated PPS source table accepted")

rows = [
    {
        "raw_masked_flux_erg_cm-2_s-1": 2.0,
        "annulus_self_subtraction_factor_1_over_r": 1.2,
        "source_mask_fill_factor_1_over_r": 1.1,
    },
    {
        "raw_masked_flux_erg_cm-2_s-1": 3.0,
        "annulus_self_subtraction_factor_1_over_r": 1.4,
        "source_mask_fill_factor_1_over_r": 1.2,
    },
]
summary = run.summarize_flux(rows)
assert set(summary) == {"flux_luminosity", "correction_factors"}
products = summary["flux_luminosity"]
assert set(products) == {
    "raw_masked",
    "annulus_self_subtraction_corrected_masked",
    "full_aperture",
}
assert np.isclose(products["raw_masked"]["flux_erg_cm-2_s-1"], 5.0)
assert np.isclose(
    products["annulus_self_subtraction_corrected_masked"]["flux_erg_cm-2_s-1"], 6.6
)
assert np.isclose(products["full_aperture"]["flux_erg_cm-2_s-1"], 7.68)
factors = summary["correction_factors"]
assert set(factors) == {
    "annulus_self_subtraction_exposure_weighted_1_over_r",
    "source_mask_fill_1_over_r",
    "raw_masked_to_full_aperture_1_over_r",
}
assert np.isclose(
    factors["annulus_self_subtraction_exposure_weighted_1_over_r"], 6.6 / 5.0
)
assert np.isclose(1 / (1 - run.APERTURE / sum(run.BACKGROUND)), 1.3894736842)
assert np.isclose(factors["source_mask_fill_1_over_r"], 7.68 / 6.6)
assert np.isclose(factors["raw_masked_to_full_aperture_1_over_r"], 7.68 / 5.0)
assert all(
    set(product) == {"flux_erg_cm-2_s-1", "luminosity_erg_s"}
    for product in products.values()
)

bands = [
    {
        "energy_eV": [350, 500],
        "source_counts": 10,
        "background_counts": 2,
        "background_scale": 0.1,
        "source_exposure_s_arcmin2": 5,
        "background_exposure_s_arcmin2": 8,
        "raw_masked_flux_erg_cm-2_s-1": 2,
        "annulus_self_subtraction_factor_1_over_r": 1.2,
        "source_mask_fill_factor_1_over_r": 1.1,
    }
]
result = run.make_result(Time("2025-12-03T21:25:46.398", scale="utc"), bands)
assert set(result) == {
    "epoch_utc",
    "distance_au",
    "selection",
    "response",
    "bands",
    "flux_luminosity",
    "correction_factors",
    "units",
}
assert set(result["selection"]) == {
    "energy_eV",
    "aperture_arcsec",
    "background_arcsec",
    "event_filter",
    "source_mask_radius_arcsec",
}
assert set(result["response"]) == {"arf", "arf_model", "rmf"}
assert set(result["bands"][0]) == set(bands[0])
assert set(result["units"]) == {"counts", "background_scale", "flux", "luminosity"}
