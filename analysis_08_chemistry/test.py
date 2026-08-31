import csv
import importlib.util
import json
import math
import os
import sys
import tempfile
from itertools import pairwise
from pathlib import Path
from types import SimpleNamespace

import chemistry
import diagnostic
import figure
import inference
import numpy as np
import products
import psf
import psf_operator
import run
from astropy.io import fits
from astropy.time import Time
from config import CONFIG
from scipy.integrate import trapezoid
from transport import ARCMIN2_SR, EV_ERG, Orbit, Transport, tau_one_radius

HERE = Path(__file__).resolve().parent
LOCAL = json.loads((HERE / "parameters.json").read_text())
assert "production_reference_r_au" not in LOCAL
assert "H2" not in LOCAL["cross_section_cm2"]
assert not (
    {
        "Qref_H2O",
        "Qref_CO2",
        "Qref_CO",
        "beta_H2O",
        "beta_CO2",
        "beta_CO",
        "heavy_ion_flux",
    }
    & set(LOCAL["priors"])
)
assert (
    CONFIG["production_reference_r_au"]
    == run.SHARED["spherex"]["reference_distance_au"]
)
assert (
    CONFIG["cross_section_cm2"]["H2"]
    == run.SHARED["h2_single_capture_cross_section_cm2"]
)
EXPECTED_PHYSICAL = (
    "Q_H2O",
    "Q_CO2",
    "Q_CO",
    "Q_H2",
    "beta_H2O",
    "beta_CO2",
    "beta_CO",
    "beta_H2",
    "log10_H2_ratio",
    "phi_heavy",
    "v_H2O_km_s",
    "v_CO2_km_s",
    "v_CO_km_s",
    "v_H2_km_s",
    "v_OH_km_s",
    "v_O_H2O_km_s",
    "v_O_OH_km_s",
    "v_H_km_s",
    "scale",
    "offset",
    "jitter",
)
assert inference.PHYSICAL_NAMES == EXPECTED_PHYSICAL
probe = np.ones((1, len(inference.PHYSICAL_NAMES)))
probe[0, :4] = (1e28, 1e27, 1e26, 1e29)
for name, value in (
    ("beta_H2", 0.4),
    ("phi_heavy", 1e5),
    ("v_H2_km_s", 0.37),
    ("scale", 1.1),
    ("offset", 2e-5),
    ("jitter", 0.08),
):
    probe[0, inference.PHYSICAL_NAMES.index(name)] = value
assert np.allclose(
    diagnostic.values(probe)[0], (28, 27, 26, 29, 0.4, 5, 0.37, 1.1, 2, 0.08)
)
assert products.orbit_output() == {
    "perihelion_q_au": CONFIG["orbit"]["perihelion_q_au"],
    "eccentricity": CONFIG["orbit"]["eccentricity"],
    "perihelion_jd_tdb": CONFIG["orbit"]["perihelion_jd_tdb"],
}
assert products.FINAL_PRODUCTS == (
    "result.json",
    "chain.npz",
    "posterior.tsv",
    "profile.tsv",
    "intrinsic.tsv",
    "chemistry.tsv",
    "figure3.png",
    "figure3.pdf",
    "corner.png",
)
assert psf.frame_detector("pn-s003-p11") == "pn"
assert (
    CONFIG["mcmc"]["walkers"],
    CONFIG["mcmc"]["burn_steps"],
    CONFIG["mcmc"]["retained_steps"],
) == (96, 500, 1500)


def test_psf_environment(work):
    captured = {}

    def invoke(*args, **kwargs):
        captured.update(kwargs["env"])
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    old = psf.subprocess.run
    psf.subprocess.run = invoke
    expected = {"SAS_CCF": "ccf", "SAS_CCFPATH": "path", "SAS_ODF": "odf"}
    try:
        psf.run_psfgen(
            work / "psf.fits", work / "exp.fits", "pn", 1, 2, 500, work, expected
        )
    finally:
        psf.subprocess.run = old
    assert all(captured[name] == value for name, value in expected.items())


with tempfile.TemporaryDirectory(dir=HERE) as temporary:
    test_psf_environment(Path(temporary))


ellbeta = np.zeros((8, 8))
ellbeta[4, 4] = 1
ellbeta_header = fits.Header(
    {
        "CRPIX1": 4.0,
        "CRPIX2": 4.0,
        "CDELT1": -1 / 3600,
        "CDELT2": 1 / 3600,
        "CUNIT1": "deg",
        "CUNIT2": "deg",
    }
)
radius, density = psf.radial_density(ellbeta, ellbeta_header)
area = np.pi * ((radius + 0.5) ** 2 - (radius - 0.5) ** 2)
assert np.isclose(np.sum(density * area), 1)
asymmetric = ellbeta.copy()
asymmetric[7, 4] = 0.1
psf.radial_density(asymmetric, ellbeta_header)
for image, header in (
    (np.roll(ellbeta, 1, axis=0), ellbeta_header),
    (ellbeta, fits.Header({**ellbeta_header, "CDELT1": -2 / 3600})),
    (ellbeta, fits.Header({**ellbeta_header, "CRPIX1": 5})),
    (ellbeta, fits.Header({**ellbeta_header, "CUNIT1": "arcsec"})),
):
    try:
        psf.radial_density(image, header)
    except RuntimeError:
        pass
    else:
        raise AssertionError("invalid ELLBETA grid accepted")


for header in (
    fits.Header({"MJDREF": 50814.0, "TIMESYS": "TT"}),
    fits.Header({"MJDREFI": 50814, "MJDREFF": 0.0, "TIMESYS": "TT"}),
):
    assert math.isclose(
        psf.event_mjd(0, header), Time(50814, format="mjd", scale="tt").utc.mjd
    )

with tempfile.TemporaryDirectory(dir=HERE) as folder:
    event = Path(folder) / "event.fits"
    events = fits.BinTableHDU.from_columns([], name="EVENTS")
    events.header.update(MJDREF=50814.0, TIMESYS="TT", ONTIME=30.0)
    gti = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="START", format="D", array=[0.0, 80.0]),
            fits.Column(name="STOP", format="D", array=[10.0, 100.0]),
        ],
        name="STDGTI",
    )
    fits.HDUList([fits.PrimaryHDU(), events, gti]).writeto(event)
    endpoint = psf.event_mjd(np.array([0.0, 100.0]), events.header)
    ephemeris = np.array([[endpoint[0], 10.0, 20.0], [endpoint[1], 11.0, 22.0]])
    position = psf.mean_position(
        *psf.comet_track(
            {"frame": "test", "start": "0", "stop": "100"}, event, ephemeris
        )
    )
    weighted_time = (5 * 10 + 90 * 20) / 30
    assert np.allclose(
        position, [10 + weighted_time / 100, 20 + 2 * weighted_time / 100], atol=1e-10
    )

probe = np.zeros((3, len(inference.COORDS)))
probe[:, inference.COORDS.index("log10_H2_ratio")] = (0, 0.2, 0.5)
high, high_best, low_best = products.mode_selection(probe, np.array([1.0, 3.0, 2.0]))
assert high.all() and high_best == 1 and low_best is None
probe[:, inference.COORDS.index("log10_H2_ratio")] = -2
try:
    products.mode_selection(probe, np.arange(3.0))
except RuntimeError as error:
    assert "high-H2" in str(error)
else:
    raise AssertionError("missing high-H2 mode accepted")


def close(left, right, relative=1e-10, absolute=0):
    assert math.isclose(left, right, rel_tol=relative, abs_tol=absolute), (left, right)


def angular_separation(left, right):
    ra1, dec1 = np.deg2rad(left)
    ra2, dec2 = np.deg2rad(right)
    dra = (ra1 - ra2 + np.pi) % (2 * np.pi) - np.pi
    return np.rad2deg(np.hypot(dra * np.cos(0.5 * (dec1 + dec2)), dec1 - dec2)) * 3600


task03_path = HERE.parent / "analysis_03_spectrum"
sys.path.insert(0, str(task03_path))
try:
    specification = importlib.util.spec_from_file_location(
        "task03_extract_contract", task03_path / "extract.py"
    )
    task03_extract = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(task03_extract)
finally:
    sys.path.pop(0)
real_frames = psf.read_tsv(run.REDUCTION / "03_events/frames.tsv")
real_events = {
    row["frame"]: row for row in psf.read_tsv(run.REDUCTION / "04_sources/events.tsv")
}
real_ephemeris_rows = psf.read_tsv(run.REDUCTION / "07_comove/ephemeris.tsv")
ephemeris_names = ("mjd", "ra_deg", "dec_deg", "r_au", "delta_au")
real_ephemeris = {
    name: np.array([float(row[name]) for row in real_ephemeris_rows])
    for name in ephemeris_names
}
real_ephemeris_array = np.column_stack(
    [real_ephemeris[name] for name in ephemeris_names[:3]]
)
old_discrepancy = []
for frame in real_frames:
    event = Path(real_events[frame["frame"]]["event"])
    actual = psf.mean_position(*psf.comet_track(frame, event, real_ephemeris_array))
    expected = task03_extract.center(
        {"frame": frame["frame"], "event": event}, real_ephemeris
    )
    expected_position = (expected["ra_deg"], expected["dec_deg"])
    assert angular_separation(actual, expected_position) < 1e-7
    header = fits.getheader(event, "EVENTS")
    midpoint = psf.event_mjd(
        0.5 * (float(frame["start"]) + float(frame["stop"])), header
    )
    old = (
        float(
            np.rad2deg(
                np.interp(
                    midpoint,
                    real_ephemeris["mjd"],
                    np.unwrap(np.deg2rad(real_ephemeris["ra_deg"])),
                )
            )
            % 360
        ),
        float(np.interp(midpoint, real_ephemeris["mjd"], real_ephemeris["dec_deg"])),
    )
    old_discrepancy.append(angular_separation(old, expected_position))
assert len(old_discrepancy) == 22
assert 26 < np.median(old_discrepancy) < 27
assert 80 < max(old_discrepancy) < 81

# Frame-exposure reconstruction.
real_settings = json.loads((run.REDUCTION.parent / "settings.json").read_text())
real_maps = psf.exposure_maps(
    run.REDUCTION, set(real_settings["expected_frames"]), run.SHARED["soft_band_ev"]
)
real_edges = np.asarray(run.SHARED["chemistry_annulus_edges_arcsec"], float)
real_frame_weights, real_centers, detector_closure = psf.frame_exposure_weights(
    run.REDUCTION,
    real_frames,
    real_events,
    real_maps,
    real_ephemeris_array,
    run.SHARED["soft_band_ev"],
    real_edges,
    float(real_settings["comove"]["map_step_s"]),
)
assert real_frame_weights.shape == (22, 10) and real_centers.shape == (22, 2)
assert detector_closure.shape == (3, 10)
assert np.max(np.abs(detector_closure)) < 1e-7
trial_scales = {"pn": 0.83, "mos1": 1.17, "mos2": 1.0}
real_stacks = psf.detector_exposure_stacks(run.REDUCTION, run.SHARED["soft_band_ev"])
combined_plane = sum(trial_scales[name] * item[1] for name, item in real_stacks.items())
combined_header = next(iter(real_stacks.values()))[2]
combined_expected = psf_operator.annulus_exposure(
    combined_plane, combined_header, real_edges
)
combined_actual = np.sum(
    real_frame_weights
    * np.asarray([trial_scales[frame["det"]] for frame in real_frames])[:, None],
    axis=0,
)
combined_closure = psf_operator.fractional_closure(
    combined_actual, combined_expected, "synthetic Task 04"
)
assert np.max(np.abs(combined_closure)) < 1e-7


rows = chemistry.reactions(HERE / "reactions.tsv", 175.0)
network = chemistry.network(rows, CONFIG["cross_section_cm2"])
assert len(rows) == 35
expected_reactions = """H|H+ + e|1.20e-8|0.45|0|0|0|0|0|0|0|0|0
H2|H + H|1.03e-8|0.36|0|0|0|2|0|0|0|0|0
H2|H + H|6.84e-9|0.50|0|0|0|2|0|0|0|0|0
H2|H2+ + e|1.36e-8|0.45|0|0|0|0|0|0|0|0|0
H2|H + H+ + e|2.16e-10|0.89|0|0|0|1|0|0|0|0|0
H2O|OH + H|3.60e-6|0.24|0|1|0|1|0|0|0|0|0
H2O|H2 + O|1.80e-7|0.36|0|0|1|0|0|0|0|1|0
H2O|O + H + H|2.40e-7|0.36|0|0|0|2|0|0|0|1|0
H2O|OH+ + H + e|2.87e-9|0.75|0|0|0|1|0|0|0|0|0
H2O|O+ + H2 + e|4.00e-11|1.15|0|0|1|0|0|0|0|0|0
H2O|H+ + OH + e|3.32e-10|0.89|0|1|0|0|0|0|0|0|0
H2O|H2O+ + e|4.23e-8|0.56|0|0|0|0|0|0|0|0|0
OH|O(3P) + H|3.97e-6|0.059|0|0|0|1|0|0|0|1|0
OH|O(1S) + H|1.74e-7|0.36|0|0|0|1|0|0|0|1|0
OH|OH+ + e|1.53e-8|0.68|0|0|0|0|0|0|0|0|0
OH|O(1D) + H|1.38e-6|0.35|0|0|0|1|0|0|0|1|0
CO|CO+ + e|3.83e-8|0.61|0|0|0|0|0|0|0|0|0
CO|C+ + O + e|3.48e-10|1.03|0|0|0|0|0|0|0|1|0
CO|O+ + C + e|2.10e-10|1.09|0|0|0|0|0|0|0|0|1
CO|C + O|4.77e-8|0.49|0|0|0|0|0|0|0|1|1
CO|C(1D) + O(1D)|1.78e-8|0.31|0|0|0|0|0|0|0|1|1
CO2|CO + O(1D)|2.22e-7|0.42|0|0|0|0|0|1|0|1|0
CO2|CO2+ + e|5.09e-8|0.66|0|0|0|0|0|0|0|0|0
CO2|CO+ + O + e|7.08e-10|1.00|0|0|0|0|0|0|0|1|0
CO2|O+ + CO + e|8.57e-10|1.01|0|0|0|0|0|1|0|0|0
CO2|C+ + O2 + e|1.98e-10|1.15|0|0|0|0|0|0|1|0|0
CO2|CO + O(3P)|8.44e-9|0.049|0|0|0|0|0|1|0|1|0
CO2|CO(a3P) + O(3P)|6.37e-8|0.47|0|0|0|0|0|1|0|1|0
O|O+ + e|1.81e-8|0.67|0|0|0|0|0|0|0|0|0
O2|O + O|3.23e-8|0.24|0|0|0|0|0|0|0|2|0
O2|O + O(1D)|1.44e-6|0.16|0|0|0|0|0|0|0|2|0
O2|O+ + O + e|2.25e-9|0.92|0|0|0|0|0|0|0|1|0
O2|O(1S) + O(1S)|7.29e-9|0.47|0|0|0|0|0|0|0|2|0
O2|O2+ + e|4.03e-8|0.63|0|0|0|0|0|0|0|0|0
C|C+ + e|9.51e-8|0.50|0|0|0|0|0|0|0|0|0""".splitlines()
actual_reactions = [
    line.replace("\t", "|")
    for line in (HERE / "reactions.tsv").read_text().splitlines()[1:]
]
assert actual_reactions == expected_reactions
expected_lifetime = {
    "H2O": 62190,
    "OH": 64820,
    "CO2": 218800,
    "CO": 584300,
    "H2": 3155000,
    "O": 1736000,
    "C": 794900,
    "H": 8156000,
    "O2": 208800,
}
for name, expected in expected_lifetime.items():
    close(1 / network[name]["rate_1au_s-1"], expected, relative=8e-4)
yield_checks = {
    ("H2O", "OH"): 0.77535,
    ("H2O", "H2"): 0.072804,
    ("H2O", "H"): 0.97352,
    ("H2O", "O"): 0.16767,
    ("CO2", "CO"): 0.61985,
    ("CO", "C"): 0.43592,
    ("H2", "H"): 1.05565,
    ("O2", "O"): 1.50983,
}
for (parent, child), expected in yield_checks.items():
    close(dict(network[parent]["children"])[child], expected, relative=2e-5)
assert [
    len(chemistry.chains(name, network)) for name in ("H2O", "CO2", "CO", "H2")
] == [8, 7, 3, 2]
assert all(
    len(chain[0]) <= 3
    for name in ("H2O", "CO2", "CO", "H2")
    for chain in chemistry.chains(name, network)
)

orbit = Orbit(
    CONFIG["orbit"]["perihelion_q_au"],
    CONFIG["orbit"]["eccentricity"],
    CONFIG["orbit"]["perihelion_jd_tdb"],
    CONFIG["orbit"]["observation_jd_tdb"],
    1e10,
)
close(orbit.exposure(2e10), orbit.exposure(1e10))
assert orbit.inverse_r2(2e10) == 0
close(orbit.current_r_au, 1.8688760328902, relative=2e-12)
close(
    orbit.radius(orbit.perihelion_age_s),
    CONFIG["orbit"]["perihelion_q_au"],
    relative=2e-12,
)
cap = orbit.activity(
    3.3, np.array([orbit.perihelion_age_s, 2 * orbit.perihelion_age_s])
)
close(cap[0], cap[1], relative=2e-12)
age = np.linspace(0, orbit.perihelion_age_s, 2000)
direct_exposure = trapezoid(orbit.radius(age) ** -2, age)
close(orbit.exposure(age[-1]), direct_exposure, relative=2e-4)

source_radius = np.geomspace(0.01, 30, 40)
q = {"H2O": 1.8e28, "CO2": 3.2e27, "CO": 8.5e27, "H2": 7e28}
beta = {"H2O": 3.3, "CO2": 0.392, "CO": 1.3, "H2": 0.5}
speed = {name: 0.37 for name in q}
speed.update(OH=1.2, O_H2O=1.5, O_OH=1.5, H=12, heavy=1.685, H2_photo=11.335)
actual_chains = {name: chemistry.chains(name, network) for name in q}
fine = Transport(orbit, network, 160, 1.3, 3e9, source_radius, 81920.5, 1, 1)
fast = Transport(orbit, network, 160, 1.3, 3e9, source_radius, 81920.5, 4, 4)
fine_tau = fine.source_profile(q, beta, speed, actual_chains, 1e5, 502.2, 1)[1]
fast_tau = fast.source_profile(q, beta, speed, actual_chains, 1e5, 502.2, 1)[1]
assert np.allclose(fast_tau, fine_tau, rtol=0.01)


class FixedOrbit:
    current_r_au = 1.0

    @staticmethod
    def exposure(age):
        value = np.asarray(age, float)
        if np.any(value < 0):
            raise ValueError
        return value

    @staticmethod
    def inverse_r2(age):
        return np.ones_like(np.asarray(age, float))

    @staticmethod
    def activity(beta, age):
        return np.ones_like(np.asarray(age, float))


synthetic_network = {
    name: {"rate_1au_s-1": 1e-4, "sigma_cm2": 1e-15} for name in ("P", "D", "G")
}
transport = Transport(
    FixedOrbit(),
    synthetic_network,
    160,
    1,
    1e5,
    np.array([1.0]),
    1.0,
    granddaughter_stride=4,
)
speed = {"P": 1.0, "D": 1.0, "G": 1.0}
x = (transport.radius - 1) * 1e-4
chains = (
    ((("P", "P"),), 1.0),
    ((("P", "P"), ("D", "D")), 1.0),
    ((("P", "P"), ("D", "D"), ("G", "G")), 1.0),
)
flux = []
for chain in chains:
    density = transport.density(chain, 0, speed)
    flux.append(4 * np.pi * (transport.radius * 1e5) ** 2 * 1e5 * density)
assert np.allclose(flux[0], np.exp(-x), rtol=2e-12)
use = x < 6
assert np.allclose(flux[1][use], (x * np.exp(-x))[use], rtol=7e-3, atol=2e-8)
assert np.allclose(flux[2][use], (0.5 * x**2 * np.exp(-x))[use], rtol=7e-2, atol=2e-8)
synthetic_network["G"]["rate_1au_s-1"] = 0.0
stable = transport.density(chains[2], 0, speed)
stable_flux = 4 * np.pi * (transport.radius * 1e5) ** 2 * 1e5 * stable
assert np.allclose(stable_flux, 1 - np.exp(-x) * (1 + x), rtol=3e-2, atol=2e-8)


class VaryingOrbit(FixedOrbit):
    @staticmethod
    def exposure(age):
        age = np.asarray(age, float)
        return age + 0.025 * age**2

    @staticmethod
    def inverse_r2(age):
        return 1 + 0.05 * np.asarray(age, float)


varying_network = {
    "P": {"rate_1au_s-1": 0.08, "sigma_cm2": 1e-15},
    "D": {"rate_1au_s-1": 0.03, "sigma_cm2": 1e-15},
}
varying = Transport(VaryingOrbit(), varying_network, 160, 1, 100, np.array([1.0]), 1.0)
daughter = varying.density(((("P", "P"), ("D", "D")), 1), 0, {"P": 1, "D": 2})
flux = 4 * np.pi * (varying.radius * 1e5) ** 2 * 2e5 * daughter
for target in (8, 16, 32):
    index = np.argmin(abs(varying.radius - target))
    birth = np.linspace(1, varying.radius[index], 200001)
    final_age = (varying.radius[index] - birth) / 2
    emission_age = final_age + birth - 1
    integrand = (
        0.08
        * VaryingOrbit.inverse_r2(final_age)
        * np.exp(
            -0.08
            * (VaryingOrbit.exposure(emission_age) - VaryingOrbit.exposure(final_age))
            - 0.03 * VaryingOrbit.exposure(final_age)
        )
    )
    close(flux[index], np.trapezoid(integrand, birth), relative=2e-4)

serial_network = {
    name: {"rate_1au_s-1": 0.08, "sigma_cm2": 1e-15} for name in ("P", "D", "G")
}
serial = Transport(
    VaryingOrbit(), serial_network, 160, 1, 100, np.array([1.0]), 1.0, 1, 1
)
granddaughter = serial.density(
    ((("P", "P"), ("D", "D"), ("G", "G")), 1), 0, {"P": 1, "D": 1, "G": 1}
)
flux = 4 * np.pi * (serial.radius * 1e5) ** 2 * 1e5 * granddaughter
exposure = 0.08 * VaryingOrbit.exposure(serial.radius - 1)
expected = 0.5 * exposure**2 * np.exp(-exposure)
for target in (8, 16, 32):
    index = np.argmin(abs(serial.radius - target))
    close(flux[index], expected[index], relative=4e-4)

# Infinite-lifetime Haser column.
haser_network = {"P": {"rate_1au_s-1": 0.0, "sigma_cm2": 3e-15}}
source_radius = np.array([100.0])
haser = Transport(FixedOrbit(), haser_network, 400, 0.01, 1e8, source_radius, 1.0)
density = 1e28 * haser.density(((("P", "P"),), 1.0), 0, {"P": 1.0})
column = float(density @ haser.chord_cm[:, 0])
expected_column = 1e28 / (4 * 1e5 * source_radius[0] * 1e5)
close(column, expected_column, relative=6e-3)

# Thin and saturated charge exchange.
chains_one = {"P": [((("P", "P"),), 1.0)]}
thin, tau = haser.source_profile(
    {"P": 1e20}, {"P": 0}, {"P": 1}, chains_one, 1e5, 500, 1
)
twice, _ = haser.source_profile(
    {"P": 2e20}, {"P": 0}, {"P": 1}, chains_one, 1e5, 500, 1
)
close(twice[0] / thin[0], 2.0, relative=2e-7)
thick, thick_tau = haser.source_profile(
    {"P": 1e40}, {"P": 0}, {"P": 1}, chains_one, 1e5, 500, 1
)
saturation = 1e5 * 500 * EV_ERG * ARCMIN2_SR / (4 * np.pi)
close(thick[0], saturation, relative=2e-12)
assert thick_tau[0] > 100
assert (
    tau_one_radius(np.array([1.0, 2.0, 4.0]), np.array([4.0, 2.0, 0.5]))[1]
    == "resolved"
)
assert tau_one_radius(np.array([1.0, 2.0]), np.array([0.2, 0.1])) == (0.0, "below_grid")

# A radial PSF preserves a constant sky.
psf_radius = np.arange(0.5, 300, 1.0)
sigma = 6.0
density = np.exp(-0.5 * (psf_radius / sigma) ** 2) / (2 * np.pi * sigma**2)
source = np.arange(0.01, 30, 0.02)
operator = psf_operator.radial_operator(
    psf_radius,
    density,
    source,
    np.array([0, 15, 30, 60, 120, 300, 720]),
    CONFIG["quadrature"]["psf_radial_order"],
    CONFIG["quadrature"]["psf_azimuth_order"],
)
assert np.all(operator >= 0)
assert np.allclose(operator.sum(axis=1), 1, atol=2e-12)
assert np.allclose(operator @ np.ones_like(source), 1, atol=2e-12)
king = (1.56 - 1) / (np.pi * 5.5**2) * (1 + (psf_radius / 5.5) ** 2) ** -1.56
trial_source = np.arange(0.01, 30, 0.04)
edges = np.array([0, 15, 30, 45, 75, 135, 240, 360, 480, 600, 720])
trial = 1 / np.maximum(trial_source, 0.05) * np.exp(-trial_source / 20)
fine = (
    psf_operator.radial_operator(psf_radius, king, trial_source, edges, 4, 1536) @ trial
)
standard = (
    psf_operator.radial_operator(psf_radius, king, trial_source, edges, 4, 384) @ trial
)
assert np.max(np.abs(standard / fine - 1)) < 0.003

# Exposure-weighted annulus quadrature against direct image integration.
shape = (381, 381)
header = fits.Header(
    {
        "NAXIS": 2,
        "NAXIS1": shape[1],
        "NAXIS2": shape[0],
        "CTYPE1": "RA---TAN",
        "CTYPE2": "DEC--TAN",
        "CRPIX1": 191.0,
        "CRPIX2": 191.0,
        "CRVAL1": 10.0,
        "CRVAL2": 20.0,
        "CDELT1": -4 / 3600,
        "CDELT2": 4 / 3600,
        "RADESYS": "FK5",
    }
)
pixel_radius, pixel_area = psf_operator.image_radii(
    shape, header, np.array([190.0, 190.0])
)
yy, xx = np.indices(shape)
exposure_image = (100 + 2 * pixel_radius) * (
    1 + 0.15 * np.cos(np.arctan2(yy - 190, xx - 190))
)
with tempfile.TemporaryDirectory(dir=HERE) as folder:
    task04_header = header.copy()
    task04_header.update(
        DETECTOR="EPIC", BAND="soft", E_MIN=200, E_MAX=1000, COMOVING=True
    )
    task04_image = exposure_image.copy()
    task04_image[0, 0] = -1e-12
    path = Path(folder) / "task04.fits"
    fits.HDUList(
        [fits.PrimaryHDU(), fits.ImageHDU(task04_image, task04_header, name="EXPOSURE")]
    ).writeto(path)
    loaded_exposure, _ = psf_operator.task04_plane(path, (200, 1000))
    assert loaded_exposure[0, 0] == 0
exposure_radius, exposure_mean = psf_operator.radial_exposure(
    exposure_image, header, 12
)
quadrature = psf_operator.annulus_quadrature(edges, 4, exposure_radius, exposure_mean)
profile_function = lambda radius: np.exp(-radius / 9) * (1 + 0.02 * radius)
for index, (lo, hi) in enumerate(pairwise(edges)):
    use = (pixel_radius >= lo / 60) & (pixel_radius < hi / 60)
    direct = np.average(
        profile_function(pixel_radius[use]), weights=exposure_image[use]
    )
    radius_nodes, radial_weights = quadrature[index]
    assert math.isclose(
        radial_weights @ profile_function(radius_nodes), direct, rel_tol=3e-3
    )
manual_ring = np.searchsorted(edges / 60, pixel_radius, side="right") - 1
manual_use = (manual_ring >= 0) & (manual_ring < len(edges) - 1)
manual_weight = np.bincount(
    manual_ring[manual_use],
    exposure_image[manual_use] * pixel_area,
    minlength=len(edges) - 1,
)
assert np.allclose(
    psf_operator.annulus_exposure(exposure_image, header, edges),
    manual_weight,
    rtol=1e-3,
)
frame_density = np.vstack((king, 0.8 * king))
frame_weight = np.vstack((np.arange(1, 11), np.arange(10, 0, -1)))
composite = psf_operator.composite_density(frame_density, frame_weight)
for index in range(10):
    assert np.allclose(
        composite[index],
        np.average(frame_density, axis=0, weights=frame_weight[:, index]),
    )
exposure_operator = psf_operator.radial_operator(
    psf_radius, composite, trial_source, edges, 4, 384, exposure_radius, exposure_mean
)
assert np.allclose(exposure_operator @ np.ones_like(trial_source), 1, atol=2e-12)

with tempfile.TemporaryDirectory(dir=HERE) as folder:
    folder = Path(folder)
    radius = np.geomspace(0.1, 18, 30)
    profile = 0.02 / (radius + 0.2)
    products.write_tsv(
        folder / "radial.tsv",
        (
            "r_arcmin",
            "net_ct_s_arcmin2",
            "error_ct_s_arcmin2",
            "comet_model_ct_s_arcmin2",
        ),
        zip(radius, profile, 0.05 * profile, profile),
    )
    edges = np.asarray(run.SHARED["chemistry_annulus_edges_arcsec"], float)
    center = (edges[:-1] + edges[1:]) / 120
    prediction = 0.02 / (center + 0.2)
    data = prediction * (1 + 0.02 * np.sin(center))
    products.write_tsv(
        folder / "chemical.tsv",
        (
            "r_lo_arcsec",
            "r_hi_arcsec",
            "data_ct_s_arcmin2",
            "total_error_ct_s_arcmin2",
            "physical_model_ct_s_arcmin2",
            "prediction_ct_s_arcmin2",
            "residual_ct_s_arcmin2",
        ),
        zip(
            edges[:-1],
            edges[1:],
            data,
            0.05 * prediction,
            prediction / 1.1,
            prediction,
            data - prediction,
        ),
    )
    intrinsic_radius = np.arange(0.01, 30, 0.02)
    intrinsic = 0.03 / (intrinsic_radius + 0.1)
    products.write_tsv(
        folder / "intrinsic.tsv",
        ("radius_arcmin", "energy_flux_erg_cm2_s_arcmin2", "count_rate_ct_s_arcmin2"),
        zip(intrinsic_radius, intrinsic * 1e-13, intrinsic),
    )
    broken = {
        "normalization_at_break_ct_s_arcmin2": 0.02,
        "alpha_inner": 0.2,
        "alpha_outer": 1.0,
        "break_radius_arcmin": 0.63,
    }
    from matplotlib.axes import Axes

    original_plot, plotted = Axes.plot, {}

    def capture_plot(self, *args, **kwargs):
        lines = original_plot(self, *args, **kwargs)
        if kwargs.get("label") == "Intrinsic MAP solution":
            plotted["intrinsic"] = np.asarray(args[1])
        return lines

    Axes.plot = capture_plot
    try:
        figure.draw(
            folder / "radial.tsv",
            folder / "chemical.tsv",
            folder / "intrinsic.tsv",
            background=0,
            tau_radius=1.0,
            scale=1.1,
            km_per_arcmin=81920.5,
            broken=broken,
            output=folder,
            ecf=1.1536344645858954e-11,
        )
    finally:
        Axes.plot = original_plot
    assert np.allclose(plotted["intrinsic"], 1.1 * intrinsic)
    assert (folder / "figure3.png").stat().st_size > 10_000
    assert (folder / "figure3.pdf").stat().st_size > 10_000

    rng = np.random.default_rng(18)
    physical = np.zeros((2400, len(inference.PHYSICAL_NAMES)))
    index = {
        name: inference.PHYSICAL_NAMES.index(name) for name in inference.PHYSICAL_NAMES
    }
    physical[:, :3] = 10 ** rng.normal((28.45, 27.5, 28.0), 0.12, (2400, 3))
    ratio = np.r_[rng.normal(-2.0, 0.15, 1200), rng.normal(0.2, 0.15, 1200)]
    physical[:, index["log10_H2_ratio"]] = ratio
    physical[:, index["Q_H2"]] = 10**ratio * physical[:, :3].sum(axis=1)
    physical[:, index["beta_H2"]] = rng.uniform(0, 1, 2400)
    physical[:, index["phi_heavy"]] = 10 ** rng.normal(4.9, 0.12, 2400)
    physical[:, index["v_H2_km_s"]] = rng.normal(0.37, 0.04, 2400)
    physical[:, index["scale"]] = np.exp(rng.normal(0, 0.08, 2400))
    physical[:, index["offset"]] = rng.normal(0, 3e-5, 2400)
    physical[:, index["jitter"]] = abs(rng.normal(0, 0.1, 2400))
    diagnostic.draw(
        physical,
        -0.5 * rng.normal(size=(2400, 4)).var(axis=1),
        folder / "corner.png",
        limit=1000,
    )
    assert (folder / "corner.png").stat().st_size > 10_000
with tempfile.TemporaryDirectory(dir=HERE) as folder:
    reduction = Path(folder)
    band = run.SHARED["soft_band_ev"]
    low, high = band
    frames = ("pn-s003-p01", "mos1-s001-p01")
    manifest_rows, expected = [], {}
    for frame in frames:
        detector, exposure, _ = frame.split("-", 2)
        path = (
            reduction
            / "05_background"
            / frame
            / "soft"
            / f"{detector}{exposure.upper()}-expimsky-{low}-{high}.fits"
        )
        path.parent.mkdir(parents=True)
        fits.PrimaryHDU(np.ones((2, 2))).writeto(path)
        expected[frame] = path
        manifest_rows.extend(
            (
                {"frame": frame, "det": detector, "band": "soft", "exposure_map": path},
                {
                    "frame": frame,
                    "det": detector,
                    "band": "hard",
                    "exposure_map": "unused",
                },
            )
        )
    manifest = reduction / "05_background/background.tsv"
    columns = tuple(manifest_rows[0])
    products.write_tsv(
        manifest, columns, ([row[name] for name in columns] for row in manifest_rows)
    )
    assert psf.exposure_maps(reduction, frames, band) == expected
    fits.PrimaryHDU(np.zeros((2, 2))).writeto(expected[frames[0]], overwrite=True)
    try:
        psf.exposure_maps(reduction, frames, band)
    except RuntimeError:
        pass
    else:
        raise AssertionError("empty Task 01 exposure map accepted")
    fits.PrimaryHDU(np.ones((2, 2))).writeto(expected[frames[0]], overwrite=True)
    event = reduction / "04_sources/events/pn-s003-p01.fits"
    event.parent.mkdir(parents=True)
    event.write_bytes(b"x")
    assert psf.event_path(reduction, frames[0], {"event": event}) == event
    try:
        psf.event_path(reduction, frames[0], {"event": "wrong.fits"})
    except RuntimeError:
        pass
    else:
        raise AssertionError("noncanonical Task 01 event path accepted")
    invalid = (
        manifest_rows[:-1],
        [{**manifest_rows[0], "det": "mos2"}, *manifest_rows[1:]],
        [{**manifest_rows[0], "exposure_map": "wrong.fits"}, *manifest_rows[1:]],
    )
    for changed in invalid:
        products.write_tsv(
            manifest, columns, ([row[name] for name in columns] for row in changed)
        )
        try:
            psf.exposure_maps(reduction, frames, band)
        except RuntimeError:
            pass
        else:
            raise AssertionError("invalid PSF exposure-map manifest accepted")

with tempfile.TemporaryDirectory(dir=HERE) as folder:
    folder = Path(folder)
    reduction, scratch = folder / "outputs", folder / "psf"
    for stage in ("03_events", "04_sources", "07_comove"):
        (reduction / stage).mkdir(parents=True)
    frame = "pn-s003-p01"
    products.write_tsv(
        reduction / "03_events/frames.tsv", ("frame", "det"), ((frame, "pn"),)
    )
    products.write_tsv(
        reduction / "04_sources/events.tsv",
        ("frame", "det", "event"),
        ((frame, "pn", "event.fits"),),
    )
    products.write_tsv(
        reduction / "07_comove/ephemeris.tsv",
        ("mjd", "ra_deg", "dec_deg"),
        ((60000, 1, 2), (60001, 1.1, 2.1)),
    )
    (folder / "settings.json").write_text(
        json.dumps(
            {
                "ccf": "ccf",
                "summary": "summary.SAS",
                "expected_frames": [frame],
                "comove": {"map_step_s": 5},
            }
        )
    )
    originals = psf.exposure_maps, psf.frame_exposure_weights, psf.run_psfgen
    psf.exposure_maps = lambda *args: {frame: folder / "exposure.fits"}
    psf.frame_exposure_weights = lambda *args: (
        np.ones((1, 10)),
        np.array([[1.0, 2.0]]),
        np.zeros((3, 10)),
    )

    def fake_psfgen(output, image, detector, ra, dec, energy, work, sas_env):
        output = Path(output)
        header = fits.Header(
            {
                "CRPIX1": 4.0,
                "CRPIX2": 4.0,
                "CDELT1": -1 / 3600,
                "CDELT2": 1 / 3600,
                "CUNIT1": "deg",
                "CUNIT2": "deg",
            }
        )
        fits.PrimaryHDU(ellbeta, header).writeto(output)
        (Path(work) / ".pfiles").mkdir()
        (Path(work) / ".pfiles/psfgen.par").touch()
        return output

    psf.run_psfgen = fake_psfgen
    try:
        generated = psf.generate(
            reduction,
            {
                "detector_K_ct_s_per_erg_cm2_s": {"pn": 2.0, "mos1": 1.0, "mos2": 1.0},
                "mean_photon_energy_eV": 500.0,
            },
            scratch,
            run.SHARED["soft_band_ev"],
            np.arange(11.0),
        )
    finally:
        psf.exposure_maps, psf.frame_exposure_weights, psf.run_psfgen = originals
    assert generated[1].shape == (1, 4)
    assert list(scratch.iterdir()) == []

with tempfile.TemporaryDirectory(dir=HERE) as folder:
    path = Path(folder) / "psf.npz"
    dependency = Path(folder) / "input"
    dependency.touch()
    current_source = np.arange(*CONFIG["source_grid_arcmin"])
    current_edges = np.asarray(run.SHARED["chemistry_annulus_edges_arcsec"], float)
    radius = np.arange(0.5, 300, 1.0)
    base_density = np.exp(-0.5 * (radius / 6) ** 2) / (2 * np.pi * 6**2)
    base_density /= np.sum(base_density * 2 * np.pi * radius)
    density = np.tile(base_density, (10, 1))
    exposure_radius = np.linspace(0, 12, 49)
    exposure_mean = np.ones(49)
    expected_operator = psf_operator.radial_operator(
        radius,
        density,
        current_source,
        current_edges,
        CONFIG["quadrature"]["psf_radial_order"],
        CONFIG["quadrature"]["psf_azimuth_order"],
        exposure_radius,
        exposure_mean,
    )
    expected_frames = json.loads(
        (run.ROOT / "analysis_01_reduction/settings.json").read_text()
    )["expected_frames"]
    np.savez(
        path,
        operator=expected_operator,
        source_radius_arcmin=current_source,
        annulus_edges_arcsec=current_edges,
        density_arcsec2=density,
        radius_arcsec=radius,
        frame_weights=np.ones((22, 10)),
        frames=np.asarray(expected_frames),
        exposure_radius_arcmin=exposure_radius,
        exposure_mean_s=exposure_mean,
        detector_names=np.asarray(("pn", "mos1", "mos2")),
        detector_closure_fraction=np.zeros((3, 10)),
        task04_closure_fraction=np.zeros(10),
    )
    original_profile_exposure = run.profile_exposure
    run.profile_exposure = lambda edges: (
        exposure_radius,
        exposure_mean,
        np.full(10, 22.0),
    )
    try:
        assert np.array_equal(
            run.load_psf(path, {"edges": current_edges}, [dependency]),
            expected_operator,
        )
        with np.load(path, allow_pickle=False) as data:
            payload = {name: data[name] for name in data.files}
        for change, message in (
            ({"extra": np.array(0)}, "augmented PSF schema accepted"),
            ({"frames": payload["frames"][::-1]}, "wrong PSF frames accepted"),
            (
                {"operator": np.roll(payload["operator"], 1, axis=1)},
                "corrupt PSF operator accepted",
            ),
        ):
            np.savez(path, **{**payload, **change})
            try:
                run.load_psf(path, {"edges": current_edges}, [dependency])
            except RuntimeError:
                pass
            else:
                raise AssertionError(message)
        np.savez(path, **payload)
        try:
            run.load_psf(path, {"edges": current_edges + 1}, [dependency])
        except RuntimeError as error:
            assert "PSF grids" in str(error)
        else:
            raise AssertionError("stale PSF accepted")
        stamp = path.stat().st_mtime_ns
        os.utime(dependency, ns=(stamp + 1_000_000_000, stamp + 1_000_000_000))
        try:
            run.load_psf(path, {"edges": current_edges}, [dependency])
        except RuntimeError as error:
            assert "PSF inputs changed" in str(error)
        else:
            raise AssertionError("newer PSF input accepted")
    finally:
        run.profile_exposure = original_profile_exposure

with tempfile.TemporaryDirectory(dir=HERE) as folder:
    saved = (
        run.OUTPUT,
        run.load_inputs,
        run.psf_product,
        run.load_psf,
        run.psf_dependencies,
        run.corner_product,
        run.shutil.which,
        list(sys.argv),
    )
    calls = []
    run.OUTPUT = Path(folder)
    run.load_inputs = lambda: {"edges": np.arange(3.0)}
    run.psf_product = lambda inputs: calls.append("generate")
    run.psf_dependencies = lambda: calls.append("dependencies") or ()
    run.load_psf = lambda path, inputs, dependencies: calls.append("validate")
    run.corner_product = lambda: calls.append("corner")
    run.shutil.which = lambda name: "/sas/psfgen"
    try:
        preserved = [run.OUTPUT / name for name in ("chain.npz", "result.json")]
        for path in preserved:
            path.write_text("keep")
        sys.argv[:] = ["run.py", "corner"]
        run.main()
        assert calls == ["corner"]
        assert all(path.read_text() == "keep" for path in preserved)
        calls.clear()
        sys.argv[:] = ["run.py", "psf"]
        run.main()
        assert calls == ["generate", "dependencies", "validate"]
        calls.clear()
        run.shutil.which = lambda name: None
        try:
            run.main()
        except RuntimeError as error:
            assert str(error) == "source SAS before psf"
        else:
            raise AssertionError("missing SAS accepted")
        assert calls == []
    finally:
        (
            run.OUTPUT,
            run.load_inputs,
            run.psf_product,
            run.load_psf,
            run.psf_dependencies,
            run.corner_product,
            run.shutil.which,
            arguments,
        ) = saved
        sys.argv[:] = arguments

with tempfile.TemporaryDirectory(dir=HERE) as folder:
    saved = run.OUTPUT, run.diagnostic.draw
    run.OUTPUT = Path(folder)
    theta = np.zeros((3, len(inference.COORDS)))
    log_probability = np.arange(3.0)
    np.savez(
        run.OUTPUT / "chain.npz",
        samples=theta,
        coordinate_names=np.asarray(inference.COORDS),
        log_probability=log_probability,
    )
    captured = {}
    run.diagnostic.draw = lambda physical, lp, path: captured.update(
        physical=physical, log_probability=lp, path=path
    )
    try:
        run.corner_product()
        assert captured["physical"].shape == (3, len(inference.PHYSICAL_NAMES))
        assert np.array_equal(captured["log_probability"], log_probability)
        assert captured["path"] == run.OUTPUT / "corner.png"
        np.savez(
            run.OUTPUT / "chain.npz",
            samples=theta,
            coordinate_names=np.asarray(inference.COORDS[::-1]),
            log_probability=log_probability,
        )
        try:
            run.corner_product()
        except RuntimeError as error:
            assert str(error) == "chain.npz does not match the current model"
        else:
            raise AssertionError("mismatched corner chain accepted")
    finally:
        run.OUTPUT, run.diagnostic.draw = saved


# Priors, current rates, count conversion, offset, and jitter.
class StubTransport:
    orbit = type("O", (), {"current_r_au": 1.8})()
    source_radius = np.arange(3.0)

    @staticmethod
    def source_profile(q, beta, speed, chains, phi, energy, photons_per_ion):
        assert photons_per_ion == 1
        return np.array([2.0, 4.0, 6.0]), np.array([0.1, 0.2, 0.3])


theta = np.zeros(20)
for i, name in enumerate(("H2O", "CO2", "CO")):
    theta[i] = math.log(CONFIG["priors"][f"Qref_{name}"][0])
theta[3:6] = [3.3, 0.392, 1.3]
theta[6:8] = [0.5, 0.0]
theta[8] = 1e5
theta[9:13] = 0.37
theta[13:17] = [1.0, 1.0, 1.0, 10.0]
theta[17:] = [math.log(1.1), 0.2, 0.1]
posterior = inference.Posterior(
    np.array([1.3, 2.4]),
    np.array([0.2, 0.3]),
    2.0,
    StubTransport(),
    {},
    np.array([[1, 0, 0], [0, 1, 0]]),
    CONFIG,
    500,
    0.05,
    0.4,
    1e5,
    3.8e4,
)
assert math.isfinite(posterior.log_prior(theta))
bad = theta.copy()
bad[7] = 2
assert posterior.log_prior(bad) == -math.inf
walkers = inference.draw_walkers(
    posterior, CONFIG["mcmc"]["walkers"], np.random.default_rng(CONFIG["mcmc"]["seed"])
)
initial_ratio = walkers[:, inference.COORDS.index("log10_H2_ratio")]
assert np.any(initial_ratio <= -1.5) and np.any(initial_ratio >= -0.5)
mode_chain = inference.run_sampler(posterior, walkers, 20, 40, 1, seed=10)[0]
retained_ratio = mode_chain[:, :, inference.COORDS.index("log10_H2_ratio")]
assert np.any(retained_ratio <= -1.5) and np.any(retained_ratio >= -0.5)
prediction, variance, model, _, _, physical = posterior.predict(theta)
assert np.allclose(model, [1, 2])
assert np.allclose(prediction, [1.3, 2.4])
assert np.allclose(variance, [0.2**2 + (0.1 * 1.1) ** 2, 0.3**2 + (0.1 * 2.2) ** 2])
logpost, loglike = posterior.evaluate(theta)
close(loglike, -0.5 * np.log(2 * np.pi * variance).sum(), relative=2e-12)
assert len(physical) == len(inference.PHYSICAL_NAMES) == 21
q = dict(zip(inference.PHYSICAL_NAMES, physical))
close(q["Q_H2"] / (q["Q_H2O"] + q["Q_CO2"] + q["Q_CO"]), 1.0)

# Count-space data equation.
events, qpb, soft_proton, oot, exposure, background = (
    120.0,
    10.0,
    5.0,
    20.0,
    1000.0,
    0.01,
)
data = (events - qpb - soft_proton - 0.063 * oot) / exposure - background
error = math.sqrt(events + qpb + soft_proton + 0.063**2 * oot) / exposure
close(data, 0.09374)
close(error, math.sqrt(135.07938) / 1000)

with tempfile.TemporaryDirectory(dir=HERE) as folder:
    folder = Path(folder)
    profile = folder / "profile"
    profile.mkdir()
    background_prior = {"mean_ct_s_arcmin2": 0.01, "sigma_ct_s_arcmin2": 0.002}
    broken = {
        "normalization_at_break_ct_s_arcmin2": 0.02,
        "alpha_inner": 0.2,
        "alpha_outer": 1.0,
        "break_radius_arcmin": 0.63,
    }
    (profile / "result.json").write_text(
        json.dumps(
            {
                "band_eV": run.SHARED["soft_band_ev"],
                "epic": {"background_prior": background_prior, "broken": broken},
            }
        )
    )
    edges = np.asarray(run.SHARED["chemistry_annulus_edges_arcsec"], float)
    products.write_tsv(
        profile / "annuli.tsv",
        (
            "r_lo_arcmin",
            "r_hi_arcmin",
            "exposure_s_arcmin2",
            "events_count",
            "qpb_count",
            "sp_count",
            "oot_count",
        ),
        [
            [edges[i] / 60, edges[i + 1] / 60, 1000, 120 + i, 10, 5, 20]
            for i in range(10)
        ],
    )
    figure_profile = profile / "profile.tsv"
    figure_profile.write_text("paper-facing radial profile\n")
    paths = [folder / f"{name}.json" for name in ("spectrum", "wind", "geometry")]
    paths[0].write_text(
        json.dumps(
            {
                "marker": "spectrum",
                "fit_band_keV": [value / 1000 for value in run.SHARED["soft_band_ev"]],
            }
        )
    )
    paths[1].write_text(
        json.dumps(
            {
                "f107": {"f107p_sfu": 175.25},
                "heavy_ion_prior": run.SHARED["heavy_ion_prior"],
            }
        )
    )
    paths[2].write_text(json.dumps({"km_per_arcmin": 81920.5}))
    image = folder / "image.fits"
    image.touch()
    stamp = max(path.stat().st_mtime_ns for path in (*paths, image)) + 2_000_000_000
    for path in (profile / "result.json", profile / "annuli.tsv", figure_profile):
        os.utime(path, ns=(stamp, stamp))
    original = (run.PROFILE, run.SPECTRUM, run.WIND, run.GEOMETRY, run.MORPHOLOGY_IMAGE)
    run.PROFILE, run.SPECTRUM, run.WIND, run.GEOMETRY, run.MORPHOLOGY_IMAGE = (
        profile,
        *paths,
        image,
    )
    try:
        loaded = run.load_inputs()
        figure_profile.unlink()
        try:
            run.load_inputs()
        except RuntimeError as error:
            assert "profile predates" in str(error)
        else:
            raise AssertionError("missing Task 05 figure profile accepted")
        figure_profile.write_text("paper-facing radial profile\n")
        input_stamp = min(path.stat().st_mtime_ns for path in (*paths, image))
        os.utime(
            figure_profile,
            ns=(input_stamp - 1_000_000_000, input_stamp - 1_000_000_000),
        )
        try:
            run.load_inputs()
        except RuntimeError as error:
            assert "profile predates" in str(error)
        else:
            raise AssertionError("stale Task 05 figure profile accepted")
        os.utime(figure_profile, ns=(stamp, stamp))
        bad_edges = edges.copy()
        bad_edges[1] += 1
        products.write_tsv(
            profile / "annuli.tsv",
            (
                "r_lo_arcmin",
                "r_hi_arcmin",
                "exposure_s_arcmin2",
                "events_count",
                "qpb_count",
                "sp_count",
                "oot_count",
            ),
            [
                [bad_edges[i] / 60, bad_edges[i + 1] / 60, 1000, 120 + i, 10, 5, 20]
                for i in range(10)
            ],
        )
        try:
            run.load_inputs()
        except RuntimeError as error:
            assert "chemistry annuli" in str(error)
        else:
            raise AssertionError("stale chemistry annuli accepted")
        stale = (
            max(
                (profile / "result.json").stat().st_mtime_ns,
                (profile / "annuli.tsv").stat().st_mtime_ns,
            )
            + 2_000_000_000
        )
        os.utime(paths[0], ns=(stale, stale))
        try:
            run.load_inputs()
        except RuntimeError as error:
            assert "profile predates" in str(error)
        else:
            raise AssertionError("stale Task 05 profile accepted")
    finally:
        run.PROFILE, run.SPECTRUM, run.WIND, run.GEOMETRY, run.MORPHOLOGY_IMAGE = (
            original
        )
    close(loaded["data"][0], 0.09374)
    close(loaded["error"][0], math.sqrt(135.07938) / 1000)
    assert loaded["f107p"] == 175.25 and loaded["geometry"]["km_per_arcmin"] == 81920.5


# Posterior and recoil-product serialization.
class ProductPosterior:
    config = json.loads(json.dumps(CONFIG))
    transport = type(
        "T",
        (),
        {
            "orbit": type("O", (), {"current_r_au": 1.8})(),
            "source_radius": np.array([0.1, 1.0, 2.0]),
            "network": {"OH": {"sigma_cm2": 4.81e-15}},
        },
    )()
    ecf = 2.0

    def predict(self, value):
        factor = 1 + value[0] - theta[0]
        model = factor * np.linspace(1, 2, 10)
        source = factor * np.array([3.0, 2.0, 1.0])
        tau = factor * np.array([2.0, 1.0, 0.5])
        physical_value = inference.unpack(value, CONFIG, self.transport.orbit)[-1]
        return model, np.ones(10), model, source, tau, physical_value


shape = (101, 101, len(inference.COORDS))
product_chain = np.broadcast_to(theta, shape).copy()
flat_product = product_chain.reshape(-1, shape[-1])
sequence = np.arange(len(flat_product))
map_index = 500
flat_product[:, 0] += sequence * 1e-7
flat_product[:, 9] = 0.3 + sequence * 1e-6
ratio_index = inference.COORDS.index("log10_H2_ratio")
flat_product[:, ratio_index] = 0.1
flat_product[0, ratio_index] = -2
product_logl = sequence.reshape(shape[:2]).astype(float)
product_logp = -abs(sequence - map_index).reshape(shape[:2]).astype(float)
product_inputs = {
    "edges": np.asarray(run.SHARED["chemistry_annulus_edges_arcsec"], float),
    "data": np.linspace(1, 2, 10),
    "error": np.full(10, 0.1),
    "background": {"mean_ct_s_arcmin2": 0.01, "sigma_ct_s_arcmin2": 0.002},
    "f107p": 175.25,
    "geometry": {"km_per_arcmin": 81920.5},
    "broken": broken,
}
with tempfile.TemporaryDirectory(dir=HERE) as folder:
    original_output, original_profile = products.OUTPUT, products.PROFILE
    original_draw = products.figure.draw
    original_corner = products.diagnostic.draw
    figure_arguments = {}
    corner_arguments = {}
    products.OUTPUT, products.PROFILE = Path(folder), Path(folder)
    products.figure.draw = lambda *args, **kwargs: figure_arguments.update(kwargs)
    products.diagnostic.draw = (
        lambda physical, probability, path: corner_arguments.update(
            physical=physical, probability=probability, path=path
        )
    )
    try:
        products.write(
            product_inputs,
            ProductPosterior(),
            rows,
            product_chain,
            product_logp,
            product_logl,
            0.25,
        )
    finally:
        products.OUTPUT, products.PROFILE = original_output, original_profile
        products.figure.draw = original_draw
        products.diagnostic.draw = original_corner
    with np.load(Path(folder) / "chain.npz", allow_pickle=False) as saved:
        assert tuple(saved["physical_names"]) == inference.PHYSICAL_NAMES
        assert len(saved["physical_draws"]) == len(flat_product)
        physical_all = np.array(
            [
                inference.unpack(value, CONFIG, ProductPosterior.transport.orbit)[-1]
                for value in flat_product
            ]
        )
        selected = physical_all.astype("f4")
        assert np.array_equal(saved["physical_draws"], selected)
        assert saved["physical_draws"][0, 8] == -2
        assert tuple(saved["best_labels"]) == (
            "overall_ml",
            "high_ml",
            "low_ml",
            "high_map",
            "low_map",
        )
        best_index = np.array(
            [len(flat_product) - 1, len(flat_product) - 1, 0, map_index, 0]
        )
        assert saved["best_coordinates"].dtype == np.float64
        assert np.array_equal(saved["best_coordinates"], flat_product[best_index])
        assert np.array_equal(
            saved["best_log_probability"], product_logp.ravel()[best_index]
        )
        assert np.array_equal(
            saved["best_log_likelihood"], product_logl.ravel()[best_index]
        )
    saved_result = json.loads((Path(folder) / "result.json").read_text())
    assert "posterior" in saved_result and "high_H2" not in saved_result
    assert saved_result["mean_acceptance_fraction"] == 0.25
    assert saved_result["sampled_delta_log_likelihood_high_minus_low"] == 10200
    assert saved_result["sampled_delta_log_probability_high_minus_low"] == 500
    assert saved_result["photons_per_captured_ion"] == 1
    assert set(saved_result["maximum_likelihood_high"]) == {
        "log_probability",
        "log_likelihood",
        "Q_total_s-1",
        "Q_H2_s-1",
        "H2_fraction",
    }
    assert saved_result["maximum_a_posteriori_high"]["log_probability"] == 0
    assert saved_result["maximum_a_posteriori_low"]["log_probability"] == -500
    assert corner_arguments["path"].name == "corner.png"
    assert corner_arguments["physical"].shape == (len(flat_product), 21)
    assert corner_arguments["probability"].shape == (len(flat_product),)
    saved_profile = figure.table(Path(folder) / "profile.tsv")
    expected_map_model = ProductPosterior().predict(flat_product[map_index])[2]
    assert np.array_equal(
        saved_profile["physical_model_ct_s_arcmin2"], expected_map_model
    )
    with (Path(folder) / "posterior.tsv").open(newline="") as stream:
        saved_posterior = list(csv.DictReader(stream, delimiter="\t"))
    expected_map = inference.unpack(
        flat_product[map_index], CONFIG, ProductPosterior.transport.orbit
    )[-1]
    assert (
        figure_arguments["scale"]
        == expected_map[inference.PHYSICAL_NAMES.index("scale")]
    )
    q_map = expected_map[:4]
    expected_map = np.r_[expected_map, q_map.sum(), q_map[3] / q_map.sum()]
    assert np.array_equal([float(row["map"]) for row in saved_posterior], expected_map)
    q_all = physical_all[:, :4]
    expected_values = np.c_[
        physical_all, q_all.sum(axis=1), q_all[:, 3] / q_all.sum(axis=1)
    ]
    for i, row in enumerate(saved_posterior):
        assert np.allclose(
            [float(row[name]) for name in ("q16", "median", "q84")],
            np.quantile(expected_values[:, i], (0.16, 0.5, 0.84)),
        )
assert CONFIG["fixed_velocity_km_s"] == {
    "heavy_daughter": 1.685,
    "photolytic_H2": 11.335,
}


# Short deterministic ensemble recovery.
class Toy:
    @staticmethod
    def log_probability(value):
        likelihood = -0.5 * float(np.sum(np.asarray(value) ** 2))
        return likelihood, likelihood


rng = np.random.default_rng(8)
start = rng.normal(0.35, 0.15, (44, len(inference.COORDS)))
chain, _, blobs, acceptance = inference.run_sampler(Toy(), start, 80, 180, 1, seed=8)
repeat = inference.run_sampler(Toy(), start, 80, 180, 1, seed=8)[0]
assert np.array_equal(chain, repeat)
means = np.mean(chain[-80:], axis=(0, 1))
assert np.mean(np.abs(means)) < 0.15
assert np.linalg.norm(means) < np.linalg.norm(np.mean(start, axis=0))
assert np.all(np.isfinite(blobs))
assert 0.1 < acceptance < 0.8

# Half-normal prior reweighting favors small jitter when its scale shrinks.
draws = np.tile(physical, (101, 1))
jitter_grid = np.linspace(0, 0.3, len(draws))
draws[:, inference.PHYSICAL_NAMES.index("jitter")] = jitter_grid
draws[:, inference.PHYSICAL_NAMES.index("Q_H2")] = 1e28 * (1 + jitter_grid)
checks = products.jitter_checks(draws)
assert all(
    set(value) == {"Q_H2_median_s-1", "H2_fraction_median", "effective_sample_size"}
    for value in checks.values()
)
for sigma in (0.05, 0.20):
    log_weight = np.log(CONFIG["priors"]["jitter_sigma"] / sigma) - 0.5 * (
        jitter_grid**2 * (sigma**-2 - CONFIG["priors"]["jitter_sigma"] ** -2)
    )
    weight = np.exp(log_weight - log_weight.max())
    effective = weight.sum() ** 2 / (weight @ weight)
    assert math.isclose(
        checks[f"sigma_{sigma:.2f}"]["effective_sample_size"], effective
    )
    assert effective < len(draws)
assert checks["sigma_0.05"]["Q_H2_median_s-1"] < checks["sigma_0.20"]["Q_H2_median_s-1"]
