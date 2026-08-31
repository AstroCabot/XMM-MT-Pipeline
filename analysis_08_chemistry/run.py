import argparse
import csv
import json
import os
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RUNTIME = ROOT / "work" / "runtime" / "analysis_08_chemistry"
for name in ("tmp", "cache", "config", "matplotlib"):
    (RUNTIME / name).mkdir(parents=True, exist_ok=True)
os.environ.update(
    TMPDIR=str(RUNTIME / "tmp"),
    XDG_CACHE_HOME=str(RUNTIME / "cache"),
    XDG_CONFIG_HOME=str(RUNTIME / "config"),
    MPLCONFIGDIR=str(RUNTIME / "matplotlib"),
)

import numpy as np

import chemistry
import diagnostic
import inference
import products
import psf
import psf_operator
from config import CONFIG, SHARED
from transport import Orbit, Transport

OUTPUT = HERE / "outputs"
WORK = ROOT / "work" / "analysis_08_chemistry"
REDUCTION = ROOT / "analysis_01_reduction" / "outputs"
PROFILE = ROOT / "analysis_05_radial_profile" / "outputs"
SPECTRUM = ROOT / "analysis_03_spectrum" / "outputs" / "result.json"
WIND = ROOT / "analysis_06_solar_wind" / "outputs" / "result.json"
MORPHOLOGY = ROOT / "analysis_04_morphology" / "outputs"
GEOMETRY = MORPHOLOGY / "result.json"
MORPHOLOGY_IMAGE = MORPHOLOGY / "image.fits"


def read_tsv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def require_current_profile():
    products = (
        PROFILE / "result.json",
        PROFILE / "annuli.tsv",
        PROFILE / "profile.tsv",
    )
    inputs = (SPECTRUM, GEOMETRY, MORPHOLOGY_IMAGE)
    if any(not path.is_file() for path in products + inputs) or max(
        path.stat().st_mtime_ns for path in inputs
    ) > min(path.stat().st_mtime_ns for path in products):
        raise RuntimeError("Task 05 profile predates its calibration inputs")


def load_inputs():
    require_current_profile()
    radial = json.loads((PROFILE / "result.json").read_text())
    spectrum = json.loads(SPECTRUM.read_text())
    wind = json.loads(WIND.read_text())
    geometry = json.loads(GEOMETRY.read_text())
    fit_band = np.asarray(spectrum.get("fit_band_keV", []), float)
    if (
        radial.get("band_eV") != SHARED["soft_band_ev"]
        or fit_band.shape != (2,)
        or not np.allclose(fit_band * 1000, SHARED["soft_band_ev"], rtol=0, atol=1e-9)
    ):
        raise RuntimeError("upstream products differ from the shared soft band")
    rows = read_tsv(PROFILE / "annuli.tsv")
    expected = np.asarray(SHARED["chemistry_annulus_edges_arcsec"], float)
    if len(rows) != len(expected) - 1:
        raise RuntimeError("Task 05 chemistry annuli do not match the shared edges")
    lower = np.array([float(row["r_lo_arcmin"]) for row in rows]) * 60
    upper = np.array([float(row["r_hi_arcmin"]) for row in rows]) * 60
    edges = np.r_[lower[0], upper]
    if not np.allclose(lower, expected[:-1]) or not np.allclose(upper, expected[1:]):
        raise RuntimeError("Task 05 chemistry annuli do not match the shared edges")
    exposure = np.array([float(row["exposure_s_arcmin2"]) for row in rows])
    events = np.array([float(row["events_count"]) for row in rows])
    qpb = np.array([float(row["qpb_count"]) for row in rows])
    sp = np.array([float(row["sp_count"]) for row in rows])
    oot = np.array([float(row["oot_count"]) for row in rows])
    background = radial["epic"]["background_prior"]
    data = (events - qpb - sp - SHARED["pn_oot_fraction"] * oot) / exposure
    data -= background["mean_ct_s_arcmin2"]
    error = np.sqrt(events + qpb + sp + SHARED["pn_oot_fraction"] ** 2 * oot) / exposure
    heavy = SHARED["heavy_ion_prior"]
    if wind["heavy_ion_prior"] != heavy:
        raise RuntimeError("Task 06 heavy-ion prior differs from shared parameters")
    return {
        "rows": rows,
        "edges": edges,
        "data": data,
        "error": error,
        "background": background,
        "spectrum": spectrum,
        "wind": wind,
        "geometry": geometry,
        "f107p": float(wind["f107"]["f107p_sfu"]),
        "heavy": heavy,
        "broken": radial["epic"]["broken"],
    }


def make_orbit():
    return Orbit(
        **{
            name: CONFIG["orbit"][key]
            for name, key in (
                ("q_au", "perihelion_q_au"),
                ("eccentricity", "eccentricity"),
                ("perihelion_jd", "perihelion_jd_tdb"),
                ("observation_jd", "observation_jd_tdb"),
            )
        },
        maximum_age_s=1.02 * CONFIG["shells"][2] / CONFIG["priors"]["O_velocity"][0],
    )


def build_model(inputs, operator):
    source = np.arange(*CONFIG["source_grid_arcmin"])
    shell_count, inner, outer = CONFIG["shells"]
    orbit = make_orbit()
    rows = chemistry.reactions(HERE / "reactions.tsv", inputs["f107p"])
    network = chemistry.network(rows, CONFIG["cross_section_cm2"])
    quadrature = CONFIG["quadrature"]
    transport = Transport(
        orbit,
        network,
        shell_count,
        inner,
        outer,
        source,
        inputs["geometry"]["km_per_arcmin"],
        quadrature["granddaughter_shell_stride"],
        quadrature["granddaughter_radial_stride"],
    )
    chains = {
        name: chemistry.chains(name, network) for name in ("H2O", "CO2", "CO", "H2")
    }
    spectrum, heavy = inputs["spectrum"], inputs["heavy"]
    posterior = inference.Posterior(
        inputs["data"],
        inputs["error"],
        spectrum["mos2_ecf_erg_cm2_per_count"],
        transport,
        chains,
        operator,
        CONFIG,
        spectrum["mean_photon_energy_eV"],
        spectrum["mos2_ecf_sigma_ln"],
        inputs["background"]["sigma_ct_s_arcmin2"],
        heavy["location"],
        heavy["scale"],
    )
    return posterior, rows


def corner_product():
    path = OUTPUT / "chain.npz"
    if not path.is_file():
        raise RuntimeError("corner plot requires chain.npz")
    with np.load(path, allow_pickle=False) as stored:
        theta = np.asarray(stored["samples"], float)
        log_probability = np.asarray(stored["log_probability"], float)
        names = tuple(stored["coordinate_names"].astype(str))
    if (
        names != inference.COORDS
        or theta.ndim != 2
        or theta.shape[1] != len(names)
        or log_probability.shape != (len(theta),)
    ):
        raise RuntimeError("chain.npz does not match the current model")
    orbit = make_orbit()
    physical = np.array([inference.unpack(row, CONFIG, orbit)[-1] for row in theta])
    diagnostic.draw(physical, log_probability, OUTPUT / "corner.png")


def profile_exposure(edges):
    exposure, header = psf_operator.task04_plane(
        MORPHOLOGY_IMAGE, SHARED["soft_band_ev"]
    )
    radius, mean = psf_operator.radial_exposure(exposure, header, edges[-1] / 60)
    annuli = psf_operator.annulus_exposure(exposure, header, edges)
    return radius, mean, annuli


def psf_product(inputs):
    radius, frame_density, weights, frames, detector_closure = psf.generate(
        REDUCTION,
        inputs["spectrum"],
        WORK / "psf",
        SHARED["soft_band_ev"],
        inputs["edges"],
    )
    density = psf_operator.composite_density(frame_density, weights)
    exposure_radius, exposure_mean, expected_weights = profile_exposure(inputs["edges"])
    task04_closure = psf_operator.fractional_closure(
        weights.sum(axis=0), expected_weights, "Task 04 MOS2-equivalent"
    )
    source = np.arange(*CONFIG["source_grid_arcmin"])
    quadrature = CONFIG["quadrature"]
    operator = psf_operator.radial_operator(
        radius,
        density,
        source,
        inputs["edges"],
        quadrature["psf_radial_order"],
        quadrature["psf_azimuth_order"],
        exposure_radius,
        exposure_mean,
    )
    np.savez_compressed(
        OUTPUT / "psf.npz",
        radius_arcsec=radius,
        density_arcsec2=density,
        operator=operator,
        source_radius_arcmin=source,
        frame_weights=weights,
        frames=frames,
        annulus_edges_arcsec=inputs["edges"],
        exposure_radius_arcmin=exposure_radius,
        exposure_mean_s=exposure_mean,
        detector_names=np.asarray(tuple(psf.INSTRUMENT)),
        detector_closure_fraction=detector_closure,
        task04_closure_fraction=task04_closure,
    )


def psf_dependencies():
    frame_rows = read_tsv(REDUCTION / "03_events/frames.tsv")
    event_rows = read_tsv(REDUCTION / "04_sources/events.tsv")
    frames = [row["frame"] for row in frame_rows]
    events = {row["frame"]: row for row in event_rows}
    settings_path = ROOT / "analysis_01_reduction/settings.json"
    settings = json.loads(settings_path.read_text())
    summary = REDUCTION / "01_init/odf" / Path(settings["summary"]).name
    expected = set(settings["expected_frames"])
    maps = psf.exposure_maps(REDUCTION, expected, SHARED["soft_band_ev"])
    stacks = psf.detector_exposure_stacks(REDUCTION, SHARED["soft_band_ev"])
    if (
        len(frames) != len(expected)
        or set(frames) != expected
        or len(event_rows) != len(expected)
        or set(events) != expected
    ):
        raise RuntimeError(
            f"PSF inputs do not contain the {len(expected)} retained frames"
        )
    for row in frame_rows:
        detector = psf.frame_detector(row["frame"])
        if row["det"] != detector or events[row["frame"]]["det"] != detector:
            raise RuntimeError(f"Task 01 detector mismatch for {row['frame']}")
    paths = [
        HERE / "run.py",
        HERE / "psf.py",
        HERE / "psf_operator.py",
        HERE / "config.py",
        HERE / "parameters.json",
        ROOT / "parameters.json",
        SPECTRUM,
        PROFILE / "annuli.tsv",
        MORPHOLOGY_IMAGE,
        settings_path,
        REDUCTION / "01_init/ccf.cif",
        summary,
        REDUCTION / "03_events/frames.tsv",
        REDUCTION / "04_sources/events.tsv",
        REDUCTION / "05_background/background.tsv",
        REDUCTION / "07_comove/ephemeris.tsv",
    ]
    paths += [
        psf.event_path(REDUCTION, row["frame"], events[row["frame"]])
        for row in frame_rows
    ]
    paths += list(maps.values())
    paths += [item[0] for item in stacks.values()]
    return paths


def load_psf(path, inputs, dependencies):
    path = Path(path)
    dependencies = [Path(item) for item in dependencies]
    if (
        not path.is_file()
        or any(not item.is_file() for item in dependencies)
        or any(
            item.stat().st_mtime_ns > path.stat().st_mtime_ns for item in dependencies
        )
    ):
        raise RuntimeError("PSF inputs changed; rerun psf")
    with np.load(path, allow_pickle=False) as data:
        source = np.arange(*CONFIG["source_grid_arcmin"])
        annuli = len(inputs["edges"]) - 1
        expected_frames = tuple(
            json.loads((ROOT / "analysis_01_reduction/settings.json").read_text())[
                "expected_frames"
            ]
        )
        frame_count = len(expected_frames)
        detector_count = len(psf.INSTRUMENT)
        required = {
            "radius_arcsec",
            "density_arcsec2",
            "operator",
            "source_radius_arcmin",
            "frame_weights",
            "frames",
            "annulus_edges_arcsec",
            "exposure_radius_arcmin",
            "exposure_mean_s",
            "detector_names",
            "detector_closure_fraction",
            "task04_closure_fraction",
        }
        if set(data.files) != required:
            raise RuntimeError("PSF cache does not match the current schema")
        operator = np.asarray(data["operator"], float)
        radius = np.asarray(data["radius_arcsec"], float)
        density = np.asarray(data["density_arcsec2"], float)
        weights = np.asarray(data["frame_weights"], float)
        frames = tuple(data["frames"].astype(str))
        exposure_radius = np.asarray(data["exposure_radius_arcmin"], float)
        exposure_mean = np.asarray(data["exposure_mean_s"], float)
        detector_closure = np.asarray(data["detector_closure_fraction"], float)
        task04_closure = np.asarray(data["task04_closure_fraction"], float)
        probability = density * (2 * np.pi * radius)[None, :]
        current_radius, current_mean, current_annuli = profile_exposure(inputs["edges"])
        quadrature = CONFIG["quadrature"]
        rebuilt = psf_operator.radial_operator(
            radius,
            density,
            source,
            inputs["edges"],
            quadrature["psf_radial_order"],
            quadrature["psf_azimuth_order"],
            exposure_radius,
            exposure_mean,
        )
        if (
            not np.array_equal(data["source_radius_arcmin"], source)
            or not np.array_equal(data["annulus_edges_arcsec"], inputs["edges"])
            or radius.ndim != 1
            or len(radius) < 2
            or not np.array_equal(radius, np.arange(len(radius)) + 0.5)
        ):
            raise RuntimeError("PSF grids do not match the current inputs")
        if (
            operator.shape != (annuli, len(source))
            or np.any(~np.isfinite(operator))
            or np.any(operator < 0)
            or not np.allclose(operator.sum(axis=1), 1, rtol=0, atol=2e-12)
            or not np.allclose(operator, rebuilt, rtol=2e-12, atol=2e-14)
        ):
            raise RuntimeError("invalid cached PSF operator")
        if (
            density.shape != (annuli, len(radius))
            or np.any(~np.isfinite(density))
            or np.any(density < 0)
            or np.any(probability.sum(axis=1) <= 0.99)
            or np.any(probability.sum(axis=1) > 1 + 1e-10)
        ):
            raise RuntimeError("invalid cached PSF density")
        if (
            weights.shape != (frame_count, annuli)
            or np.any(~np.isfinite(weights))
            or np.any(weights <= 0)
            or frames != expected_frames
        ):
            raise RuntimeError("invalid cached PSF frame weights")
        if (
            exposure_radius.ndim != 1
            or len(exposure_radius) < 2
            or exposure_mean.shape != exposure_radius.shape
            or np.any(~np.isfinite(exposure_radius))
            or np.any(np.diff(exposure_radius) <= 0)
            or np.any(~np.isfinite(exposure_mean))
            or np.any(exposure_mean < 0)
            or not np.allclose(exposure_radius, current_radius, rtol=0, atol=1e-12)
            or not np.allclose(exposure_mean, current_mean, rtol=2e-12, atol=2e-12)
        ):
            raise RuntimeError("invalid cached PSF exposure profile")
        if (
            not np.allclose(
                weights.sum(axis=0),
                current_annuli,
                rtol=psf_operator.ANNULUS_CLOSURE_RTOL,
                atol=psf_operator.ANNULUS_CLOSURE_RTOL,
            )
            or not np.array_equal(
                data["detector_names"], np.asarray(tuple(psf.INSTRUMENT))
            )
            or detector_closure.shape != (detector_count, annuli)
            or task04_closure.shape != (annuli,)
            or np.any(~np.isfinite(detector_closure))
            or np.any(~np.isfinite(task04_closure))
            or np.max(np.abs(detector_closure)) > psf_operator.ANNULUS_CLOSURE_RTOL
            or np.max(np.abs(task04_closure)) > psf_operator.ANNULUS_CLOSURE_RTOL
        ):
            raise RuntimeError("cached PSF exposure closure failed")
        return operator


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "stage", choices=("all", "psf", "fit", "corner"), nargs="?", default="all"
    )
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "corner":
        corner_product()
        return
    for name in products.FINAL_PRODUCTS:
        (OUTPUT / name).unlink(missing_ok=True)
    if args.stage in ("all", "psf"):
        (OUTPUT / "psf.npz").unlink(missing_ok=True)
    inputs = load_inputs()
    if args.stage in ("all", "psf"):
        if shutil.which("psfgen") is None:
            raise RuntimeError("source SAS before psf")
        psf_product(inputs)
    operator = load_psf(OUTPUT / "psf.npz", inputs, psf_dependencies())
    if args.stage == "psf":
        return
    posterior, reactions = build_model(inputs, operator)
    settings = CONFIG["mcmc"]
    rng = np.random.default_rng(settings["seed"])
    start = inference.draw_walkers(posterior, settings["walkers"], rng)
    chain, logp, logl, acceptance = inference.run_sampler(
        posterior,
        start,
        settings["burn_steps"],
        settings["retained_steps"],
        settings["processes"],
        settings["seed"],
    )
    products.write(inputs, posterior, reactions, chain, logp, logl, acceptance)


if __name__ == "__main__":
    main()
