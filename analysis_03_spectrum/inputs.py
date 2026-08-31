import csv
import json
from pathlib import Path

import numpy as np

import extract

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REDUCTION = ROOT / "analysis_01_reduction" / "outputs"
OUTPUT = HERE / "outputs"
WORK = ROOT / "work" / "analysis_03_spectrum"
PARAMS = json.loads((HERE / "parameters.json").read_text())
SHARED = json.loads((ROOT / "parameters.json").read_text())
SOFTWARE = json.loads((ROOT / "analysis_01_reduction/settings.json").read_text())
FIT_BAND_KEV = tuple(float(value) / 1000 for value in SHARED["soft_band_ev"])
REGIONS = {
    "inner": (0, float(SHARED["science_aperture_radius_arcsec"])),
    "outer": tuple(map(float, SHARED["spectrum_outer_annulus_arcsec"])),
}


def read_tsv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def exact_index(rows, key, expected, label):
    keys = [key(row) for row in rows]
    if len(keys) != len(expected) or set(keys) != expected:
        raise RuntimeError(f"Task 01 {label} table is incomplete or duplicated")
    return dict(zip(keys, rows))


def frame_detector(frame):
    detector = frame.split("-", 1)[0]
    if detector not in ("mos1", "mos2", "pn"):
        raise RuntimeError(f"invalid EPIC frame name: {frame}")
    return detector


def load_inputs():
    needed = ("05_background", "06_soft_protons", "07_comove")
    if any(not (REDUCTION / name / ".done").exists() for name in needed):
        raise RuntimeError(
            "Task 01 background, soft-proton, and comet products must finish first"
        )
    expected = set(SOFTWARE["expected_frames"])
    frame_rows = read_tsv(REDUCTION / "03_events/frames.tsv")
    exact_index(frame_rows, lambda row: row["frame"], expected, "frame")
    events = exact_index(
        read_tsv(REDUCTION / "04_sources/events.tsv"),
        lambda row: row["frame"],
        expected,
        "event",
    )
    pairs = {(frame, band) for frame in expected for band in ("soft", "hard")}
    backgrounds = exact_index(
        read_tsv(REDUCTION / "05_background/background.tsv"),
        lambda row: (row["frame"], row["band"]),
        pairs,
        "background",
    )
    protons = exact_index(
        read_tsv(REDUCTION / "06_soft_protons/soft_protons.tsv"),
        lambda row: (row["frame"], row["band"]),
        pairs,
        "soft-proton",
    )
    ephemeris_rows = read_tsv(REDUCTION / "07_comove/ephemeris.tsv")
    names = ("mjd", "ra_deg", "dec_deg", "r_au", "delta_au")
    ephemeris = {
        name: np.array([float(row[name]) for row in ephemeris_rows]) for name in names
    }
    if (
        len(ephemeris_rows) < 2
        or np.any(~np.isfinite(np.column_stack([ephemeris[name] for name in names])))
        or np.any(np.diff(ephemeris["mjd"]) <= 0)
    ):
        raise RuntimeError("Task 01 ephemeris is invalid")
    frames = []
    for row in frame_rows:
        detector = frame_detector(row["frame"])
        masked, background, proton = (
            events[row["frame"]],
            backgrounds[row["frame"], "soft"],
            protons[row["frame"], "soft"],
        )
        if any(item["det"] != detector for item in (row, masked, background, proton)):
            raise RuntimeError(f"Task 01 detector mismatch for {row['frame']}")
        item = dict(row)
        item.update(
            event=Path(masked["event"]),
            oot=Path(masked["oot"]) if masked["oot"] else None,
            qpb_pha=Path(background["qpb_pha"]),
            qpb_map=Path(background["qpb_map"]),
            exposure_map=Path(background["exposure_map"]),
            chips=background["chips"],
            sp_pha=Path(proton["sp_pha"]),
            sp_map=Path(proton["sp_map"]),
        )
        item["position"] = extract.center(item, ephemeris)
        frames.append(item)
    source_rows = read_tsv(REDUCTION / "04_sources/sources.tsv")
    sources = [
        {
            "ra": float(row["ra"]),
            "dec": float(row["dec"]),
            "radius": float(row["radius_arcsec"]),
        }
        for row in source_rows
    ]
    identities = {(item["ra"], item["dec"]) for item in sources}
    if (
        len(sources) != 15
        or len(identities) != 15
        or not np.isfinite(
            [[item["ra"], item["dec"], item["radius"]] for item in sources]
        ).all()
        or any(
            item["radius"] != SHARED["source_mask_radius_arcsec"] for item in sources
        )
    ):
        raise RuntimeError("Task 01 source catalog is invalid")
    for i, left in enumerate(sources):
        for right in sources[i + 1 :]:
            if (
                extract.separation(left["ra"], left["dec"], right["ra"], right["dec"])
                <= left["radius"] + right["radius"]
            ):
                raise RuntimeError(
                    "overlapping source holes require union-area accounting"
                )
    return frames, sources


def reload_groups():
    rows = read_tsv(OUTPUT / "spectra.tsv")
    expected = {
        (region, detector) for region in REGIONS for detector in ("pn", "mos1", "mos2")
    }
    actual = [(row["region"], row["det"]) for row in rows]
    if len(actual) != len(expected) or set(actual) != expected:
        raise RuntimeError("Task 03 spectra table is incomplete or duplicated")
    groups = []
    for row in rows:
        stem = f"{row['det']}-{row['region']}"
        names = {
            "source_pha": f"{stem}-source.pi",
            "background_pha": f"{stem}-background.pi",
            "rmf": f"{stem}.rmf",
            "arf": f"{stem}.arf",
        }
        if any(row[name] != value for name, value in names.items()):
            raise RuntimeError("Task 03 spectra filenames are invalid")
        group = {
            name: float(row[name])
            for name in (
                "exposure_s",
                "omega_arcmin2",
                "full_area_arcmin2",
                "response_throughput",
                "source_counts",
                "background_counts",
                "fit_low_keV",
                "fit_high_keV",
                "pn_oot_fraction",
            )
        }
        group.update(
            region=row["region"], det=row["det"], n_frames=int(row["n_frames"])
        )
        group["group_min_counts"] = int(row["group_min_counts"])
        group["channel_max"] = int(row["channel_max"])
        inner, outer = REGIONS[group["region"]]
        expected = (
            *FIT_BAND_KEV,
            np.pi * (outer**2 - inner**2) / 3600,
            SHARED["pn_oot_fraction"],
        )
        actual = (
            group["fit_low_keV"],
            group["fit_high_keV"],
            group["full_area_arcmin2"],
            group["pn_oot_fraction"],
        )
        if (
            group["group_min_counts"] != PARAMS["group_min_counts"]
            or group["channel_max"] != PARAMS["channel_max"][group["det"]]
            or not np.allclose(actual, expected, rtol=0, atol=1e-12)
        ):
            raise RuntimeError("Task 03 extraction parameters changed; rerun extract")
        group.update(
            {
                name: OUTPUT / "stack" / row[name]
                for name in ("source_pha", "background_pha", "rmf", "arf")
            }
        )
        groups.append(group)
    return groups


def validate_extraction_frames(frames):
    expected = {
        (frame["frame"], frame["det"], region) for frame in frames for region in REGIONS
    }
    rows = read_tsv(OUTPUT / "extraction.tsv")
    actual = [(row["frame"], row["det"], row["region"]) for row in rows]
    if len(actual) != len(expected) or set(actual) != expected:
        raise RuntimeError("Task 03 extraction frames changed; rerun extract")


def require_not_newer(paths, stamp):
    paths = [Path(path) for path in paths if path]
    if any(not path.is_file() for path in paths):
        raise RuntimeError("Task 03 extraction input is missing; rerun extract")
    if any(path.stat().st_mtime_ns > stamp.stat().st_mtime_ns for path in paths):
        raise RuntimeError("Task 03 extraction input changed; rerun extract")


def validate_fresh_inputs(frames, groups):
    paths = [
        HERE / name
        for name in (
            "extract.py",
            "inputs.py",
            "spectrum.py",
            "run.py",
            "parameters.json",
        )
    ]
    paths += [
        ROOT / "parameters.json",
        ROOT / "analysis_01_reduction/settings.json",
        REDUCTION / "01_init/ccf.cif",
        REDUCTION / "01_init/odf" / Path(SOFTWARE["summary"]).name,
    ]
    paths += [
        REDUCTION / name
        for name in (
            "03_events/frames.tsv",
            "04_sources/events.tsv",
            "04_sources/sources.tsv",
            "05_background/background.tsv",
            "06_soft_protons/soft_protons.tsv",
            "07_comove/ephemeris.tsv",
        )
    ]
    paths += [
        frame[name]
        for frame in frames
        for name in (
            "event",
            "oot",
            "qpb_pha",
            "qpb_map",
            "exposure_map",
            "sp_pha",
            "sp_map",
        )
    ]
    paths += [
        group[name]
        for group in groups
        for name in ("source_pha", "background_pha", "rmf", "arf")
    ]
    require_not_newer(paths, OUTPUT / "spectra.tsv")
