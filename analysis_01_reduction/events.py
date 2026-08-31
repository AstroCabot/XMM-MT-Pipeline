import shutil
from pathlib import Path

import numpy as np
from astropy.io import fits

from common import CFG, directory, finish, sas, valid_fits, write_tsv
from reprocess import filter_exposures, run_reprocessing

STAGE = "03_events"
FLAG_EXPRESSION = "FLAG==0"


def pointing_attitude(gti, atthk):
    with fits.open(gti) as hdul:
        data = hdul[1].data
        start, stop = float(data["START"].min()), float(data["STOP"].max())
    with fits.open(atthk) as hdul:
        data = hdul["ATTHK"].data
        use = (data["TIME"] >= start) & (data["TIME"] <= stop)
        if not use.any():
            use[np.argmin(abs(data["TIME"] - (start + stop) / 2))] = True
        values = tuple(
            float(np.median(data[name][use])) for name in ("AHFRA", "AHFDEC", "AHFPA")
        )
    return start, stop, values


def rewrite_pointing(path, start, stop, attitude):
    with fits.open(path, memmap=False) as source:
        hdul = fits.HDUList([hdu.copy() for hdu in source])
    for hdu in (hdul[0], hdul["EVENTS"]):
        hdu.header["TSTART"] = start
        hdu.header["TSTOP"] = stop
    for key, value in zip(("RA_PNT", "DEC_PNT", "PA_PNT"), attitude):
        hdul[0].header[key] = value
    temporary = path.with_suffix(path.suffix + ".new")
    hdul.writeto(temporary, overwrite=True)
    temporary.replace(path)


def anomaly_drops(canonical, band):
    if canonical.name.startswith("pn"):
        return set()
    header = fits.getheader(canonical)
    states = {ccd: header.get(f"ANOMFL{ccd}", "U") for ccd in range(2, 8)}
    bad = {"O"} if band == "hard" else {"O", "B", "I"}
    return {ccd for ccd, state in states.items() if state in bad}


def band_mask(data, det, band, disabled=()):
    config = CFG["bands"][band]
    lo, hi = config["pi"]
    selection = config["pn" if det == "pn" else "mos"]
    good = (
        (data["PI"] >= lo)
        & (data["PI"] <= hi)
        & (data["PATTERN"] <= selection["pattern"])
        & (data["FLAG"] == 0)
    )
    for ccd in disabled:
        good &= data["CCDNR"] != ccd
    return good


def split_exposure(row, atthk, canonical):
    base = canonical.stem
    root = directory(STAGE) / "work" / base
    root.mkdir(parents=True, exist_ok=True)
    event = root / Path(row["event"]).name
    if not event.exists():
        event.symlink_to(Path(row["event"]))
    params = {
        "atthkfile": atthk,
        {"mos1": "mos1evtfile", "mos2": "mos2evtfile", "pn": "pnevtfile"}[
            row["det"]
        ]: event.name,
    }
    if not (root / ".prepared").exists():
        for partial in root.glob("prep_mosaic_*"):
            shutil.rmtree(partial)
        sas("emosaic_prep", root, **params)
        (root / ".prepared").touch()
    output = []
    for prep in sorted(root.glob("prep_mosaic_*")):
        pid = int(prep.name.rsplit("_", 1)[1])
        frame = f"{base}-p{pid:02d}"
        gti = prep / "gti_pointing_position.ds"
        start, stop, attitude = pointing_attitude(gti, atthk)
        frame_dir = directory(STAGE) / "unmasked" / frame
        frame_dir.mkdir(parents=True, exist_ok=True)
        science = frame_dir / f"{frame}.fits"
        oot_out = ""
        oot_in = None
        if row["det"] == "pn" and row.get("oot_event"):
            oot_in = root / Path(row["oot_event"]).name
            if not oot_in.exists():
                oot_in.symlink_to(Path(row["oot_event"]))
            oot_out = frame_dir / f"{frame}-oot.fits"
        products = [science] + ([oot_out] if oot_out else [])
        marker = frame_dir / ".done"
        ready = marker.exists() and all(valid_fits(path) for path in products)
        if not ready:
            marker.unlink(missing_ok=True)
            for path in products:
                path.unlink(missing_ok=True)
            sas(
                "evselect",
                prep,
                table=f"{event}:EVENTS",
                withfilteredset=True,
                filteredset=science,
                destruct=True,
                keepfilteroutput=True,
                updateexposure=True,
                filterexposure=True,
                expression=f"gti({gti},TIME)",
            )
            rewrite_pointing(science, start, stop, attitude)
            sas(
                "attcalc",
                frame_dir,
                eventset=science,
                attitudelabel="ahf",
                withatthkset=True,
                atthkset=atthk,
                refpointlabel="user",
                nominalra=f"{attitude[0]:.6f}",
                nominaldec=f"{attitude[1]:.6f}",
                setpnttouser=True,
            )
            # attcalc changes PA_PNT; retain the per-pointing AHF attitude.
            rewrite_pointing(science, start, stop, attitude)
            if oot_out:
                sas(
                    "evselect",
                    prep,
                    table=f"{oot_in}:EVENTS",
                    withfilteredset=True,
                    filteredset=oot_out,
                    destruct=True,
                    keepfilteroutput=True,
                    updateexposure=True,
                    filterexposure=True,
                    expression=f"gti({gti},TIME)",
                )
                rewrite_pointing(oot_out, start, stop, attitude)
                sas(
                    "attcalc",
                    frame_dir,
                    eventset=oot_out,
                    attitudelabel="ahf",
                    withatthkset=True,
                    atthkset=atthk,
                    refpointlabel="user",
                    nominalra=f"{attitude[0]:.6f}",
                    nominaldec=f"{attitude[1]:.6f}",
                    setpnttouser=True,
                )
                rewrite_pointing(oot_out, start, stop, attitude)
            if not all(valid_fits(path) for path in products):
                raise RuntimeError(f"invalid frame products for {frame}")
            marker.touch()
        with fits.open(science, memmap=True) as hdul:
            data = hdul["EVENTS"].data
            counts = {
                band: int(
                    band_mask(
                        data, row["det"], band, anomaly_drops(canonical, band)
                    ).sum()
                )
                for band in CFG["bands"]
            }
            ontime = float(hdul["EVENTS"].header.get("ONTIME", 0))
        output.append(
            {
                "frame": frame,
                "det": row["det"],
                "exposure": base,
                "pointing": pid,
                "event": science,
                "oot": oot_out,
                "filter_gti": row["gti"],
                "pointing_gti": gti,
                "start": f"{start:.3f}",
                "stop": f"{stop:.3f}",
                "ontime": f"{ontime:.1f}",
                "soft": counts["soft"],
                "hard": counts["hard"],
                "ra": f"{attitude[0]:.6f}",
                "dec": f"{attitude[1]:.6f}",
                "pa": f"{attitude[2]:.6f}",
            }
        )
    return output


def selection_rows(rows, kept):
    output = []
    for exposure in sorted({row["exposure"] for row in rows}):
        group = [row for row in rows if row["exposure"] == exposure]
        accepted = [row for row in group if row["frame"] in kept]
        max_time = max(float(row["ontime"]) for row in group)
        max_soft = max(int(row["soft"]) for row in group)
        if accepted:
            decision = "kept"
        elif max_time < CFG["frame_cut"]["ontime_s"]:
            decision = (
                f"real espfilt GTI; maximum frame ONTIME {max_time:.0f} s < "
                f"{CFG['frame_cut']['ontime_s']:.0f} s"
            )
        elif max_soft < CFG["frame_cut"]["soft_events"]:
            decision = f"maximum frame soft events {max_soft} < {CFG['frame_cut']['soft_events']}"
        else:
            decision = "no frame met both cuts"
        output.append(
            {
                "exposure": exposure,
                "det": group[0]["det"],
                "pointings": ",".join(str(row["pointing"]) for row in accepted),
                "decision": decision,
            }
        )
    output += [{**row, "pointings": ""} for row in CFG["excluded_exposures"]]
    return sorted(output, key=lambda row: (row["det"], row["exposure"]))


def run():
    exposures, atthk = run_reprocessing()
    filtered = filter_exposures(exposures)
    canonical = {
        path.stem: path
        for path in (directory("02_reprocess") / "events").glob("*.fits")
        if not path.stem.endswith("-oot")
    }
    rows = []
    for exposure in filtered:
        rows += split_exposure(exposure, atthk, canonical[Path(exposure["path"]).stem])
    cut = CFG["frame_cut"]
    kept = [
        row
        for row in rows
        if int(row["soft"]) >= cut["soft_events"]
        and float(row["ontime"]) >= cut["ontime_s"]
    ]
    names = {row["frame"] for row in kept}
    if names != set(CFG["expected_frames"]):
        raise RuntimeError(f"unexpected retained frames: {sorted(names)}")
    selection = selection_rows(rows, names)
    for row in rows:
        if row["frame"] not in names:
            shutil.rmtree(Path(row["event"]).parent)
    used = {Path(row["pointing_gti"]).parent.resolve() for row in kept}
    for root in (directory(STAGE) / "work").iterdir():
        preps = list(root.glob("prep_mosaic_*"))
        root_used = any(prep.resolve() in used for prep in preps)
        for prep in preps:
            if prep.resolve() not in used:
                shutil.rmtree(prep)
        if not root_used:
            shutil.rmtree(root)
    rows = kept
    write_tsv(
        directory(STAGE) / "selection.tsv",
        ("exposure", "det", "pointings", "decision"),
        selection,
    )
    columns = tuple(rows[0])
    write_tsv(directory(STAGE) / "frames.tsv", columns, rows)
    finish(STAGE)
