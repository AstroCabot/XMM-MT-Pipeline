import re
import shutil
from pathlib import Path

import numpy as np
from astropy.io import fits

from common import CFG, directory, sas, valid_fits, write_tsv

EVENT_NAME = re.compile(
    r"^\d+_\d+_(EMOS1|EMOS2|EPN)_([SU]\d+)_(?:OoT)?ImagingEvts\.ds$"
)
DETECTOR = {"EMOS1": "mos1", "EMOS2": "mos2", "EPN": "pn"}


def event_info(path, oot=False):
    match = EVENT_NAME.match(path.name)
    if not match:
        return None
    with fits.open(path) as hdul:
        header = hdul["EVENTS"].header
        n_events = int(header.get("NAXIS2", 0))
        start = float(header.get("TSTART", 0))
        stop = float(header.get("TSTOP", 0))
        filt = str(header.get("FILTER", "Unknown")).strip()
    det, expid = DETECTOR[match[1]], match[2].lower()
    return {
        "det": det,
        "expid": expid,
        "oot": oot,
        "name": f"{det}-{expid}{'-oot' if oot else ''}.fits",
        "filter": filt,
        "events": n_events,
        "start": start,
        "stop": stop,
        "source": path,
    }


def run_task(name, command):
    root = directory("02_reprocess") / name
    marker = root / ".done"
    if marker.exists():
        return root
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    sas(command[0], root, **command[1])
    marker.touch()
    return root


def real_gti(gti, event):
    with fits.open(gti, memmap=True) as hdul:
        data = hdul[1].data
        start, stop = np.asarray(data["START"]), np.asarray(data["STOP"])
    with fits.open(event, memmap=True) as hdul:
        header = hdul["EVENTS"].header
        first, last = float(header["TSTART"]), float(header["TSTOP"])
    overlap = np.maximum(np.minimum(stop, last) - np.maximum(start, first), 0)
    sentinel = len(start) == 1 and start[0] <= first and stop[0] >= last * 1000
    return bool(overlap.sum() > 0 and not sentinel)


def run_reprocessing():
    mos_ids = "M1S001 M1U009 M2S002 M2U009"
    pn_ids = "PNS003"
    mos = run_task(
        "emproc", ("emproc", {"withinstexpids": True, "instexpids": mos_ids})
    )
    pn = run_task("epproc", ("epproc", {"withinstexpids": True, "instexpids": pn_ids}))
    oot = run_task(
        "epproc_oot",
        (
            "epproc",
            {"withoutoftime": True, "withinstexpids": True, "instexpids": pn_ids},
        ),
    )
    found = []
    for root, is_oot in ((mos, False), (pn, False), (oot, True)):
        for path in root.glob("*ImagingEvts.ds"):
            row = event_info(path, is_oot)
            if row:
                found.append(row)

    science = {(row["det"], row["expid"]): row for row in found if not row["oot"]}
    selected = []
    output = directory("02_reprocess") / "events"
    output.mkdir(exist_ok=True)
    for row in sorted(found, key=lambda x: (x["det"], x["expid"], x["oot"])):
        parent = science.get((row["det"], row["expid"]))
        keep = bool(
            parent
            and parent["filter"] in {"Thin1", "Thin2", "Medium", "Thick"}
            and parent["events"] > 0
            and parent["stop"] - parent["start"] >= 100
        )
        if keep:
            destination = output / row["name"]
            if not valid_fits(destination):
                partial = destination.with_suffix(".partial")
                shutil.copyfile(row["source"], partial)
                if not valid_fits(partial):
                    raise RuntimeError(f"invalid copied events: {partial}")
                partial.replace(destination)
            selected.append({**row, "path": destination})

    atthk = directory("02_reprocess") / "atthk.fits"
    if not atthk.exists():
        sas("atthkgen", atthk.parent, atthkset=atthk)
    for row in selected:
        if row["oot"] or row["det"] == "pn":
            continue
        path = row["path"]
        if "ANOMFL2" not in fits.getheader(path):
            sas(
                "emanom",
                path.parent,
                eventfile=path,
                cornerfile=path.parent / f"{path.stem}-corners.fits",
                keepcorner=False,
                writekeys=True,
                writelog=False,
            )

    columns = ("det", "expid", "oot", "filter", "events", "start", "stop", "path")
    rows = [{key: row[key] for key in columns} for row in selected]
    write_tsv(directory("02_reprocess") / "exposures.tsv", columns, rows)
    return rows, atthk


def filter_exposures(rows):
    science = [row for row in rows if not row["oot"]]
    oot = {(row["det"], row["expid"]): row for row in rows if row["oot"]}
    cfg = CFG["espfilt"]
    result = []
    for row in science:
        base = Path(row["path"]).stem
        work = directory("02_reprocess") / "filter" / base
        if work.exists() and not (work / ".done").exists():
            shutil.rmtree(work)
        work.mkdir(parents=True, exist_ok=True)
        event = work / Path(row["path"]).name
        if not event.exists():
            event.symlink_to(Path(row["path"]))
        params = {"eventfile": event.name}
        companion = oot.get((row["det"], row["expid"]))
        if companion:
            oot_link = work / Path(companion["path"]).name
            if not oot_link.exists():
                oot_link.symlink_to(Path(companion["path"]))
            params.update(withoot="Y", ootfile=oot_link.name)
        params.update(
            method="histogram",
            withsmoothing=True,
            smooth=cfg["smooth"],
            elow=cfg["elow"],
            ehigh=cfg["ehigh"],
            rangescale=(
                cfg["rangescale_pn"] if row["det"] == "pn" else cfg["rangescale_mos"]
            ),
            allowsigma=cfg["allowsigma"],
            keepinterfiles=False,
        )
        if not (work / ".done").exists():
            sas("espfilt", work, **params)
            (work / ".done").touch()
        gti = next(iter(sorted(work.glob("*-gti.fits"))), None)
        allevc = next(
            (
                p
                for p in sorted(work.glob("*-allevc.fits"))
                if "allevcoot" not in p.name
            ),
            None,
        )
        if gti and allevc and real_gti(gti, Path(row["path"])):
            item = {**row, "gti": gti, "event": allevc}
            item["oot_event"] = next(iter(sorted(work.glob("*-allevcoot.fits"))), None)
            item["core"] = next(iter(sorted(work.glob("*-corevc.fits"))), None)
            result.append(item)
    return result
