import os
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs.utils import proj_plane_pixel_scales

import pha
from common import CFG, directory, finish, read_tsv, sas, sky_wcs, write_tsv
from events import anomaly_drops
from sources import FOV, PN_BOX

STAGE = "05_background"


def mask_sky_map(path, sources):
    with fits.open(path, memmap=False) as hdul:
        output = fits.HDUList([hdu.copy() for hdu in hdul])
    data, header = output[0].data, output[0].header
    wcs = sky_wcs(header)
    scale = float(proj_plane_pixel_scales(wcs)[0]) * 3600
    y, x = np.indices(data.shape)
    for source in sources:
        sx, sy = wcs.world_to_pixel_values(float(source["ra"]), float(source["dec"]))
        radius = float(source["radius_arcsec"]) / scale
        data[(x - sx) ** 2 + (y - sy) ** 2 <= radius**2] = 0
    output.writeto(path, overwrite=True)


def chip_string(frame, band, canonical):
    if frame["det"] == "pn":
        return "T T T T"
    disabled = anomaly_drops(canonical, band)
    if frame["exposure"] == "mos1-u009" and band == "hard":
        root = directory("02_reprocess") / "filter" / frame["exposure"]
        core = next(iter(sorted(root.glob("*-corevc.fits"))), None)
        if core is None:
            raise RuntimeError("missing MOS1 U009 corner events")
        with fits.open(core, memmap=True) as hdul:
            data = hdul["EVENTS"].data
            use = (
                (data["CCDNR"] == 5)
                & (data["TIME"] >= float(frame["start"]))
                & (data["TIME"] <= float(frame["stop"]))
            )
            low = np.count_nonzero(use & (data["PI"] >= 500) & (data["PI"] <= 800))
            high = np.count_nonzero(use & (data["PI"] >= 2500) & (data["PI"] <= 5000))
        if low and high:
            raise RuntimeError("MOS1 U009 CCD5 now has usable corner data")
        disabled.add(5)
    return " ".join("F" if ccd in disabled else "T" for ccd in range(1, 8))


def inputs_stamp(frame, event, band, chips, lo, hi, pattern):
    region = Path(event["region"]).read_text().strip()
    values = (
        band,
        f"{lo}-{hi}",
        str(pattern),
        chips,
        frame["soft"],
        frame["hard"],
        region,
    )
    if event.get("oot"):
        values += (event["oot"],)
    return "\n".join(values) + "\n"


def make_oot_detector(work, oot, region, output, pattern, lo, hi):
    box = ",".join(map(str, PN_BOX))
    fov = ",".join(map(str, FOV["EPN"]))
    ccds = "||".join(f"(CCDNR == {ccd})" for ccd in range(1, 13))
    region = region.read_text().strip()
    expression = (
        f"(PATTERN<={pattern})&&(FLAG==0)&&(PI in [{lo}:{hi}])"
        f"&&({ccds})&&((DETX,DETY) in BOX({box},0)){region}"
        f"&&((DETX,DETY) IN circle({fov}))"
    )
    output.unlink(missing_ok=True)
    sas(
        "evselect",
        work,
        table=f"{oot.name}:EVENTS",
        expression=expression,
        withimageset=True,
        imageset=output.name,
        xcolumn="DETX",
        ycolumn="DETY",
        squarepixels=True,
        imagebinning="imageSize",
        ximagesize=780,
        yimagesize=780,
        withxranges=True,
        ximagemin=-19499,
        ximagemax=19500,
        withyranges=True,
        yimagemin=-19499,
        yimagemax=19500,
        withimagedatatype=True,
        imagedatatype="Int32",
        ignorelegallimits=True,
    )


def one_job(frame, event, band, sources):
    work = directory(STAGE) / frame["frame"] / band
    canonical = directory("02_reprocess") / "events" / f"{frame['exposure']}.fits"
    chips = chip_string(frame, band, canonical)
    lo, hi = CFG["bands"][band]["pi"]
    selection = CFG["bands"][band]["pn" if frame["det"] == "pn" else "mos"]
    pattern = selection["pattern"]
    stamp = inputs_stamp(frame, event, band, chips, lo, hi, pattern)
    marker = work / ".done"
    if marker.exists() and marker.read_text() == stamp:
        return result_row(frame, band, work, chips)
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    # ESAS applies the DET exclusion to FOV, response, and FWC products.
    local_event = work / Path(frame["event"]).name
    local_event.symlink_to(Path(frame["event"]))
    local_oot = None
    if frame["oot"]:
        local_oot = work / Path(frame["oot"]).name
        local_oot.symlink_to(Path(frame["oot"]))
    region = work / "sources.txt"
    region.symlink_to(Path(event["region"]))
    atthk = directory("02_reprocess") / "atthk.fits"
    link = work / "atthk.fits"
    link.symlink_to(atthk)
    if frame["det"] == "pn":
        params = {
            "eventfile": local_event.name,
            "quads": chips,
            "pattern": pattern,
            "elow": lo,
            "ehigh": hi,
            "withregion": True,
            "regionfile": region.name,
        }
        if local_oot:
            params["ootevtfile"] = local_oot.name
        sas("pnspectra", work, **params)
    else:
        sas(
            "mosspectra",
            work,
            eventfile=local_event.name,
            ccds=chips,
            pattern=pattern,
            elow=lo,
            ehigh=hi,
            withregion=True,
            regionfile=region.name,
        )
    fov = next(iter(sorted(work.glob("*-fovt.pi"))))
    prefix = fov.name.removesuffix("-fovt.pi")
    if frame["det"] == "pn":
        make_oot_detector(
            work,
            local_oot,
            region,
            work / f"{prefix}-fovimootdet-{lo}-{hi}.fits",
            pattern,
            lo,
            hi,
        )
    if frame["det"] == "pn":
        params = {
            "inspecfile": fov.name,
            "rmffile": f"{prefix}.rmf",
            "quads": chips,
            "elow": lo,
            "ehigh": hi,
            "withplotfiles": False,
        }
        oot = work / f"{prefix}-fovtoot.pi"
        if oot.exists():
            params["inspecoot"] = oot.name
        sas("pnback", work, **params)
    else:
        sas(
            "mosback",
            work,
            inspecfile=fov.name,
            rmffile=f"{prefix}.rmf",
            ccds=chips,
            elow=lo,
            ehigh=hi,
            withplotfiles=False,
        )
    tag = f"{lo}-{hi}"
    qpb = work / f"{prefix}-qpb-sky-{tag}.fits"
    detector_qpb = work / f"{prefix}-bkgimdet-{tag}.fits"
    sas(
        "rotdet2sky",
        work,
        intemplate=work / f"{prefix}-expimsky-{tag}.fits",
        inimage=detector_qpb,
        outimage=qpb,
    )
    before = float(np.nansum(fits.getdata(detector_qpb)))
    after = float(np.nansum(fits.getdata(qpb)))
    if not np.isclose(after, before, rtol=0.01, atol=0.05):
        raise RuntimeError(f"rotdet2sky changed the QPB integral in {work}")
    sky_maps = [
        work / f"{prefix}-fovimsky-{tag}.fits",
        work / f"{prefix}-expimsky-{tag}.fits",
        qpb,
    ]
    oot_sky = work / f"{prefix}-fovimootsky-{tag}.fits"
    if oot_sky.exists():
        sky_maps.append(oot_sky)
    for path in sky_maps:
        mask_sky_map(path, sources)
    result = result_row(frame, band, work, chips)
    marker.write_text(stamp)
    return result


def result_row(frame, band, work, chips):
    fovs = sorted(work.glob("*-fovt.pi"))
    if len(fovs) != 1:
        raise RuntimeError(f"expected one source PHA in {work}, found {len(fovs)}")
    fov = fovs[0]
    prefix = fov.name.removesuffix("-fovt.pi")
    lo, hi = CFG["bands"][band]["pi"]
    tag = f"{lo}-{hi}"
    paths = {
        "qpb_pha": work / f"{prefix}-bkg.pi",
        "rmf": work / f"{prefix}.rmf",
        "arf": work / f"{prefix}.arf",
        "events_image": work / f"{prefix}-fovimsky-{tag}.fits",
        "exposure_map": work / f"{prefix}-expimsky-{tag}.fits",
        "qpb_map": work / f"{prefix}-qpb-sky-{tag}.fits",
    }
    for path in (fov, *paths.values()):
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError(f"missing background product: {path}")
        with fits.open(path):
            pass
    qpb_counts = pha.band_counts(pha.read(paths["qpb_pha"]), lo, hi)
    return {
        "frame": frame["frame"],
        "det": frame["det"],
        "band": band,
        "chips": chips,
        "dir": work,
        "source_pha": fov,
        **paths,
        "source_counts": f"{pha.band_counts(pha.read(fov), lo, hi):.1f}",
        "qpb_pha_counts": f"{qpb_counts:.1f}",
    }


def run():
    frames = {
        row["frame"]: row for row in read_tsv(directory("03_events") / "frames.tsv")
    }
    events = read_tsv(directory("04_sources") / "events.tsv")
    sources = read_tsv(directory("04_sources") / "sources.tsv")
    jobs = [
        (frames[row["frame"]], row, band, sources)
        for row in events
        for band in CFG["bands"]
    ]
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    with ThreadPoolExecutor(max_workers=2) as mos, ThreadPoolExecutor(
        max_workers=1
    ) as pn:
        futures = [
            (pn if job[0]["det"] == "pn" else mos).submit(one_job, *job) for job in jobs
        ]
        rows = [future.result() for future in futures]
    write_tsv(directory(STAGE) / "background.tsv", tuple(rows[0]), rows)
    finish(STAGE)
