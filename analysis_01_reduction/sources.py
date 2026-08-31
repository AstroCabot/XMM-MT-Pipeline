import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from common import CFG, directory, finish, read_tsv, sas, valid_fits, write_tsv
from events import FLAG_EXPRESSION

STAGE = "04_sources"
PHYSICAL_PIXEL_ARCSEC = 0.05
CATALOG = json.loads((Path(__file__).with_name("source_catalog.json")).read_text())
FOV = {
    "EPN": (-2200, -1110, 17980),
    "EMOS1": (97, -172, 17152),
    "EMOS2": (-61, -228, 17085),
}
PN_BOX = (-2196, -1110, 16060, 15510)


def event_wcs(path):
    header = fits.getheader(path, "EVENTS")
    wcs = WCS(naxis=2)
    wcs.wcs.crval = [header["REFXCRVL"], header["REFYCRVL"]]
    wcs.wcs.crpix = [header["REFXCRPX"], header["REFYCRPX"]]
    wcs.wcs.cdelt = [header["REFXCDLT"], header["REFYCDLT"]]
    wcs.wcs.ctype = [header["REFXCTYP"], header["REFYCTYP"]]
    return wcs


def intersects_fov(instrument, x, y, radius):
    fov_x, fov_y, fov_radius = FOV[instrument]
    if np.hypot(x - fov_x, y - fov_y) > fov_radius + radius:
        return False
    if instrument == "EPN":
        box_x, box_y, half_x, half_y = PN_BOX
        distance = np.maximum(np.abs([x - box_x, y - box_y]) - [half_x, half_y], 0)
        return np.hypot(*distance) <= radius
    return True


def apertures(path, sources):
    wcs = event_wcs(path)
    for source in sources:
        x, y = wcs.wcs_world2pix(float(source["ra"]), float(source["dec"]), 1)
        radius = float(source["radius_arcsec"]) / PHYSICAL_PIXEL_ARCSEC
        yield x, y, radius


def source_expression(path, sources):
    terms = [
        f"!((X,Y) IN circle({x:.3f},{y:.3f},{radius:.3f}))"
        for x, y, radius in apertures(path, sources)
    ]
    return "&&".join(terms) if terms else "true"


def detector_region(path, sources):
    with fits.open(path, memmap=True) as hdul:
        data = hdul["EVENTS"].data
        use = slice(None, None, max(1, len(data) // 10000))
        xy = np.column_stack((data["X"][use], data["Y"][use])).astype(float)
        det = np.column_stack((data["DETX"][use], data["DETY"][use])).astype(float)
        instrument = hdul["EVENTS"].header["INSTRUME"]
    center, det_center = xy.mean(axis=0), det.mean(axis=0)
    transform = np.linalg.lstsq(xy - center, det - det_center, rcond=None)[0]
    # SKY and DET coordinates share 0.05-arcsec units and differ by a rotation.
    residual = det - det_center - (xy - center) @ transform
    if np.sqrt(np.mean(residual**2)) > 20:
        raise RuntimeError(f"unstable SKY-to-DET transform in {path}")
    scale = np.linalg.svd(transform, compute_uv=False)
    if np.any(abs(scale - 1) > 0.02):
        raise RuntimeError(f"nonphysical SKY-to-DET scale in {path}: {scale}")
    terms = []
    for x, y, radius in apertures(path, sources):
        dx, dy = det_center + (np.array([x, y]) - center) @ transform
        if not intersects_fov(instrument, dx, dy, radius):
            continue
        terms.append(f"&&!((DETX,DETY) IN circle({dx:.1f},{dy:.1f},{radius:.1f}))")
    # ESAS appends region-file text to its own event expression.
    return "".join(terms)


def assert_excluded(path, sources):
    with fits.open(path, memmap=True) as hdul:
        data = hdul["EVENTS"].data
        x, y = np.asarray(data["X"]), np.asarray(data["Y"])
    for sx, sy, radius in apertures(path, sources):
        if np.any((x - sx) ** 2 + (y - sy) ** 2 <= radius**2):
            raise RuntimeError(f"source events remain in {path}")


def run():
    frames = read_tsv(directory("03_events") / "frames.tsv")
    rows = [
        {
            "ra": f"{ra:.6f}",
            "dec": f"{dec:.6f}",
            "peak_band": band,
            "radius_arcsec": f"{CFG['sources']['radius_arcsec']:.1f}",
        }
        for ra, dec, band in CATALOG
    ]
    columns = ("ra", "dec", "peak_band", "radius_arcsec")
    write_tsv(directory(STAGE) / "sources.tsv", columns, rows)
    region = ["# Region file format: DS9 version 4.1", "fk5"]
    region += [f"circle({r['ra']},{r['dec']},{r['radius_arcsec']}\")" for r in rows]
    (directory(STAGE) / "sources.reg").write_text("\n".join(region) + "\n")

    output = directory(STAGE) / "events"
    output.mkdir(exist_ok=True)
    regions = directory(STAGE) / "regions"
    regions.mkdir(exist_ok=True)
    masked = []
    for row in frames:
        region = regions / f"{row['frame']}.txt"
        region.write_text(detector_region(Path(row["event"]), rows) + "\n")
        item = {"frame": row["frame"], "det": row["det"], "region": region}
        for key in ("event", "oot"):
            if not row[key]:
                item[key] = ""
                continue
            source = Path(row[key])
            destination = (
                output / f"{row['frame']}{'-oot' if key == 'oot' else ''}.fits"
            )
            if not valid_fits(destination):
                temporary = destination.with_name(destination.stem + "-partial.fits")
                temporary.unlink(missing_ok=True)
                sas(
                    "evselect",
                    output,
                    table=f"{source}:EVENTS",
                    withfilteredset=True,
                    filteredset=temporary,
                    destruct=True,
                    keepfilteroutput=True,
                    updateexposure=True,
                    filterexposure=True,
                    expression=(
                        f"({FLAG_EXPRESSION})&&" f"({source_expression(source, rows)})"
                    ),
                )
                if not valid_fits(temporary):
                    raise RuntimeError(f"invalid masked events: {temporary}")
                temporary.replace(destination)
            assert_excluded(destination, rows)
            item[key] = destination
        masked.append(item)
    write_tsv(directory(STAGE) / "events.tsv", tuple(masked[0]), masked)
    finish(STAGE)
