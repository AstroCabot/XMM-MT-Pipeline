import json
import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
from scipy.signal import fftconvolve
from common import CFG, PARAMS, directory, finish, read_tsv, sky_wcs, write_tsv
from events import band_mask
from sources import event_wcs

STAGE = "07_comove"
PIXEL_ARCSEC, HALF_WIDTH_ARCMIN = 4.0, 25.0
XMM_MJDREF = 50814.0
PN_OOT_SCALE = PARAMS["pn_oot_fraction"]
KINDS = ("events", "oot", "exposure", "qpb", "sp")


def mjd(times, header):
    ref = float(
        header.get(
            "MJDREF", header.get("MJDREFI", XMM_MJDREF) + header.get("MJDREFF", 0)
        )
    )
    scale = str(header.get("TIMESYS", "TT")).strip().lower()
    offset = TimeDelta(times, format="sec")
    return np.asarray((Time(ref, format="mjd", scale=scale) + offset).utc.mjd)


def load_ephemeris(path, atthk):
    if not path.exists():
        from astroquery.jplhorizons import Horizons

        with fits.open(atthk, memmap=True) as hdul:
            table = next(
                h
                for h in hdul
                if h.data is not None
                and getattr(h.data, "names", None)
                and "TIME" in h.data.names
            )
            times = np.array([table.data["TIME"].min(), table.data["TIME"].max()])
            limits = mjd(times, table.header)
        pad = CFG["comove"]["cadence_s"] / 86400
        epochs = {
            "start": Time(limits[0] - pad, format="mjd").iso,
            "stop": Time(limits[1] + pad, format="mjd").iso,
            "step": f"{round(CFG['comove']['cadence_s'] / 60)}m",
        }
        query = Horizons(
            id=CFG["comove"]["target"], location="500@-125989", epochs=epochs
        )
        query = query.ephemerides(quantities="1,19,20", extra_precision=True)
        rows = [
            {
                "mjd": f"{float(r['datetime_jd']) - 2400000.5:.8f}",
                "ra_deg": f"{float(r['RA']):.7f}",
                "dec_deg": f"{float(r['DEC']):.7f}",
                "r_au": f"{float(r['r']):.6f}",
                "delta_au": f"{float(r['delta']):.6f}",
            }
            for r in query
        ]
        if not rows:
            raise RuntimeError("Horizons returned no XMM ephemeris")
        write_tsv(path, tuple(rows[0]), rows)
    rows = read_tsv(path)
    return {
        key: np.array([float(row[key]) for row in rows])
        for key in ("mjd", "ra_deg", "dec_deg", "r_au", "delta_au")
    }


def position(times, eph):
    times = np.asarray(times)
    if times.min() < eph["mjd"].min() or times.max() > eph["mjd"].max():
        raise RuntimeError("ephemeris does not cover an input GTI")
    ra = (
        np.rad2deg(np.interp(times, eph["mjd"], np.unwrap(np.deg2rad(eph["ra_deg"]))))
        % 360
    )
    return ra, np.interp(times, eph["mjd"], eph["dec_deg"])


def output_wcs(ra, dec):
    n = int(np.ceil(2 * HALF_WIDTH_ARCMIN * 60 / PIXEL_ARCSEC))
    n += (n + 1) % 2
    wcs = WCS(naxis=2)
    wcs.wcs.crval = [ra, dec]
    wcs.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    wcs.wcs.cdelt = [-PIXEL_ARCSEC / 3600, PIXEL_ARCSEC / 3600]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return wcs, (n, n)


def event_gti(path):
    with fits.open(path, memmap=True) as hdul:
        ontime = float(hdul["EVENTS"].header["ONTIME"])
        candidates = []
        for hdu in hdul:
            if (
                hdu.data is None
                or not getattr(hdu.data, "names", None)
                or not {"START", "STOP"} <= set(hdu.data.names)
            ):
                continue
            spans = np.column_stack((hdu.data["START"], hdu.data["STOP"])).astype(float)
            duration = float(np.sum(spans[:, 1] - spans[:, 0]))
            candidates.append((abs(duration - ontime), spans))
    if not candidates:
        raise RuntimeError(f"no event GTI: {path}")
    difference, spans = min(candidates, key=lambda value: value[0])
    if (
        not len(spans)
        or np.any(spans[:, 1] <= spans[:, 0])
        or np.any(spans[1:, 0] < spans[:-1, 1])
        or difference > 0.01
    ):
        raise RuntimeError(f"invalid event GTI: {path}")
    return spans, ontime


def samples(frame, event=None):
    spans, ontime = event_gti(event or frame["event"])
    if abs(ontime - float(frame["ontime"])) > 0.051:
        raise RuntimeError(f"ONTIME changed for {frame['frame']}")
    times, weights = [], []
    step = float(CFG["comove"]["map_step_s"])
    for start, stop in spans:
        count = max(1, int(np.ceil((stop - start) / step)))
        edges = np.linspace(start, stop, count + 1)
        times.extend((edges[:-1] + edges[1:]) / 2)
        weights.extend(np.diff(edges))
    return np.asarray(times), np.asarray(weights)


def deposit(shape, x, y, values):
    image = np.zeros(shape, float)
    ix, iy = np.floor(x).astype(int), np.floor(y).astype(int)
    fx, fy = x - ix, y - iy
    for dx, dy, weight in (
        (0, 0, (1 - fx) * (1 - fy)),
        (1, 0, fx * (1 - fy)),
        (0, 1, (1 - fx) * fy),
        (1, 1, fx * fy),
    ):
        xx, yy = ix + dx, iy + dy
        use = (xx >= 0) & (xx < shape[1]) & (yy >= 0) & (yy < shape[0])
        np.add.at(image, (yy[use], xx[use]), values[use] * weight[use])
    return image


def event_image(path, det, band, chips, out_wcs, shape, eph):
    with fits.open(path, memmap=True) as hdul:
        table = hdul["EVENTS"]
        data, header = table.data, table.header.copy()
        disabled = [i + 1 for i, flag in enumerate(chips.split()) if flag == "F"]
        use = band_mask(data, det, band, disabled)
        use &= (
            np.isfinite(data["X"]) & np.isfinite(data["Y"]) & np.isfinite(data["TIME"])
        )
        x, y, times = data["X"][use], data["Y"][use], data["TIME"][use]
    if not len(times):
        return np.zeros(shape)
    ra, dec = event_wcs(path).wcs_pix2world(x, y, 1)
    sx, sy = out_wcs.wcs_world2pix(ra, dec, 0)
    cra, cdec = position(mjd(times, header), eph)
    cx, cy = out_wcs.wcs_world2pix(cra, cdec, 0)
    center = (np.array(shape[::-1]) - 1) / 2
    x, y = sx - cx + center[0], sy - cy + center[1]
    image = np.histogram2d(
        y, x, bins=shape, range=((-0.5, shape[0] - 0.5), (-0.5, shape[1] - 0.5))
    )[0]
    if image.sum() != len(times):
        raise RuntimeError(f"event loss in comet reprojection: {path}")
    return image


def map_image(path, kind, times, weights, event_header, out_wcs, shape, eph):
    with fits.open(path, memmap=True) as hdul:
        hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
        data, source_wcs = np.asarray(hdu.data, float), sky_wcs(hdu.header)
    keep = np.isfinite(data) & (data != 0)
    y, x = np.nonzero(keep)
    ra, dec = source_wcs.pixel_to_world_values(x, y)
    sx, sy = out_wcs.world_to_pixel_values(ra, dec)
    tmjd = mjd(times, event_header)
    cra, cdec = position(tmjd, eph)
    cx, cy = out_wcs.world_to_pixel_values(cra, cdec)
    base_x, base_y = np.average(cx, weights=weights), np.average(cy, weights=weights)
    center = (np.array(shape[::-1]) - 1) / 2
    # Exposure is intensive; count maps are extensive under resampling.
    factor = 1
    if kind == "exposure":
        factor = abs(np.linalg.det(source_wcs.pixel_scale_matrix)) / abs(
            np.linalg.det(out_wcs.pixel_scale_matrix)
        )
    base = deposit(
        shape, sx - base_x + center[0], sy - base_y + center[1], data[keep] * factor
    )
    dx, dy = -(cx - base_x), -(cy - base_y)
    radius = int(np.ceil(max(np.max(abs(dx)), np.max(abs(dy))))) + 2
    kernel = deposit(
        (2 * radius + 1,) * 2, radius + dx, radius + dy, weights / weights.sum()
    )
    output = fftconvolve(base, kernel, mode="same")
    expected = float(np.nansum(data) * factor)
    if abs(output.sum() - expected) > 1e-5 * max(abs(expected), 1):
        raise RuntimeError(f"map integral changed in comet reprojection: {path}")
    return output


def write_image(path, data, wcs, det, band, kind):
    low, high = CFG["bands"][band]["pi"]
    header = wcs.to_header()
    header.update(
        BUNIT="s" if kind == "exposure" else "count",
        DETECTOR=det,
        BAND=band,
        E_MIN=low,
        E_MAX=high,
        COMOVING=True,
    )
    fits.PrimaryHDU(np.asarray(data, np.float32), header).writeto(path, overwrite=True)


def run():
    work, stack = directory(STAGE), directory(STAGE) / "stack"
    stack.mkdir(exist_ok=True)
    frames = read_tsv(directory("03_events") / "frames.tsv")
    for exposure in {row["exposure"] for row in frames if row["det"] == "pn"}:
        event = directory("02_reprocess") / "events" / f"{exposure}.fits"
        if fits.getheader(event).get("SUBMODE") != "PrimeFullWindow":
            raise RuntimeError(f"unsupported PN OoT mode: {exposure}")
    masked = {
        row["frame"]: row for row in read_tsv(directory("04_sources") / "events.tsv")
    }
    back = {
        (row["frame"], row["band"]): row
        for row in read_tsv(directory("05_background") / "background.tsv")
    }
    protons = {
        (row["frame"], row["band"]): row
        for row in read_tsv(directory("06_soft_protons") / "soft_protons.tsv")
    }
    eph = load_ephemeris(
        work / "ephemeris.tsv", directory("02_reprocess") / "atthk.fits"
    )
    ra0 = float(np.rad2deg(np.median(np.unwrap(np.deg2rad(eph["ra_deg"]))))) % 360
    dec0 = float(np.median(eph["dec_deg"]))
    out_wcs, shape = output_wcs(ra0, dec0)
    detectors = ("pn", "mos1", "mos2")
    stacks = {
        (det, band, kind): np.zeros(shape)
        for det in detectors
        for band in CFG["bands"]
        for kind in KINDS
        if det == "pn" or kind != "oot"
    }
    duration = {(det, band): 0.0 for det in detectors for band in CFG["bands"]}
    for frame in frames:
        det, fid = frame["det"], frame["frame"]
        times, weights = samples(frame, masked[fid]["event"])
        with fits.open(masked[fid]["event"], memmap=True) as hdul:
            event_header = hdul["EVENTS"].header.copy()
        for band in CFG["bands"]:
            background = back[fid, band]
            stacks[det, band, "events"] += event_image(
                masked[fid]["event"],
                det,
                band,
                background["chips"],
                out_wcs,
                shape,
                eph,
            )
            if masked[fid]["oot"]:
                stacks[det, band, "oot"] += event_image(
                    masked[fid]["oot"],
                    det,
                    band,
                    background["chips"],
                    out_wcs,
                    shape,
                    eph,
                )
            paths = {
                "exposure": background["exposure_map"],
                "qpb": background["qpb_map"],
                "sp": protons[fid, band]["sp_map"],
            }
            for kind, path in paths.items():
                stacks[det, band, kind] += map_image(
                    path, kind, times, weights, event_header, out_wcs, shape, eph
                )
            duration[det, band] += weights.sum()
    products = []
    pixel_area = (PIXEL_ARCSEC / 60) ** 2
    for det in detectors:
        for band in CFG["bands"]:
            for kind in KINDS:
                if det != "pn" and kind == "oot":
                    continue
                write_image(
                    stack / f"{det}-{kind}-{band}.fits",
                    stacks[det, band, kind],
                    out_wcs,
                    det,
                    band,
                    kind,
                )
            oot = stacks.get((det, band, "oot"))
            products.append(
                {
                    "det": det,
                    "band": band,
                    "gti_s": round(duration[det, band], 1),
                    "events": round(stacks[det, band, "events"].sum()),
                    "oot": round(oot.sum()) if oot is not None else 0,
                    "qpb_image_counts": round(stacks[det, band, "qpb"].sum(), 2),
                    "sp_counts": round(stacks[det, band, "sp"].sum(), 2),
                    "exposure_s_arcmin2": round(
                        stacks[det, band, "exposure"].sum() * pixel_area, 1
                    ),
                }
            )
    result = {
        "observer": "XMM (500@-125989)",
        "pixel_arcsec": PIXEL_ARCSEC,
        "pn_oot_scale": PN_OOT_SCALE,
        "reference_ra_deg": ra0,
        "reference_dec_deg": dec0,
        "ephemeris": "ephemeris.tsv",
        "products": products,
    }
    (work / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    finish(STAGE)
