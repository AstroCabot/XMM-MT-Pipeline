import csv
import math

import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta

EVENT_CENTER = (24904.3189645, 23850.7229654)
MAP_CENTER = (310.7977371, 297.6277871)
EVENT_PIXEL = 0.05
MAP_PIXEL = 4.0


def band_mask(energy, low, high, include_high=False):
    energy = np.asarray(energy)
    upper = energy <= high if include_high else energy < high
    return (energy >= low) & upper


def read_gtis(path):
    with fits.open(path, memmap=True) as hdul:
        return [
            np.column_stack(
                (
                    hdul[f"GTI{i:03d}07"].data["START"],
                    hdul[f"GTI{i:03d}07"].data["STOP"],
                )
            )
            for i in range(12)
        ]


def gti_expression():
    terms = [f"((CCDNR=={i + 1})&&gti(gti.FTZ:GTI{i:03d}07,TIME))" for i in range(12)]
    return "(" + "||".join(terms) + ")"


def time_samples(gtis, step=60.0):
    times, weights = [], []
    for start, stop in np.concatenate(gtis):
        n = max(1, int(np.ceil((stop - start) / step)))
        edges = np.linspace(start, stop, n + 1)
        times.extend((edges[:-1] + edges[1:]) / 2)
        weights.extend(np.diff(edges))
    return np.asarray(times), np.asarray(weights)


def event_mjd(times, header):
    ref = float(header.get("MJDREF", 50814.0))
    scale = str(header.get("TIMESYS", "TT")).strip().lower()
    return np.asarray(
        (Time(ref, format="mjd", scale=scale) + TimeDelta(times, format="sec")).utc.mjd
    )


def read_table(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def comet_position(mjd, ephemeris):
    time = ephemeris[:, 0]
    if np.min(mjd) < time[0] or np.max(mjd) > time[-1]:
        raise ValueError("ephemeris does not cover the PPS GTIs")
    ra = np.rad2deg(np.interp(mjd, time, np.unwrap(np.deg2rad(ephemeris[:, 1])))) % 360
    return ra, np.interp(mjd, time, ephemeris[:, 2])


def source_offsets(source, comet_ra, comet_dec):
    dra = (source[0] - comet_ra + 180) % 360 - 180
    east = dra * np.cos(np.deg2rad((source[1] + comet_dec) / 2)) * 3600
    return east, (source[1] - comet_dec) * 3600


def source_keep(x, y, mjd, sources, ephemeris, radius):
    east = -(np.asarray(x) - EVENT_CENTER[0]) * EVENT_PIXEL
    north = (np.asarray(y) - EVENT_CENTER[1]) * EVENT_PIXEL
    comet_ra, comet_dec = comet_position(mjd, ephemeris)
    keep = np.ones(len(east), bool)
    for source in sources:
        sx, sy = source_offsets(source, comet_ra, comet_dec)
        keep &= (east - sx) ** 2 + (north - sy) ** 2 > radius**2
    return keep


def exposure_coverage(gtis, event_header, sources, ephemeris, shape, radius_arcsec):
    times, weights = time_samples(gtis)
    comet_ra, comet_dec = comet_position(event_mjd(times, event_header), ephemeris)
    lost = np.zeros(shape, float)
    radius = radius_arcsec / MAP_PIXEL
    for source in sources:
        east, north = source_offsets(source, comet_ra, comet_dec)
        cx = MAP_CENTER[0] - east / MAP_PIXEL
        cy = MAP_CENTER[1] + north / MAP_PIXEL
        for x, y, weight in zip(cx, cy, weights, strict=True):
            x0, x1 = max(0, math.floor(x - radius)), min(
                shape[1], math.ceil(x + radius) + 1
            )
            y0, y1 = max(0, math.floor(y - radius)), min(
                shape[0], math.ceil(y + radius) + 1
            )
            if x0 >= x1 or y0 >= y1:
                continue
            yy, xx = np.ogrid[y0:y1, x0:x1]
            hole = (xx - x) ** 2 + (yy - y) ** 2 <= radius**2
            lost[y0:y1, x0:x1][hole] += weight
    return 1.0 - np.minimum(lost / weights.sum(), 1.0)


def correction(exposure, radius, source, background):
    source_mean = np.sum(exposure[source] / radius[source]) / exposure[source].sum()
    background_mean = (
        np.sum(exposure[background] / radius[background]) / exposure[background].sum()
    )
    return 1.0 / (1.0 - background_mean / source_mean)
