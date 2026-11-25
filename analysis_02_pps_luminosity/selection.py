import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta


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
