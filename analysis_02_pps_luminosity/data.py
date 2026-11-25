import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

from .selection import event_mjd, read_gtis, time_samples

HERE = Path(__file__).resolve().parent
REPRO = HERE.parent
SHARED = json.loads((REPRO / "parameters.json").read_text())
SETTINGS = json.loads((REPRO / "analysis_01_reduction/settings.json").read_text())
PPS = Path(SETTINGS["mt_pps"])
PREFIX = "P0963720201PNS003"
NATIVE = ((200, 500), (500, 1000), (1000, 2000))


def pattern_max(band):
    return 0 if band[1] <= 500 else 4


def gti_mask(events, gtis):
    keep = np.zeros(len(events), bool)
    for ccd, intervals in enumerate(gtis, 1):
        rows = np.flatnonzero(events["CCDNR"] == ccd)
        times = events["TIME"][rows]
        accepted = np.zeros(len(rows), bool)
        for start, stop in intervals:
            accepted |= (times >= start) & (times <= stop)
        keep[rows] = accepted
    return keep


def selected(energy, low, high):
    # Native PPS images use (low, high], with no overlap between adjacent bands.
    return (energy > low) & (energy <= high)


class Dataset:
    def __init__(self, edges):
        self.edges = np.asarray(edges)
        self.inputs = [
            PPS / (PREFIX + name) for name in ("PIEVLI0000.FTZ", "IMAGE_2000.FTZ")
        ]
        gtis = read_gtis(self.inputs[1])
        with fits.open(self.inputs[0], memmap=True) as hdul:
            head = hdul[0].header
            if not head.get("MTFLAG"):
                raise ValueError("PPS events are not in the comet frame")
            events = hdul["EVENTS"].data
            header = hdul["EVENTS"].header.copy()
            keep = (
                (events["RAWY"] > 12)
                & (events["PATTERN"] <= 4)
                & (events["FLAG"] == 0)
                & gti_mask(events, gtis)
            )
            x0 = (
                head["REFXCRPX"]
                + (head["MRA_OBJ"] - head["REFXCRVL"])
                * np.cos(np.deg2rad(head["MDEC_OBJ"]))
                / head["REFXCDLT"]
            )
            y0 = (
                head["REFYCRPX"]
                + (head["MDEC_OBJ"] - head["REFYCRVL"]) / head["REFYCDLT"]
            )
            self.radius = np.hypot(events["X"][keep] - x0, events["Y"][keep] - y0)
            self.radius *= abs(head["REFYCDLT"]) * 60
            self.pi = np.asarray(events["PI"][keep], int)
            self.pattern = np.asarray(events["PATTERN"][keep], int)
            center = head["MRA_OBJ"], head["MDEC_OBJ"]
        times, weights = time_samples(gtis)
        epoch = Time(
            np.average(event_mjd(times, header), weights=weights),
            format="mjd",
            scale="utc",
        )
        self.epoch = epoch.isot
        self.good_time = float(np.mean([np.sum(g[:, 1] - g[:, 0]) for g in gtis]))
        hard = selected(self.pi, *SHARED["hard_band_ev"])
        self.hard = np.histogram(self.radius[hard], edges)[0]
        self.maps = []
        grid = None
        for number in range(1, 4):
            path = PPS / f"{PREFIX}EXPMAP{number}000.FTZ"
            self.inputs.append(path)
            image, header = fits.getdata(path, header=True)
            if "RADECSYS" in header:
                header["RADESYS"] = header.pop("RADECSYS")
            wcs = WCS(header, fix=False)
            position = wcs.world_to_pixel_values(*center)
            current = (
                image.shape,
                tuple(position),
                tuple(wcs.pixel_scale_matrix.ravel()),
            )
            if grid is not None and current != grid:
                raise ValueError("PPS exposure maps use different grids")
            grid = current
            yy, xx = np.indices(image.shape)
            self.map_radius = np.hypot(xx - position[0], yy - position[1])
            self.map_radius *= abs(header["CDELT2"]) * 60
            self.pixel_area = abs(np.linalg.det(wcs.pixel_scale_matrix)) * 3600
            self.maps.append(np.asarray(image, float))

    def band(self, low, high, singles=False):
        chosen = selected(self.pi, low, high) & (
            (self.pattern == 0) if singles else True
        )
        counts = np.histogram(self.radius[chosen], self.edges)[0]
        rates = []
        inside = self.map_radius < self.edges[-1]
        for (a, b), image in zip(NATIVE, self.maps):
            lo, hi = max(a, low), min(b, high)
            count = (
                np.count_nonzero(
                    chosen & selected(self.pi, lo, hi) & (self.radius < self.edges[-1])
                )
                if lo < hi
                else 0
            )
            rates.append(count / (image[inside].sum() * self.pixel_area))
        weights = np.asarray(rates) / np.sum(rates)
        exposure = sum(w * image for w, image in zip(weights, self.maps))
        ring = np.searchsorted(self.edges, self.map_radius, side="right") - 1
        use = (ring >= 0) & (ring < len(counts)) & (exposure > 0)
        weight = exposure[use] * self.pixel_area
        return {
            "soft": counts,
            "hard": self.hard.copy(),
            "area": np.bincount(ring[use], weight, minlength=len(counts)),
            "radius": self.map_radius[use],
            "pixel_area_arcmin2": self.pixel_area,
            "weight": weight,
            "ring": ring[use],
            "exposure_weights": weights,
        }
