from pathlib import Path
from typing import NamedTuple

import numpy as np
from astropy.io import fits

from common import sas

PHYSICAL_PIXEL_ARCSEC = 0.05


class Spectrum(NamedTuple):
    channel: np.ndarray
    counts: np.ndarray
    variance: np.ndarray
    exposure: float
    backscal: float
    channel_ev: float


def read(path):
    with fits.open(path, memmap=False) as hdul:
        hdu = hdul["SPECTRUM"] if "SPECTRUM" in hdul else hdul[1]
        data, header = hdu.data, hdu.header
        names = set(data.dtype.names or ())
        exposure = float(header.get("EXPOSURE", 0))
        factor = exposure if "RATE" in names else 1.0
        field = "RATE" if "RATE" in names else "COUNTS"
        counts = np.asarray(data[field], float) * factor
        if "STAT_ERR" in names:
            errors = np.array(data["STAT_ERR"], float, copy=True) * factor
            channel = np.asarray(data["CHANNEL"], int)
            errors[(channel == 0) & (counts == 0)] = 0
            if np.any(errors < 0):
                raise ValueError(f"negative PHA error: {path}")
            variance = errors**2
        else:
            variance = np.abs(counts)
        return Spectrum(
            np.asarray(data["CHANNEL"], int),
            counts,
            variance,
            exposure,
            float(header.get("BACKSCAL", 1)),
            float(header.get("SPECDELT", 5)),
        )


def net(source, qpb):
    source, qpb = read(source), read(qpb)
    if not np.array_equal(source.channel, qpb.channel):
        raise ValueError("source and QPB channels differ")
    if (
        not np.isclose(source.exposure, qpb.exposure, rtol=1e-6, atol=1e-3)
        or not np.isclose(source.backscal, qpb.backscal, rtol=1e-6)
        or not np.isclose(source.channel_ev, qpb.channel_ev)
    ):
        raise ValueError("source and QPB normalization differs")
    counts = source.counts - qpb.counts
    variance = source.variance + qpb.variance
    return source._replace(counts=counts, variance=variance)


def stack(spectra):
    spectra = list(spectra)
    if not spectra:
        raise ValueError("no spectra")
    first = spectra[0]
    if any(not np.array_equal(first.channel, item.channel) for item in spectra[1:]):
        raise ValueError("stacked channels differ")
    if any(not np.isclose(first.channel_ev, item.channel_ev) for item in spectra[1:]):
        raise ValueError("stacked channel widths differ")
    exposure = sum(item.exposure for item in spectra)
    if exposure <= 0:
        raise ValueError("stacked exposure is zero")
    backscal = sum(item.exposure * item.backscal for item in spectra) / exposure
    return first._replace(
        counts=sum((item.counts for item in spectra), np.zeros_like(first.counts)),
        variance=sum(
            (item.variance for item in spectra), np.zeros_like(first.variance)
        ),
        exposure=exposure,
        backscal=backscal,
    )


def grouping(values, minimum=25):
    values = np.asarray(values, float)
    groups = np.full(values.size, -1, dtype=np.int16)
    if not groups.size:
        return groups
    groups[0] = 1
    total = 0.0
    last = 0
    for index, value in enumerate(values):
        total += max(value, 0)
        if total >= minimum and index + 1 < values.size:
            last = index + 1
            groups[last] = 1
            total = 0.0
    if total < minimum and last:
        groups[last] = -1
    return groups


def write(path, template, spectrum, rate=False, group_counts=None):
    with fits.open(template, memmap=False) as hdul:
        primary = fits.PrimaryHDU(header=hdul[0].header)
        old = hdul["SPECTRUM"] if "SPECTRUM" in hdul else hdul[1]
        header = old.header.copy()
    factor = spectrum.exposure if rate else 1.0
    name = "RATE" if rate else "COUNTS"
    columns = [
        fits.Column(name="CHANNEL", format="I", array=spectrum.channel),
        fits.Column(
            name=name, format="E", array=(spectrum.counts / factor).astype("f4")
        ),
        fits.Column(
            name="STAT_ERR",
            format="E",
            array=(np.sqrt(np.maximum(spectrum.variance, 0)) / factor).astype("f4"),
        ),
    ]
    if group_counts is not None:
        columns.append(
            fits.Column(name="GROUPING", format="I", array=grouping(group_counts))
        )
    hdu = fits.BinTableHDU.from_columns(columns, header=header)
    hdu.header.update(
        EXTNAME="SPECTRUM",
        EXPOSURE=spectrum.exposure,
        BACKSCAL=spectrum.backscal,
        BACKFILE="NONE",
        RESPFILE="NONE",
        ANCRFILE="NONE",
        POISSERR=False,
        HDUCLAS3="RATE" if rate else "COUNT",
        DETCHANS=spectrum.channel.size,
        TLMIN1=int(spectrum.channel.min()),
        TLMAX1=int(spectrum.channel.max()),
        SPECDELT=spectrum.channel_ev,
    )
    fits.HDUList([primary, hdu]).writeto(path, overwrite=True)


def add_response(inputs, weights, output, tool):
    total = float(sum(weights))
    if len(inputs) != len(weights) or total <= 0:
        raise ValueError("invalid response weights")
    listing = Path(output).with_suffix(".list")
    listing.write_text(
        "".join(
            f"{path} {weight / total:.10f}\n" for path, weight in zip(inputs, weights)
        )
    )
    argument = "rmffile" if tool == "addrmf" else "out_ARF"
    try:
        sas(
            tool,
            Path(output).parent,
            list=f"@{listing}",
            weights="-",
            **{argument: output, "clobber": True},
        )
    finally:
        listing.unlink(missing_ok=True)


def area_arcmin2(spectrum):
    return spectrum.backscal * (PHYSICAL_PIXEL_ARCSEC / 60) ** 2


def band_counts(spectrum, low_ev, high_ev):
    energy = spectrum.channel * spectrum.channel_ev
    return float(spectrum.counts[(energy >= low_ev) & (energy <= high_ev)].sum())


def response_counts(counts, response, template):
    with fits.open(response, memmap=True) as hdul:
        bounds = hdul["EBOUNDS"].data
        low = np.asarray(bounds["E_MIN"], float) * 1000
        high = np.asarray(bounds["E_MAX"], float) * 1000
    if (
        len(counts) != len(low)
        or not np.allclose(low[1:], high[:-1])
        or not np.all(np.diff(template.channel) == 1)
    ):
        raise ValueError("invalid response grid")
    source_edges = np.r_[low, high[-1]]
    target_edges = np.r_[
        template.channel * template.channel_ev,
        (template.channel[-1] + 1) * template.channel_ev,
    ]
    if (
        target_edges[0] < source_edges[0] - 1e-3
        or target_edges[-1] > source_edges[-1] + 1e-3
    ):
        raise ValueError("response does not cover the spectrum grid")
    cumulative = np.r_[0.0, np.cumsum(counts)]
    output = np.diff(np.interp(target_edges, source_edges, cumulative))
    if not np.isclose(output.sum(), np.sum(counts), rtol=1e-8, atol=1e-5):
        raise RuntimeError("response rebinning changed total counts")
    return output
