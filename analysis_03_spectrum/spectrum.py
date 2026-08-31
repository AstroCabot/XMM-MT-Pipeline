import os
from pathlib import Path
from typing import NamedTuple

import numpy as np
from astropy.io import fits


class Pha(NamedTuple):
    channel: np.ndarray
    counts: np.ndarray
    variance: np.ndarray
    exposure: float
    backscal: float
    width_ev: float


def read(path, poisson=True):
    with fits.open(path, memmap=False) as hdul:
        hdu = hdul["SPECTRUM"] if "SPECTRUM" in hdul else hdul[1]
        data, header = hdu.data, hdu.header
        names = set(data.names)
        exposure = float(header["EXPOSURE"])
        factor = exposure if "RATE" in names else 1.0
        counts = (
            np.asarray(data["RATE" if "RATE" in names else "COUNTS"], float) * factor
        )
        channel = np.asarray(data["CHANNEL"], int)
        if "STAT_ERR" in names:
            errors = np.array(data["STAT_ERR"], float, copy=True) * factor
            errors[(channel == 0) & (counts == 0)] = 0
            if np.any(errors < 0):
                raise ValueError(f"negative PHA error: {path}")
            variance = errors**2
        else:
            variance = np.abs(counts) if poisson else np.zeros_like(counts)
        return Pha(
            channel,
            counts,
            variance,
            exposure,
            float(header["BACKSCAL"]),
            float(header.get("SPECDELT", 5.0)),
        )


def write(path, template, pha, grouping=None, quality=None, links=None):
    with fits.open(template, memmap=False) as hdul:
        primary = hdul[0].copy()
        old = hdul["SPECTRUM"] if "SPECTRUM" in hdul else hdul[1]
        header = old.header.copy()
    columns = [
        fits.Column(name="CHANNEL", format="I", array=pha.channel),
        fits.Column(name="COUNTS", format="E", array=pha.counts.astype("f4")),
        fits.Column(
            name="STAT_ERR",
            format="E",
            array=np.sqrt(np.maximum(pha.variance, 0)).astype("f4"),
        ),
    ]
    if grouping is not None:
        columns += [
            fits.Column(name="GROUPING", format="I", array=grouping),
            fits.Column(name="QUALITY", format="I", array=quality),
        ]
    hdu = fits.BinTableHDU.from_columns(columns, header=header, name="SPECTRUM")
    hdu.header.update(
        EXPOSURE=pha.exposure,
        BACKSCAL=pha.backscal,
        AREASCAL=1.0,
        POISSERR=False,
        HDUCLAS3="COUNT",
        DETCHANS=len(pha.channel),
        TLMIN1=int(pha.channel.min()),
        TLMAX1=int(pha.channel.max()),
        SPECDELT=pha.width_ev,
        BACKFILE="NONE",
        RESPFILE="NONE",
        ANCRFILE="NONE",
    )
    if links:
        hdu.header.update(**links)
    fits.HDUList([primary, hdu]).writeto(path, overwrite=True)


def band_counts(path, low_ev, high_ev):
    pha = read(path)
    energy = pha.channel * pha.width_ev
    return float(pha.counts[(energy >= low_ev) & (energy <= high_ev)].sum())


def grouped_range(path, low_ev, high_ev):
    with fits.open(path, memmap=True) as hdul:
        hdu = hdul["SPECTRUM"]
        channel = np.asarray(hdu.data["CHANNEL"], int)
        grouping = np.asarray(hdu.data["GROUPING"], int)
        width = int(round(float(hdu.header["SPECDELT"])))
    if width <= 0:
        raise ValueError("invalid PHA channel width")
    first, stop = round(low_ev) // width, round(high_ev) // width + 1
    positions = {value: index for index, value in enumerate(channel)}
    if (
        len(positions) != len(channel)
        or first not in positions
        or stop not in positions
    ):
        raise ValueError("image band is outside the PHA channel grid")
    first, stop = positions[first], positions[stop]
    if grouping[first] != 1 or grouping[stop] != 1:
        raise ValueError("image band cuts a grouped spectral bin")
    ordinal = np.cumsum(grouping == 1)
    return int(ordinal[first]), int(ordinal[stop] - 1)


def background(path, source_path, qpb_path, sp_path, fractions, oot_path, oot_scale):
    source = read(source_path)
    total = np.zeros_like(source.counts)
    variance = np.zeros_like(source.counts)
    for component_path, fraction in ((qpb_path, fractions[0]), (sp_path, fractions[1])):
        component = read(component_path, poisson=False)
        if not np.array_equal(component.channel, source.channel):
            raise ValueError("background channel grid differs from source")
        scale = fraction * source.exposure / component.exposure
        total += scale * component.counts
        variance += scale**2 * component.variance
    if oot_path:
        oot = read(oot_path)
        if not np.array_equal(oot.channel, source.channel):
            raise ValueError("OoT channel grid differs from source")
        total += oot_scale * oot.counts
        variance += oot_scale**2 * oot.variance
    write(path, source_path, source._replace(counts=total, variance=variance))


def grouping(source, rmf, band_keV, minimum):
    pha = read(source)
    with fits.open(rmf, memmap=True) as hdul:
        bounds = hdul["EBOUNDS"].data
        channels = np.asarray(bounds["CHANNEL"], int)
        energy_min = np.asarray(bounds["E_MIN"])
        energy_max = np.asarray(bounds["E_MAX"])
    if not np.array_equal(channels, pha.channel):
        raise ValueError("RMF and PHA channels differ")
    low, high = np.asarray(band_keV, dtype=energy_min.dtype)
    # Match discrete PI <= high at the channel's lower edge.
    use = (energy_max > low) & (energy_min <= high)
    quality = np.where(use, 0, 5).astype(np.int16)
    groups = np.ones(len(pha.channel), np.int16)
    indices = np.flatnonzero(use)
    start, running = indices[0], 0.0
    for index in indices:
        groups[index] = 1 if index == start else -1
        running += max(pha.counts[index], 0)
        if running >= minimum and index != indices[-1]:
            start, running = index + 1, 0.0
    if running < minimum and start > indices[0]:
        groups[start] = -1
    return groups, quality


def combine_responses(paths, weights, output, tool, run_tool):
    if len(paths) != len(weights) or sum(weights) <= 0:
        raise ValueError("invalid response weights")
    listing = Path(output).with_suffix(".list")
    total = sum(weights)
    listing.write_text(
        "".join(
            f"{os.path.relpath(path, listing.parent)} {weight / total:.12g}\n"
            for path, weight in zip(paths, weights)
        )
    )
    key = "rmffile" if tool == "addrmf" else "out_ARF"
    try:
        run_tool(
            tool,
            listing.parent,
            list=f"@{listing.name}",
            weights="-",
            **{key: Path(output).name, "clobber": "yes"},
        )
    finally:
        listing.unlink(missing_ok=True)


def stack(rows, output, stem, band_keV, minimum, run_tool):
    output.mkdir(parents=True, exist_ok=True)
    source = [read(row["source_pha"]) for row in rows]
    back = [read(row["background_pha"], poisson=False) for row in rows]
    first = source[0]
    if any(not np.array_equal(first.channel, item.channel) for item in source + back):
        raise ValueError("stacked channel grids differ")
    exposure = sum(item.exposure for item in source)
    backscal = sum(item.exposure * item.backscal for item in source) / exposure
    summed_source = first._replace(
        counts=sum((item.counts for item in source), np.zeros_like(first.counts)),
        variance=sum((item.variance for item in source), np.zeros_like(first.variance)),
        exposure=exposure,
        backscal=backscal,
    )
    summed_back = summed_source._replace(
        counts=sum((item.counts for item in back), np.zeros_like(first.counts)),
        variance=sum((item.variance for item in back), np.zeros_like(first.variance)),
    )
    raw = output / f"{stem}-raw.pi"
    write(raw, rows[0]["source_pha"], summed_source)
    bkg = output / f"{stem}-background.pi"
    rmf, arf = output / f"{stem}.rmf", output / f"{stem}.arf"
    write(bkg, rows[0]["source_pha"], summed_back)
    weights = [float(row["exposure_s"]) * float(row["omega_arcmin2"]) for row in rows]
    combine_responses([row["rmf"] for row in rows], weights, rmf, "addrmf", run_tool)
    combine_responses([row["arf"] for row in rows], weights, arf, "addarf", run_tool)
    groups, quality = grouping(raw, rmf, band_keV, minimum)
    source_path = output / f"{stem}-source.pi"
    write(
        source_path,
        raw,
        summed_source,
        groups,
        quality,
        {"BACKFILE": bkg.name, "RESPFILE": rmf.name, "ANCRFILE": arf.name},
    )
    raw.unlink()
    omega = (
        sum(float(r["exposure_s"]) * float(r["omega_arcmin2"]) for r in rows) / exposure
    )
    return source_path, bkg, rmf, arf, exposure, backscal, omega
