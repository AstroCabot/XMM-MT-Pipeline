#!/usr/bin/env python3
"""
PPS / MT_PPS analytics — extract and summarize the official pipeline products
for comparison with our custom comet pipeline.

Produces:
  <outdir>/pps_summary.txt       — text report
  <outdir>/pps_gti_timeline.png  — GTI intervals per instrument overlaid on flare rates
  <outdir>/pps_attitude.png      — attitude RA/DEC vs time
  <outdir>/pps_flare_rates.png   — per-instrument flare background light curves + thresholds
  <outdir>/pps_source_map.png    — detected source positions on soft-band image
  <outdir>/pps_image_bands.png   — multi-band image mosaic from PPS products

Usage:
  python3 pps_analytics.py /path/to/PPS [--mt /path/to/MT_PPS] [--workdir /path/to/pipeline/output] [-o outdir]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def glob_ftz(directory, pattern):
    """Glob for .FTZ files, ignoring .gpg signatures."""
    return sorted(directory.glob(pattern))


def inst_label(fname):
    """Extract instrument label from PPS filename."""
    name = fname.name
    if "PNS" in name or "PNU" in name or "PNX" in name:
        return "PN"
    if "M1S" in name or "M1U" in name or "M1X" in name:
        return "MOS1"
    if "M2S" in name or "M2U" in name or "M2X" in name:
        return "MOS2"
    return "?"


def is_scheduled(fname):
    """True if the exposure is scheduled (S) rather than unscheduled (U)."""
    name = fname.name
    for tag in ("PNS", "M1S", "M2S", "R1S", "R2S"):
        if tag in name:
            return True
    return False


def mjd_to_label(t_mjd):
    """Convert MJD to 'HH:MM' string."""
    tt = Time(t_mjd, format="mjd")
    return tt.iso[11:16]


# ---------------------------------------------------------------------------
# GTI extraction
# ---------------------------------------------------------------------------
def extract_gtis(evtfile):
    """Return dict of {extname: (start_arr, stop_arr)} for all GTI extensions."""
    gtis = {}
    with fits.open(evtfile) as hdul:
        for ext in hdul:
            if "GTI" in ext.name.upper() or ext.name.startswith("STDGTI"):
                if hasattr(ext, "data") and ext.data is not None and len(ext.data) > 0:
                    gtis[ext.name] = (ext.data["START"].copy(), ext.data["STOP"].copy())
    return gtis


def extract_event_time_range(evtfile):
    """Return (tmin, tmax) from the EVENTS extension."""
    with fits.open(evtfile) as hdul:
        t = hdul["EVENTS"].data["TIME"]
        return t.min(), t.max()


# ---------------------------------------------------------------------------
# flare rate extraction
# ---------------------------------------------------------------------------
def extract_flare_rate(fbkfile):
    """Return dict with time, rate, error, and SRC_GTIS if present."""
    result = {}
    with fits.open(fbkfile) as hdul:
        rate = hdul["RATE"].data
        result["time"] = rate["TIME"].copy()
        result["rate"] = rate["RATE"].copy()
        result["error"] = rate["ERROR"].copy()
        if "SRC_GTIS" in hdul:
            gti = hdul["SRC_GTIS"].data
            result["gti_start"] = gti["START"].copy()
            result["gti_stop"] = gti["STOP"].copy()
            result["gti_clean_ks"] = np.sum(gti["STOP"] - gti["START"]) / 1000
    return result


# ---------------------------------------------------------------------------
# attitude extraction
# ---------------------------------------------------------------------------
def extract_attitude(attfile):
    """Return time, AHFRA, AHFDEC from the ATTHK extension."""
    with fits.open(attfile) as hdul:
        data = hdul["ATTHK"].data
        return {
            "time": data["TIME"].copy(),
            "ra": data["AHFRA"].copy(),
            "dec": data["AHFDEC"].copy(),
            "pa": data["AHFPA"].copy(),
        }


# ---------------------------------------------------------------------------
# source list extraction
# ---------------------------------------------------------------------------
def extract_source_list(srcfile):
    """Return source list table columns."""
    with fits.open(srcfile) as hdul:
        data = hdul["SRCLIST"].data
        cols = data.columns.names
        result = {"columns": cols}
        for c in (
            "RA",
            "DEC",
            "SRC_NUM",
            "PN_CTS",
            "M1_CTS",
            "M2_CTS",
            "EP_EXTENT",
            "EP_EXTENT_ERR",
            "EP_DET_ML",
            "EP_8_DET_ML",
        ):
            if c in cols:
                result[c] = data[c].copy()
        return result


# ---------------------------------------------------------------------------
# image extraction
# ---------------------------------------------------------------------------
def extract_image(imgfile):
    """Return (data_2d, wcs_header_dict) from a PPS IMAGE file."""
    with fits.open(imgfile) as hdul:
        hdr = hdul[0].header
        data = hdul[0].data
        wcs_info = {}
        for k in (
            "CRVAL1",
            "CRVAL2",
            "CDELT1",
            "CDELT2",
            "CRPIX1",
            "CRPIX2",
            "NAXIS1",
            "NAXIS2",
        ):
            if k in hdr:
                wcs_info[k] = hdr[k]
        return data, wcs_info


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------
def plot_gti_timeline(evtfiles, fbkfiles, tref, outpath):
    """Plot GTI bars per instrument overlaid on flare rates."""
    fig, axes = plt.subplots(
        2, 1, figsize=(14, 7), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
    )

    # TOP: flare light curves
    ax = axes[0]
    colors_inst = {"PN": "C0", "MOS1": "C1", "MOS2": "C2"}
    for fbk in fbkfiles:
        inst = inst_label(fbk)
        if not is_scheduled(fbk):
            continue
        fr = extract_flare_rate(fbk)
        t_rel = (fr["time"] - tref) / 1000
        good = np.isfinite(fr["rate"])
        ax.plot(
            t_rel[good],
            fr["rate"][good],
            lw=0.5,
            alpha=0.8,
            color=colors_inst.get(inst, "gray"),
            label=inst,
        )
    ax.set_ylabel("Background rate (ct/s)")
    ax.set_title("PPS flare background rates (scheduled exposures)")
    ax.legend(fontsize=8)
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)

    # BOTTOM: GTI bars
    ax = axes[1]
    y_positions = {}
    y_idx = 0
    for evtf in evtfiles:
        inst = inst_label(evtf)
        sched = "S" if is_scheduled(evtf) else "U"
        label = f"{inst} ({sched})"
        if label not in y_positions:
            y_positions[label] = y_idx
            y_idx += 1
        y = y_positions[label]

        gtis = extract_gtis(evtf)
        # Use first GTI extension (usually STDGTI01, the main one for CCD 1)
        for gti_name, (starts, stops) in gtis.items():
            for s, e in zip(starts, stops):
                ax.barh(
                    y,
                    (e - s) / 1000,
                    left=(s - tref) / 1000,
                    height=0.6,
                    color=colors_inst.get(inst, "gray"),
                    alpha=0.4,
                    edgecolor="none",
                )
            break  # only first GTI ext

    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(list(y_positions.keys()), fontsize=8)
    ax.set_xlabel("Time from first event (ks)")
    ax.set_ylabel("Instrument")
    ax.set_title("PPS Good Time Intervals (first CCD GTI per event list)")
    ax.grid(True, axis="x", alpha=0.3)

    fig.tight_layout()
    fig.savefig(str(outpath), dpi=150)
    plt.close(fig)


def plot_attitude(attdata, tref, outpath):
    """Plot RA/DEC vs time from attitude timeseries."""
    t_rel = (attdata["time"] - tref) / 1000

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 5), sharex=True)
    ax1.plot(t_rel, attdata["ra"], lw=0.5, color="C0")
    ax1.set_ylabel("RA (deg)")
    ax1.grid(True, alpha=0.3)
    ax1.set_title("PPS attitude timeseries (ATTHK)")

    ax2.plot(t_rel, attdata["dec"], lw=0.5, color="C1")
    ax2.set_ylabel("DEC (deg)")
    ax2.set_xlabel("Time from first event (ks)")
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(str(outpath), dpi=150)
    plt.close(fig)


def plot_flare_rates(fbkfiles, tref, outpath):
    """Detailed per-instrument flare rate plots with clean-GTI shading."""
    scheduled = [f for f in fbkfiles if is_scheduled(f)]
    n = len(scheduled)
    if n == 0:
        return

    fig, axes = plt.subplots(n, 1, figsize=(14, 3 * n), sharex=True)
    if n == 1:
        axes = [axes]

    for ax, fbk in zip(axes, scheduled):
        inst = inst_label(fbk)
        fr = extract_flare_rate(fbk)
        t_rel = (fr["time"] - tref) / 1000
        good = np.isfinite(fr["rate"])

        ax.plot(t_rel[good], fr["rate"][good], lw=0.5, color="C0", label="rate")

        # Shade clean GTI intervals
        if "gti_start" in fr:
            for s, e in zip(fr["gti_start"], fr["gti_stop"]):
                ax.axvspan(
                    (s - tref) / 1000, (e - tref) / 1000, alpha=0.1, color="green"
                )
            ax.text(
                0.99,
                0.95,
                f"Clean: {fr['gti_clean_ks']:.1f} ks",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=8,
                bbox=dict(facecolor="white", alpha=0.7),
            )

        finite_rates = fr["rate"][good]
        if len(finite_rates) > 0:
            med = np.median(finite_rates)
            ax.axhline(med, color="k", ls="--", lw=0.7, alpha=0.6)

        ax.set_ylabel("Rate (ct/s)")
        ax.set_title(f"{inst} — {fbk.name}", fontsize=9)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Time from first event (ks)")
    fig.suptitle(
        "PPS flare background rates with clean GTI shading (green)", fontsize=11
    )
    fig.tight_layout()
    fig.savefig(str(outpath), dpi=150)
    plt.close(fig)


def plot_source_map(srcdata, imgfile, outpath):
    """Plot detected sources on the soft-band image."""
    if imgfile is None or not imgfile.exists():
        return

    img_data, wcs = extract_image(imgfile)
    if img_data is None:
        return

    # Build coordinate grid
    nx = wcs.get("NAXIS1", img_data.shape[1])
    ny = wcs.get("NAXIS2", img_data.shape[0])
    x_grid = wcs["CRVAL1"] + wcs["CDELT1"] * (
        np.arange(nx) - wcs.get("CRPIX1", nx // 2) + 1
    )
    y_grid = wcs["CRVAL2"] + wcs["CDELT2"] * (
        np.arange(ny) - wcs.get("CRPIX2", ny // 2) + 1
    )

    fig, ax = plt.subplots(figsize=(10, 8))
    vmax = np.nanpercentile(img_data[img_data > 0], 99) if np.any(img_data > 0) else 1
    ax.pcolormesh(x_grid, y_grid, img_data, vmin=0, vmax=vmax, cmap="inferno")

    if "RA" in srcdata and "DEC" in srcdata:
        ax.scatter(
            srcdata["RA"],
            srcdata["DEC"],
            s=60,
            facecolors="none",
            edgecolors="cyan",
            lw=1,
            label=f'{len(srcdata["RA"])} sources',
        )
        for i in range(len(srcdata["RA"])):
            label = str(srcdata.get("SRC_NUM", np.arange(len(srcdata["RA"])))[i])
            ax.annotate(
                label,
                (srcdata["RA"][i], srcdata["DEC"][i]),
                color="cyan",
                fontsize=6,
                xytext=(3, 3),
                textcoords="offset points",
            )

    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("DEC (deg)")
    ax.set_title("PPS detected sources on broad-band image")
    ax.legend(fontsize=8)
    ax.invert_xaxis()
    fig.tight_layout()
    fig.savefig(str(outpath), dpi=150)
    plt.close(fig)


def plot_image_bands(pps_dir, inst_prefix, outpath):
    """Plot multi-band image mosaic for a given instrument."""
    band_map = {
        "1000": "0.2-0.5 keV",
        "2000": "0.5-2.0 keV",
        "3000": "2.0-4.5 keV",
        "4000": "4.5-7.5 keV",
        "5000": "7.5-12 keV",
        "8000": "0.2-12 keV",
    }

    images = []
    labels = []
    wcs_list = []
    for band_code, band_label in band_map.items():
        pattern = f"*{inst_prefix}*IMAGE_{band_code}.FTZ"
        files = glob_ftz(pps_dir, pattern)
        if files:
            data, wcs = extract_image(files[0])
            if data is not None:
                images.append(data)
                labels.append(band_label)
                wcs_list.append(wcs)

    if not images:
        return

    ncols = min(3, len(images))
    nrows = int(np.ceil(len(images) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows))
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes[np.newaxis, :]
    elif ncols == 1:
        axes = axes[:, np.newaxis]

    for idx in range(nrows * ncols):
        r, c = divmod(idx, ncols)
        ax = axes[r, c]
        if idx >= len(images):
            ax.axis("off")
            continue

        img = images[idx]
        wcs = wcs_list[idx]
        nx = wcs.get("NAXIS1", img.shape[1])
        ny = wcs.get("NAXIS2", img.shape[0])
        x_grid = wcs["CRVAL1"] + wcs["CDELT1"] * (
            np.arange(nx) - wcs.get("CRPIX1", nx // 2) + 1
        )
        y_grid = wcs["CRVAL2"] + wcs["CDELT2"] * (
            np.arange(ny) - wcs.get("CRPIX2", ny // 2) + 1
        )

        vmax = np.nanpercentile(img[img > 0], 99) if np.any(img > 0) else 1
        ax.pcolormesh(x_grid, y_grid, img, vmin=0, vmax=vmax, cmap="inferno")
        ax.set_title(labels[idx], fontsize=9)
        ax.set_xlabel("RA (deg)", fontsize=7)
        ax.set_ylabel("DEC (deg)", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.invert_xaxis()

    fig.suptitle(f"PPS {inst_prefix} images by energy band", fontsize=11)
    fig.tight_layout()
    fig.savefig(str(outpath), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# text summary
# ---------------------------------------------------------------------------
def write_summary(outpath, pps_dir, mt_dir, evtfiles, fbkfiles, attdata, srcdata):
    """Write a text summary of all extracted PPS information."""
    lines = []
    lines.append("=" * 72)
    lines.append("PPS / MT_PPS Analytics Summary")
    lines.append("=" * 72)
    lines.append(f"PPS dir:    {pps_dir}")
    if mt_dir:
        lines.append(f"MT_PPS dir: {mt_dir}")
    lines.append("")

    # Event lists
    lines.append("--- Event Lists ---")
    for evtf in evtfiles:
        inst = inst_label(evtf)
        sched = "scheduled" if is_scheduled(evtf) else "unscheduled"
        tmin, tmax = extract_event_time_range(evtf)
        gtis = extract_gtis(evtf)
        lines.append(f"  {evtf.name}")
        lines.append(f"    Instrument: {inst} ({sched})")
        lines.append(f"    Time range: {(tmax - tmin)/1000:.1f} ks")
        with fits.open(evtf) as hdul:
            nevt = len(hdul["EVENTS"].data)
            lines.append(f"    Events: {nevt:,}")
        for gname, (starts, stops) in gtis.items():
            total_ks = np.sum(stops - starts) / 1000
            lines.append(
                f"    {gname}: {len(starts)} intervals, {total_ks:.1f} ks clean"
            )
        lines.append("")

    # Flare background
    lines.append("--- Flare Background (scheduled exposures) ---")
    for fbk in fbkfiles:
        if not is_scheduled(fbk):
            continue
        inst = inst_label(fbk)
        fr = extract_flare_rate(fbk)
        good = np.isfinite(fr["rate"])
        finite_rates = fr["rate"][good]
        lines.append(f"  {fbk.name} ({inst})")
        if len(finite_rates) > 0:
            lines.append(
                f"    Rate range: {finite_rates.min():.2f} - {finite_rates.max():.2f} ct/s"
            )
            lines.append(f"    Median rate: {np.median(finite_rates):.2f} ct/s")
        if "gti_clean_ks" in fr:
            lines.append(f"    Clean exposure (SRC_GTIS): {fr['gti_clean_ks']:.1f} ks")
        lines.append("")

    # Attitude
    if attdata:
        lines.append("--- Attitude ---")
        dt_ks = (attdata["time"][-1] - attdata["time"][0]) / 1000
        lines.append(f"  Duration: {dt_ks:.1f} ks")
        lines.append(
            f"  RA range:  {attdata['ra'].min():.4f} - {attdata['ra'].max():.4f} deg "
            f"(span {(attdata['ra'].max() - attdata['ra'].min()) * 3600:.1f}\")"
        )
        lines.append(
            f"  DEC range: {attdata['dec'].min():.4f} - {attdata['dec'].max():.4f} deg "
            f"(span {(attdata['dec'].max() - attdata['dec'].min()) * 3600:.1f}\")"
        )
        lines.append(
            f"  PA range:  {attdata['pa'].min():.4f} - {attdata['pa'].max():.4f} deg"
        )
        lines.append("")

    # Source list
    if srcdata and "RA" in srcdata:
        lines.append("--- Detected Sources (EPIC) ---")
        n = len(srcdata["RA"])
        lines.append(f"  Total sources: {n}")
        if "EP_DET_ML" in srcdata:
            ml = srcdata["EP_DET_ML"]
            lines.append(f"  Detection ML range: {ml.min():.1f} - {ml.max():.1f}")
        if "PN_CTS" in srcdata:
            cts = srcdata["PN_CTS"]
            lines.append(f"  PN counts range: {cts.min():.0f} - {cts.max():.0f}")
        lines.append(f"  Source positions (RA, DEC):")
        for i in range(n):
            num = srcdata.get("SRC_NUM", np.arange(n))[i]
            ra = srcdata["RA"][i]
            dec = srcdata["DEC"][i]
            extra = ""
            if "PN_CTS" in srcdata:
                extra += f", PN_cts={srcdata['PN_CTS'][i]:.0f}"
            if "EP_DET_ML" in srcdata:
                extra += f", ML={srcdata['EP_DET_ML'][i]:.1f}"
            lines.append(f"    #{num:3d}  RA={ra:.5f}  DEC={dec:.5f}{extra}")
        lines.append("")

    text = "\n".join(lines)
    outpath.write_text(text)
    return text


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="PPS / MT_PPS analytics")
    parser.add_argument("pps_dir", type=Path, help="Path to PPS directory")
    parser.add_argument(
        "--mt", type=Path, default=None, help="Path to MT_PPS directory"
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=None,
        help="Pipeline WORKDIR for comparison (optional)",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        required=True,
        help="Output directory (must NOT be inside PPS or MT_PPS)",
    )
    args = parser.parse_args()

    pps_dir = args.pps_dir.resolve()
    mt_dir = args.mt.resolve() if args.mt else None
    outdir = args.outdir.resolve()

    # Safety: refuse to write inside any input data directory
    for data_dir in [pps_dir, mt_dir]:
        if data_dir and (outdir == data_dir or data_dir in outdir.parents):
            print(
                f"ERROR: output directory {outdir} is inside data directory {data_dir}",
                file=sys.stderr,
            )
            print(
                "Refusing to write into original data. Use a different -o path.",
                file=sys.stderr,
            )
            return 1

    outdir.mkdir(parents=True, exist_ok=True)

    # Use MT_PPS preferentially for event lists and source lists (richer products)
    primary_dir = mt_dir if mt_dir and mt_dir.exists() else pps_dir

    print(f"PPS dir:    {pps_dir}")
    print(f"MT_PPS dir: {mt_dir}")
    print(f"Primary:    {primary_dir}")
    print(f"Output:     {outdir}")
    print()

    # ---- Gather files ----
    evtfiles = glob_ftz(primary_dir, "*PIEVLI*.FTZ") + glob_ftz(
        primary_dir, "*MIEVLI*.FTZ"
    )
    fbkfiles_primary = glob_ftz(primary_dir, "*FBKTSR*.FTZ")
    # Also get from PPS if MT is primary (PPS may have merged-instrument FBKTSR)
    fbkfiles = fbkfiles_primary
    if mt_dir and mt_dir != pps_dir:
        fbkfiles += glob_ftz(pps_dir, "*FBKTSR*.FTZ")
    # Deduplicate by filename
    seen = set()
    deduped = []
    for f in fbkfiles:
        if f.name not in seen:
            seen.add(f.name)
            deduped.append(f)
    fbkfiles = deduped

    attfiles = glob_ftz(primary_dir, "*ATTTSR*.FTZ")
    srcfiles = glob_ftz(primary_dir, "*OBSMLI*.FTZ")
    # Filter to EPIC source list (EPX)
    epic_srcfiles = [f for f in srcfiles if "EPX" in f.name]

    print(f"Event lists:  {len(evtfiles)}")
    print(f"FBKTSR files: {len(fbkfiles)}")
    print(f"Attitude:     {len(attfiles)}")
    print(f"Source lists:  {len(srcfiles)} ({len(epic_srcfiles)} EPIC)")
    print()

    # ---- Reference time ----
    tref = None
    if evtfiles:
        tmin, _ = extract_event_time_range(evtfiles[0])
        tref = tmin

    # ---- Attitude ----
    attdata = None
    if attfiles:
        attdata = extract_attitude(attfiles[0])
        if tref is None:
            tref = attdata["time"][0]
        print("Attitude extracted")
        plot_attitude(attdata, tref, outdir / "pps_attitude.png")
        print(f"  -> {outdir / 'pps_attitude.png'}")

    # ---- Source list ----
    srcdata = None
    if epic_srcfiles:
        srcdata = extract_source_list(epic_srcfiles[0])
        print(f"Source list: {len(srcdata.get('RA', []))} EPIC sources")

    # ---- GTI timeline ----
    if evtfiles:
        plot_gti_timeline(evtfiles, fbkfiles, tref, outdir / "pps_gti_timeline.png")
        print(f"  -> {outdir / 'pps_gti_timeline.png'}")

    # ---- Flare rates ----
    if fbkfiles:
        plot_flare_rates(fbkfiles, tref, outdir / "pps_flare_rates.png")
        print(f"  -> {outdir / 'pps_flare_rates.png'}")

    # ---- Source map on image ----
    if srcdata:
        # Find a broad-band image for background
        broad_imgs = glob_ftz(primary_dir, "*IMAGE_8000.FTZ")
        pn_imgs = [f for f in broad_imgs if "PNS" in f.name or "PNU" in f.name]
        bg_img = pn_imgs[0] if pn_imgs else (broad_imgs[0] if broad_imgs else None)
        if bg_img:
            plot_source_map(srcdata, bg_img, outdir / "pps_source_map.png")
            print(f"  -> {outdir / 'pps_source_map.png'}")

    # ---- Multi-band images ----
    for prefix in ("PNS", "M1S", "M2S"):
        imgs = glob_ftz(primary_dir, f"*{prefix}*IMAGE_*.FTZ")
        if imgs:
            plot_image_bands(
                primary_dir, prefix, outdir / f"pps_image_bands_{prefix}.png"
            )
            print(f"  -> {outdir / f'pps_image_bands_{prefix}.png'}")

    # ---- Text summary ----
    summary = write_summary(
        outdir / "pps_summary.txt",
        pps_dir,
        mt_dir,
        evtfiles,
        fbkfiles,
        attdata,
        srcdata,
    )
    print(f"\n  -> {outdir / 'pps_summary.txt'}")
    print()
    print(summary)

    return 0


if __name__ == "__main__":
    sys.exit(main())
