#!/usr/bin/env python3
"""Fast comet-frame diagnostic: time-binned images and light curve.

Reads comet-frame event files directly and produces quick-look
images in time bins to spot flare contamination and centring issues.
"""
import argparse
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def load_events(evt_path, pimin, pimax):
    with fits.open(evt_path) as hdul:
        data = hdul["EVENTS"].data
        pi = data["PI"]
        mask = (pi >= pimin) & (pi <= pimax)
        return (data["TIME"][mask], data["X"][mask], data["Y"][mask])


def make_image(x, y, xmin, xmax, ymin, ymax, binsize):
    nx = int(np.ceil((xmax - xmin) / binsize))
    ny = int(np.ceil((ymax - ymin) / binsize))
    img, _, _ = np.histogram2d(x, y, bins=[nx, ny], range=[[xmin, xmax], [ymin, ymax]])
    return img.T


def main():
    parser = argparse.ArgumentParser(description="Fast comet-frame diagnostic")
    parser.add_argument("workdir", help="Pipeline WORKDIR (or /tmp/ shortlink)")
    parser.add_argument(
        "--nbins", type=int, default=10, help="Number of time bins (default: 10)"
    )
    parser.add_argument("--pimin", type=int, default=200, help="PI min (default: 200)")
    parser.add_argument(
        "--pimax", type=int, default=2000, help="PI max (default: 2000)"
    )
    parser.add_argument(
        "--bin-arcsec",
        type=float,
        default=8.0,
        help="Image pixel size in arcsec (default: 8)",
    )
    parser.add_argument(
        "--inst", nargs="+", default=["PN", "M1", "M2"], help="Instruments to use"
    )
    parser.add_argument(
        "--radius-arcsec", type=float, default=900.0, help="Image half-width in arcsec"
    )
    parser.add_argument(
        "--src-r-arcsec",
        type=float,
        default=500.0,
        help="Source aperture radius for overlay",
    )
    parser.add_argument(
        "--lc-binsize",
        type=float,
        default=1000.0,
        help="Light curve bin size in seconds",
    )
    args = parser.parse_args()
    workdir = Path(args.workdir)
    outdir = workdir / "qc" / "quick_diag"
    outdir.mkdir(parents=True, exist_ok=True)
    PHYS_PER_ARCSEC = 1.0 / 0.05
    binsize_phys = args.bin_arcsec * PHYS_PER_ARCSEC
    all_t, all_x, all_y = ([], [], [])
    per_inst = {}
    for inst in args.inst:
        evt = workdir / "comet" / f"{inst}_comet.fits"
        if not evt.exists():
            print(f"  {inst}: not found, skipping")
            continue
        t, x, y = load_events(str(evt), args.pimin, args.pimax)
        print(f"  {inst}: {len(t):,} events in PI [{args.pimin}:{args.pimax}]")
        all_t.append(t)
        all_x.append(x)
        all_y.append(y)
        per_inst[inst] = (t, x, y)
    if not all_t:
        print("No events found!")
        return 1
    all_t = np.concatenate(all_t)
    all_x = np.concatenate(all_x)
    all_y = np.concatenate(all_y)
    tmin, tmax = (all_t.min(), all_t.max())
    duration_ks = (tmax - tmin) / 1000
    print(f"\nTotal: {len(all_t):,} events, {duration_ks:.1f} ks")
    xcen = np.median(all_x)
    ycen = np.median(all_y)
    half_phys = args.radius_arcsec * PHYS_PER_ARCSEC
    xmin, xmax = (xcen - half_phys, xcen + half_phys)
    ymin, ymax = (ycen - half_phys, ycen + half_phys)
    print("\nGenerating light curve...")
    lc_bins = np.arange(tmin, tmax + args.lc_binsize, args.lc_binsize)
    lc_counts, _ = np.histogram(all_t, bins=lc_bins)
    lc_rates = lc_counts / args.lc_binsize
    lc_midtimes = 0.5 * (lc_bins[:-1] + lc_bins[1:])
    lc_t_rel = (lc_midtimes - tmin) / 1000
    median_rate = np.median(lc_rates)
    flare_thresh = 3 * median_rate
    fig, ax = plt.subplots(figsize=(12, 3))
    colors = ["#d62728" if r > flare_thresh else "#1f77b4" for r in lc_rates]
    ax.bar(
        lc_t_rel, lc_rates, width=args.lc_binsize / 1000, color=colors, edgecolor="none"
    )
    ax.axhline(
        median_rate, color="k", ls="--", lw=0.8, label=f"median = {median_rate:.1f}"
    )
    ax.axhline(
        flare_thresh,
        color="red",
        ls=":",
        lw=0.8,
        label=f"3× median = {flare_thresh:.1f}",
    )
    ax.set_xlabel("Time from start (ks)")
    ax.set_ylabel(f"EPIC rate (ct/s, PI {args.pimin}-{args.pimax})")
    ax.set_title("Count rate — red bins = flare (>3× median)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(str(outdir / "lightcurve.png"), dpi=150)
    plt.close(fig)
    print(f"  -> {outdir / 'lightcurve.png'}")
    bin_edges = np.linspace(tmin, tmax, args.nbins + 1)
    print(f"\nGenerating {args.nbins} time-bin counts images...")
    images = []
    bin_labels = []
    bin_rates_for_label = []
    for i in range(args.nbins):
        t0, t1 = (bin_edges[i], bin_edges[i + 1])
        sel = (all_t >= t0) & (all_t < t1)
        img = make_image(all_x[sel], all_y[sel], xmin, xmax, ymin, ymax, binsize_phys)
        images.append(img)
        dt = t1 - t0
        n = sel.sum()
        rate = n / dt if dt > 0 else 0
        is_flare = rate > flare_thresh
        label = f"Bin {i}: {(t0 - tmin) / 1000:.0f}-{(t1 - tmin) / 1000:.0f} ks\n{n:,} ct, {rate:.0f} ct/s"
        if is_flare:
            label += " [FLARE]"
        bin_labels.append(label)
        bin_rates_for_label.append(rate)
    ncols = min(5, args.nbins)
    nrows = int(np.ceil(args.nbins / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes[np.newaxis, :]
    elif ncols == 1:
        axes = axes[:, np.newaxis]
    src_r_pix = args.src_r_arcsec * PHYS_PER_ARCSEC / binsize_phys
    cx_pix = (xcen - xmin) / binsize_phys
    cy_pix = (ycen - ymin) / binsize_phys
    quiet_rates = [r for r in bin_rates_for_label if r <= flare_thresh]
    if quiet_rates:
        ref_idx = [j for j, r in enumerate(bin_rates_for_label) if r <= flare_thresh]
        ref_img = np.median([images[j] for j in ref_idx], axis=0)
    else:
        ref_img = np.median(images, axis=0)
    vmax = np.percentile(ref_img[ref_img > 0], 99) if np.any(ref_img > 0) else 1
    for idx in range(nrows * ncols):
        row, col = divmod(idx, ncols)
        ax = axes[row, col]
        if idx >= args.nbins:
            ax.axis("off")
            continue
        img = images[idx]
        extent_arcsec = [
            -args.radius_arcsec,
            args.radius_arcsec,
            -args.radius_arcsec,
            args.radius_arcsec,
        ]
        ax.imshow(
            img,
            origin="lower",
            extent=extent_arcsec,
            cmap="magma",
            vmin=0,
            vmax=vmax,
            interpolation="nearest",
        )
        circ = plt.Circle(
            (0, 0), args.src_r_arcsec, ec="cyan", fc="none", lw=1, ls="--"
        )
        ax.add_patch(circ)
        is_flare = bin_rates_for_label[idx] > flare_thresh
        color = "red" if is_flare else "white"
        ax.set_title(bin_labels[idx], fontsize=7, color=color)
        ax.set_xlabel("East offset (arcsec)", fontsize=6)
        ax.set_ylabel("North offset (arcsec)", fontsize=6)
        ax.tick_params(labelsize=5)
    fig.suptitle(
        f'Quick comet-frame diagnostic — PI [{args.pimin}:{args.pimax}], {args.bin_arcsec}" pixels, {args.nbins} time bins',
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(str(outdir / "mosaic.png"), dpi=150)
    plt.close(fig)
    print(f"  -> {outdir / 'mosaic.png'}")
    print("\nGenerating full and quiet-only images...")
    full_img = make_image(all_x, all_y, xmin, xmax, ymin, ymax, binsize_phys)
    quiet_sel = np.zeros(len(all_t), dtype=bool)
    for i in range(args.nbins):
        if bin_rates_for_label[i] <= flare_thresh:
            t0, t1 = (bin_edges[i], bin_edges[i + 1])
            quiet_sel |= (all_t >= t0) & (all_t < t1)
    quiet_img = make_image(
        all_x[quiet_sel], all_y[quiet_sel], xmin, xmax, ymin, ymax, binsize_phys
    )
    n_quiet = quiet_sel.sum()
    n_total = len(all_t)
    print(f"  Quiet events: {n_quiet:,} / {n_total:,} ({100 * n_quiet / n_total:.0f}%)")
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    extent_arcsec = [
        -args.radius_arcsec,
        args.radius_arcsec,
        -args.radius_arcsec,
        args.radius_arcsec,
    ]
    for ax, img, title in [
        (ax1, full_img, f"All data ({n_total:,} events)"),
        (
            ax2,
            quiet_img,
            f"Quiet only ({n_quiet:,} events, {100 * n_quiet / n_total:.0f}%)",
        ),
        (
            ax3,
            full_img - quiet_img * (n_total / n_quiet) if n_quiet > 0 else full_img,
            "Flare excess (all − scaled quiet)",
        ),
    ]:
        v = np.percentile(img[img > 0], 99) if np.any(img > 0) else 1
        ax.imshow(
            img,
            origin="lower",
            extent=extent_arcsec,
            cmap="magma",
            vmin=0,
            vmax=v,
            interpolation="nearest",
        )
        circ = plt.Circle(
            (0, 0), args.src_r_arcsec, ec="cyan", fc="none", lw=1, ls="--"
        )
        ax.add_patch(circ)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("East offset (arcsec)")
        ax.set_ylabel("North offset (arcsec)")
    fig.suptitle(
        f"Comet-frame comparison — PI [{args.pimin}:{args.pimax}]", fontsize=11
    )
    fig.tight_layout()
    fig.savefig(str(outdir / "comparison.png"), dpi=150)
    plt.close(fig)
    print(f"  -> {outdir / 'comparison.png'}")
    print(f"\nAll outputs in: {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
