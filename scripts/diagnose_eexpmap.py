#!/usr/bin/env python3
"""Diagnose eexpmap coordinate mismatch in comet-frame imaging.

Reads the moved ATTHK, comet event file WCS, and counts/exposure FITS
to show where eexpmap thinks the FOV is vs where the events actually are.

Usage:
    python3 scripts/diagnose_eexpmap.py output/
"""
import argparse
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="eexpmap coordinate diagnostic")
    parser.add_argument("workdir", type=Path)
    args = parser.parse_args()
    w = args.workdir
    outdir = w / "qc" / "eexpmap_diag"
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- 1. Read moved ATTHK ----
    atthk_path = w / "comet" / "moved_atthk.dat"
    with fits.open(atthk_path) as hdul:
        data = hdul[1].data
        ahf_ra = data["AHFRA"]
        ahf_dec = data["AHFDEC"]
        time_col = data["TIME"]
    print(f"Moved ATTHK: {len(ahf_ra)} rows")
    print(f"  AHFRA  range: {ahf_ra.min():.6f} .. {ahf_ra.max():.6f} deg")
    print(f"  AHFDEC range: {ahf_dec.min():.6f} .. {ahf_dec.max():.6f} deg")
    print(f"  AHFRA  median: {np.median(ahf_ra):.6f} deg")
    print(f"  AHFDEC median: {np.median(ahf_dec):.6f} deg")

    # ---- 2. Read ref_attitude ----
    ref_env = w / "track" / "ref_attitude.env"
    ref_ra = ref_dec = None
    for line in ref_env.read_text().splitlines():
        if line.startswith("COMET_REF_RA="):
            ref_ra = float(line.split("=", 1)[1].strip('"'))
        elif line.startswith("COMET_REF_DEC="):
            ref_dec = float(line.split("=", 1)[1].strip('"'))
    print(f"\nReference position: RA={ref_ra:.6f}, DEC={ref_dec:.6f}")

    # offset of moved boresight from REF, in arcsec
    dra = (ahf_ra - ref_ra) * 3600.0 * np.cos(np.radians(ref_dec))
    ddec = (ahf_dec - ref_dec) * 3600.0
    offset = np.sqrt(dra**2 + ddec**2)
    print(f"\nMoved boresight offset from REF:")
    print(f"  median: {np.median(offset):.1f} arcsec ({np.median(offset)/60:.1f} arcmin)")
    print(f"  max:    {offset.max():.1f} arcsec ({offset.max()/60:.1f} arcmin)")
    print(f"  min:    {offset.min():.1f} arcsec ({offset.min()/60:.1f} arcmin)")
    print(f"  dRA  range: {dra.min():.1f} .. {dra.max():.1f} arcsec")
    print(f"  dDEC range: {ddec.min():.1f} .. {ddec.max():.1f} arcsec")

    # ---- 3. Read comet event file WCS ----
    for inst in ["PN", "M1", "M2"]:
        evt_path = w / "comet" / f"{inst}_comet.fits"
        if not evt_path.exists():
            continue
        with fits.open(evt_path) as hdul:
            ehdr = hdul["EVENTS"].header
            xi = yi = None
            for i in range(1, 100):
                if ehdr.get(f"TTYPE{i}") == "X":
                    xi = i
                if ehdr.get(f"TTYPE{i}") == "Y":
                    yi = i
                if xi and yi:
                    break
            x_crpx = ehdr.get(f"TCRPX{xi}")
            x_crvl = ehdr.get(f"TCRVL{xi}")
            x_cdlt = ehdr.get(f"TCDLT{xi}")
            y_crpx = ehdr.get(f"TCRPX{yi}")
            y_crvl = ehdr.get(f"TCRVL{yi}")
            y_cdlt = ehdr.get(f"TCDLT{yi}")
            ra_pnt = ehdr.get("RA_PNT", ehdr.get("REFXCRVL"))
            dec_pnt = ehdr.get("DEC_PNT", ehdr.get("REFYCRVL"))
            ra_nom = ehdr.get("RA_NOM")
            dec_nom = ehdr.get("DEC_NOM")

            xvals = hdul["EVENTS"].data["X"]
            yvals = hdul["EVENTS"].data["Y"]

            print(f"\n{inst} event file WCS:")
            print(f"  TCRPX(X)={x_crpx}  TCRVL(X)={x_crvl}  TCDLT(X)={x_cdlt}")
            print(f"  TCRPX(Y)={y_crpx}  TCRVL(Y)={y_crvl}  TCDLT(Y)={y_cdlt}")
            print(f"  RA_PNT={ra_pnt}  DEC_PNT={dec_pnt}")
            print(f"  RA_NOM={ra_nom}  DEC_NOM={dec_nom}")
            print(f"  X range: {xvals.min():.0f} .. {xvals.max():.0f}, median={np.median(xvals):.0f}")
            print(f"  Y range: {yvals.min():.0f} .. {yvals.max():.0f}, median={np.median(yvals):.0f}")
            print(f"  Grid center X={43201}, Y={43201}")
        break  # just check first available

    # ---- 4. Read per-instrument and combined exposure maps ----
    for band in ["soft", "broad", "hard"]:
        band_dir = w / "images" / band
        if not band_dir.exists():
            continue
        print(f"\n--- {band} band ---")
        for tag in ["PN_counts", "M1_counts", "M2_counts", "EPIC_counts",
                     "PN_exp", "M1_exp", "M2_exp", "EPIC_exp",
                     "EPIC_rate", "EPIC_rate_masked"]:
            fpath = band_dir / f"{tag}.fits"
            if not fpath.exists():
                continue
            with fits.open(fpath) as hdul:
                img = hdul[0].data
                if img is None:
                    continue
                img = img.astype(float)
                nonzero = (img > 0).sum()
                frac = nonzero / img.size * 100
                total = img.sum()
                mx = img.max()
                mn = img[img > 0].min() if nonzero > 0 else 0
                print(f"  {tag}: shape={img.shape}  nonzero={frac:.1f}%  "
                      f"min(>0)={mn:.4g}  max={mx:.4g}  sum={total:.4g}")

    # ---- 5. Diagnostic plot ----
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Panel 1: Moved boresight trajectory relative to image center
    ax = axes[0, 0]
    sc = ax.scatter(dra, ddec, c=(time_col - time_col.min()) / 1000,
                    s=0.3, cmap="viridis", rasterized=True)
    ax.axhline(0, color="red", ls="--", lw=1)
    ax.axvline(0, color="red", ls="--", lw=1)
    circ_img = plt.Circle((0, 0), 900, fill=False, color="cyan", lw=2,
                           label="Image radius (900\")")
    ax.add_patch(circ_img)
    fov = plt.Circle((np.median(dra), np.median(ddec)), 900, fill=False,
                      color="orange", lw=2, ls="--",
                      label=f"FOV at median bore ({np.median(offset):.0f}\")")
    ax.add_patch(fov)
    lim = max(2500, offset.max() * 1.2)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("dRA from REF (arcsec)")
    ax.set_ylabel("dDEC from REF (arcsec)")
    ax.set_title("Moved ATTHK boresight vs image center")
    ax.legend(fontsize=8)
    plt.colorbar(sc, ax=ax, label="Time (ks)")

    # Panel 2: Boresight offset vs time
    ax = axes[0, 1]
    t_ks = (time_col - time_col.min()) / 1000.0
    ax.plot(t_ks, offset / 60.0, "k.", ms=0.5, rasterized=True)
    ax.axhline(15, color="red", ls="--", label="Image radius (15')")
    ax.set_xlabel("Time (ks)")
    ax.set_ylabel("Boresight offset from image center (arcmin)")
    ax.set_title("Boresight distance from image center vs time")
    ax.legend()

    # Panels 3-6: counts vs exposure for available bands
    panel_idx = 2
    for band in ["soft", "broad", "hard"]:
        if panel_idx >= 6:
            break
        counts_path = w / "images" / band / "EPIC_counts.fits"
        exp_path = w / "images" / band / "EPIC_exp.fits"
        if not counts_path.exists():
            continue
        with fits.open(counts_path) as hdul:
            cimg = hdul[0].data.astype(float)
            chdr = hdul[0].header
            cdelt1 = chdr.get("CDELT1", -0.000555556)

        ax = axes[panel_idx // 3, panel_idx % 3]
        vmax = max(np.percentile(cimg[cimg > 0], 99) if (cimg > 0).any() else 1, 1)
        ax.imshow(cimg, origin="lower", cmap="magma", vmin=0, vmax=vmax,
                  extent=[-900, 900, -900, 900])
        ax.set_title(f"{band} counts")
        ax.set_xlabel("arcsec")
        panel_idx += 1

        if exp_path.exists() and panel_idx < 6:
            with fits.open(exp_path) as hdul:
                eimg = hdul[0].data.astype(float)
            ax = axes[panel_idx // 3, panel_idx % 3]
            ax.imshow(eimg, origin="lower", cmap="viridis",
                      extent=[-900, 900, -900, 900])
            ax.set_title(f"{band} exposure (max={eimg.max():.0f}s)")
            ax.set_xlabel("arcsec")
            # Overlay where the median boresight is relative to image center
            ax.plot(np.median(dra), np.median(ddec), "r+", ms=15, mew=3,
                    label=f"Median boresight ({np.median(offset):.0f}\" off)")
            ax.legend(fontsize=8)
            panel_idx += 1

    fig.suptitle("eexpmap coordinate diagnostic — is the exposure map\n"
                 "aligned with the counts image?", fontsize=14)
    fig.tight_layout()
    out_png = outdir / "eexpmap_diagnostic.png"
    fig.savefig(out_png, dpi=150)
    print(f"\nDiagnostic plot saved to: {out_png}")
    plt.close(fig)


if __name__ == "__main__":
    main()
