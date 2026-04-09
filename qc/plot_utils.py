#!/usr/bin/env python3
"""Shared plotting utilities for QC diagnostic pages.

Provides consistent formatting: matched colorbars, standardised norms,
clean axis labelling, and reusable figure/panel helpers.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import PowerNorm, Normalize, TwoSlopeNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
FIGDPI = 180
PAD_INCHES = 0.2
INSTS = ["PN", "M1", "M2"]
INST_COLORS = {"PN": "#2196F3", "M1": "#4CAF50", "M2": "#FF9800", "EPIC": "#9C27B0"}
NET_CMAP = "RdBu_r"
COUNTS_CMAP = "magma"
RATE_CMAP = "inferno"
EXPO_CMAP = "viridis"

# ---------------------------------------------------------------------------
# Grid / extent helpers
# ---------------------------------------------------------------------------

def load_grid(grid_path: str | Path) -> dict[str, Any]:
    """Load images/grid.json and return the dict."""
    return json.loads(Path(grid_path).read_text(encoding="utf-8"))


def grid_extent(grid: dict[str, Any]) -> tuple[float, float, float, float]:
    """Return (east_min, east_max, north_min, north_max) in arcsec from grid dict."""
    cx = float(grid["center_x_phys"])
    cy = float(grid["center_y_phys"])
    sc = float(grid["scale_arcsec_per_phys"])
    return (
        (float(grid["x_min_phys"]) - cx) * sc,
        (float(grid["x_max_phys"]) - cx) * sc,
        (float(grid["y_min_phys"]) - cy) * sc,
        (float(grid["y_max_phys"]) - cy) * sc,
    )


# ---------------------------------------------------------------------------
# Colorbar
# ---------------------------------------------------------------------------

def matched_colorbar(
    ax: plt.Axes,
    im: matplotlib.image.AxesImage,
    label: str = "",
    width: str = "4%",
    pad: float = 0.06,
) -> matplotlib.colorbar.Colorbar:
    """Add a colorbar that exactly matches the axes height.

    Uses ``make_axes_locatable`` so the colorbar is always the same height
    as the image, avoiding the ``shrink`` parameter which never quite matches.
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=width, pad=pad)
    cb = plt.colorbar(im, cax=cax)
    if label:
        cb.set_label(label, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    return cb


# ---------------------------------------------------------------------------
# Norm helpers
# ---------------------------------------------------------------------------

def rate_norm(
    data: np.ndarray,
    exposure: np.ndarray | None = None,
    gamma: float = 0.5,
    pct: float = 99.5,
    exp_frac: float = 0.1,
) -> PowerNorm:
    """PowerNorm sqrt-like stretch, limited to well-exposed pixels."""
    if exposure is not None and np.any(np.isfinite(exposure)):
        thresh = np.nanmax(exposure) * exp_frac
        good = np.isfinite(data) & (exposure > thresh)
    else:
        good = np.isfinite(data)
    if np.any(good):
        vmax = float(np.nanpercentile(data[good], pct))
        vmax = max(vmax, 1e-30)
    else:
        vmax = 1.0
    return PowerNorm(gamma=gamma, vmin=0.0, vmax=vmax)


def symmetric_norm(
    data: np.ndarray,
    exposure: np.ndarray | None = None,
    pct: float = 99.0,
    exp_frac: float = 0.1,
) -> Normalize:
    """Symmetric ±absmax linear norm for net-rate / residual images."""
    if exposure is not None and np.any(np.isfinite(exposure)):
        thresh = np.nanmax(exposure) * exp_frac
        good = np.isfinite(data) & (exposure > thresh)
    else:
        good = np.isfinite(data)
    if np.any(good):
        absmax = float(np.nanpercentile(np.abs(data[good]), pct))
        absmax = max(absmax, 1e-30)
    else:
        absmax = 1.0
    return Normalize(vmin=-absmax, vmax=absmax)


def sqrt_norm(data: np.ndarray, pct: float = 99.5) -> PowerNorm:
    """Sqrt stretch with percentile clipping."""
    finite = np.isfinite(data)
    if np.any(finite):
        vmax = float(np.nanpercentile(data[finite], pct))
        vmax = max(vmax, 1e-30)
    else:
        vmax = 1.0
    return PowerNorm(gamma=0.5, vmin=0.0, vmax=vmax)


# ---------------------------------------------------------------------------
# Figure creation
# ---------------------------------------------------------------------------

def make_figure(
    nrows: int = 1,
    ncols: int = 1,
    figscale: float = 4.0,
    left_pad: float = 0.7,
    right_pad: float = 0.5,
    suptitle: str = "",
    wspace: float = 0.35,
    hspace: float = 0.40,
) -> tuple[plt.Figure, np.ndarray]:
    """Create a figure with ``constrained_layout`` and uniform panel sizing.

    Returns (fig, axes) where *axes* always has shape ``(nrows, ncols)``.
    """
    w = left_pad + ncols * figscale + (ncols - 1) * wspace + right_pad + ncols * 0.5
    h = 0.5 + nrows * figscale + (nrows - 1) * hspace + (0.5 if suptitle else 0.0)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(w, h),
        constrained_layout=True,
        squeeze=False,
    )
    if suptitle:
        fig.suptitle(suptitle, fontsize=13, fontweight="bold")
    return fig, axes


# ---------------------------------------------------------------------------
# Image display
# ---------------------------------------------------------------------------

def imshow_arcsec(
    ax: plt.Axes,
    data: np.ndarray,
    extent: tuple[float, float, float, float],
    norm=None,
    cmap: str | None = None,
    title: str = "",
    cbar_label: str = "",
    mask_bad: bool = True,
    show_xlabel: bool = True,
    show_ylabel: bool = True,
) -> matplotlib.image.AxesImage:
    """Display a 2-D array with E/N arcsec extent and matched colorbar.

    Parameters
    ----------
    mask_bad : bool
        If True, grey-out NaN / zero-exposure pixels.
    """
    plot_data = np.where(np.isfinite(data), data, np.nan) if mask_bad else data
    if cmap is not None:
        cmap_obj = plt.get_cmap(cmap).copy()
    else:
        cmap_obj = plt.get_cmap().copy()
    cmap_obj.set_bad("0.85")  # light grey for NaN

    im = ax.imshow(plot_data, origin="lower", extent=extent, norm=norm, cmap=cmap_obj,
                   interpolation="nearest", aspect="equal")
    matched_colorbar(ax, im, label=cbar_label)
    if title:
        ax.set_title(title, fontsize=9, fontweight="bold")
    if show_xlabel:
        ax.set_xlabel("East offset (arcsec)", fontsize=8)
    else:
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelbottom=False)
    if show_ylabel:
        ax.set_ylabel("North offset (arcsec)", fontsize=8)
    else:
        ax.set_ylabel("")
        ax.tick_params(axis="y", labelleft=False)
    ax.tick_params(labelsize=7)
    return im


# ---------------------------------------------------------------------------
# Aperture overlays
# ---------------------------------------------------------------------------

def add_apertures(
    ax: plt.Axes,
    env: dict[str, str],
    color: str = "k",
    lw_src: float = 1.3,
    lw_bkg: float = 1.0,
) -> None:
    """Draw source circle + background annulus/circle on *ax*."""
    from matplotlib.patches import Circle

    def _ef(key: str, default: float) -> float:
        v = env.get(key, "")
        return float(v) if v else float(default)

    src_dx = _ef("SRC_DX_ARCSEC", 0.0)
    src_dy = _ef("SRC_DY_ARCSEC", 0.0)
    src_r = _ef("SRC_R_ARCSEC", 60.0)
    ax.add_patch(Circle((src_dx, src_dy), src_r, fill=False,
                         edgecolor=color, lw=lw_src, ls="-"))
    bkg_mode = env.get("BKG_MODE", "annulus").strip().lower()
    bkg_dx = _ef("BKG_DX_ARCSEC", 0.0)
    bkg_dy = _ef("BKG_DY_ARCSEC", 0.0)
    if bkg_mode == "circle":
        bkg_r = _ef("BKG_R_ARCSEC", 150.0)
        ax.add_patch(Circle((bkg_dx, bkg_dy), bkg_r, fill=False,
                             edgecolor=color, lw=lw_bkg, ls="--"))
    else:
        rin = _ef("BKG_RIN_ARCSEC", 90.0)
        rout = _ef("BKG_ROUT_ARCSEC", 150.0)
        ax.add_patch(Circle((bkg_dx, bkg_dy), rin, fill=False,
                             edgecolor=color, lw=lw_bkg, ls="--"))
        ax.add_patch(Circle((bkg_dx, bkg_dy), rout, fill=False,
                             edgecolor=color, lw=lw_bkg, ls="--"))


# ---------------------------------------------------------------------------
# Text panel
# ---------------------------------------------------------------------------

def summary_text_panel(ax: plt.Axes, lines: list[str], fontsize: float = 8.0) -> None:
    """Render key-value text lines in a panel with axis hidden."""
    ax.set_axis_off()
    text = "\n".join(lines)
    ax.text(0.05, 0.95, text, transform=ax.transAxes, va="top", ha="left",
            fontsize=fontsize, family="monospace",
            bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.9, pad=5))


# ---------------------------------------------------------------------------
# Scalebar
# ---------------------------------------------------------------------------

def add_scalebar(
    ax: plt.Axes,
    length_arcsec: float,
    label: str | None = None,
    loc: str = "lower right",
    color: str = "white",
    fontsize: float = 7,
    pad: float = 0.05,
) -> None:
    """Draw a horizontal scale bar inside the axes."""
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    dx = xlim[1] - xlim[0]
    dy = ylim[1] - ylim[0]
    bar_len = length_arcsec
    if label is None:
        if length_arcsec >= 60:
            label = f"{length_arcsec / 60:.0f}'"
        else:
            label = f'{length_arcsec:.0f}"'
    if "right" in loc:
        x0 = xlim[1] - pad * dx - bar_len
    else:
        x0 = xlim[0] + pad * dx
    if "lower" in loc:
        y0 = ylim[0] + pad * dy
    else:
        y0 = ylim[1] - pad * dy - 0.02 * dy
    ax.plot([x0, x0 + bar_len], [y0, y0], lw=2.5, color=color, solid_capstyle="butt")
    ax.text(x0 + 0.5 * bar_len, y0 + 0.015 * dy, label, ha="center", va="bottom",
            fontsize=fontsize, color=color, fontweight="bold")


# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

def savefig(fig: plt.Figure, path: str | Path) -> None:
    """Save figure at standard DPI and close."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out), dpi=FIGDPI, bbox_inches="tight", pad_inches=PAD_INCHES)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Mask-out low-exposure edges
# ---------------------------------------------------------------------------

def mask_low_exposure(
    data: np.ndarray,
    exposure: np.ndarray,
    frac: float = 0.05,
) -> np.ndarray:
    """Return a copy of *data* with pixels below *frac* * max(exposure) set to NaN."""
    out = np.array(data, dtype=float, copy=True)
    if exposure is None or not np.any(np.isfinite(exposure)):
        return out
    thresh = np.nanmax(exposure) * frac
    out[(~np.isfinite(exposure)) | (exposure < thresh)] = np.nan
    return out
