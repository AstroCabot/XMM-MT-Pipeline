from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

from common import CFG, directory, finish, read_tsv, sky_wcs
from events import band_mask
from sources import event_wcs

STAGE = "08_figure"
PIXEL_ARCSEC = 4
STYLE = {
    "figure.facecolor": "w",
    "axes.facecolor": "w",
    "savefig.facecolor": "w",
    "text.color": "k",
    "axes.labelcolor": "k",
    "axes.titlecolor": "k",
    "axes.edgecolor": "k",
    "xtick.color": "w",
    "ytick.color": "w",
    "xtick.labelcolor": "k",
    "ytick.labelcolor": "k",
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Nimbus Sans", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def event_bounds(paths, center):
    trial = WCS(naxis=2)
    trial.wcs.crval = center
    trial.wcs.crpix = [1, 1]
    trial.wcs.cdelt = [-PIXEL_ARCSEC / 3600, PIXEL_ARCSEC / 3600]
    trial.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    bounds = [np.inf, np.inf, -np.inf, -np.inf]
    for path in paths:
        with fits.open(path, memmap=True) as hdul:
            data = hdul["EVENTS"].data
            sky = event_wcs(path)
            for start in range(0, len(data), 1_000_000):
                block = data[start : start + 1_000_000]
                good = np.isfinite(block["X"]) & np.isfinite(block["Y"])
                if not np.any(good):
                    continue
                ra, dec = sky.wcs_pix2world(block["X"][good], block["Y"][good], 1)
                px, py = trial.wcs_world2pix(ra, dec, 0)
                bounds = [
                    min(bounds[0], px.min()),
                    min(bounds[1], py.min()),
                    max(bounds[2], px.max()),
                    max(bounds[3], py.max()),
                ]
    if not np.all(np.isfinite(bounds)):
        raise RuntimeError("no finite PN events for the reduction mosaic")
    xmin, ymin = np.floor(bounds[:2]).astype(int) - 2
    xmax, ymax = np.ceil(bounds[2:]).astype(int) + 2
    trial.wcs.crpix = [1 - xmin, 1 - ymin]
    return trial, (int(ymax - ymin + 1), int(xmax - xmin + 1))


def event_image(paths, band, wcs, shape, quality):
    image = np.zeros(shape)
    selected = 0
    lo, hi = CFG["bands"][band]["pi"]
    for path in paths:
        sky = event_wcs(path)
        with fits.open(path, memmap=True) as hdul:
            data = hdul["EVENTS"].data
            for start in range(0, len(data), 1_000_000):
                block = data[start : start + 1_000_000]
                use = (
                    band_mask(block, "pn", band)
                    if quality
                    else (block["PI"] >= lo) & (block["PI"] <= hi)
                )
                use &= np.isfinite(block["X"]) & np.isfinite(block["Y"])
                selected += np.count_nonzero(use)
                ra, dec = sky.wcs_pix2world(block["X"][use], block["Y"][use], 1)
                x, y = wcs.wcs_world2pix(ra, dec, 0)
                image += np.histogram2d(
                    y,
                    x,
                    bins=shape,
                    range=((-0.5, shape[0] - 0.5), (-0.5, shape[1] - 0.5)),
                )[0]
    if image.sum() != selected:
        raise RuntimeError(f"event loss in reduction mosaic: {band}")
    return image


def exposure_image(paths, wcs, shape):
    y, x = np.indices(shape)
    ra, dec = wcs.wcs_pix2world(x, y, 0)
    total = np.zeros(shape)
    for path in paths:
        with fits.open(path, memmap=True) as hdul:
            hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
            data, source = np.asarray(hdu.data, float), sky_wcs(hdu.header)
        sx, sy = source.wcs_world2pix(ra, dec, 0)
        sampled = map_coordinates(
            data, [sy, sx], order=1, mode="constant", cval=0, prefilter=False
        )
        outside = (
            (sx < 0) | (sx > data.shape[1] - 1) | (sy < 0) | (sy > data.shape[0] - 1)
        )
        sampled[outside] = 0
        total += sampled
    return total


def _coordinates(axis, show_ra, show_dec):
    axis.coords[0].set_major_formatter("hh:mm:ss")
    axis.coords[1].set_major_formatter("dd:mm")
    axis.coords[0].set_ticks(number=3)
    axis.coords[1].set_ticks(number=3)
    for coordinate in (axis.coords[0], axis.coords[1]):
        coordinate.set_ticks_position("all")
        coordinate.display_minor_ticks(True)
        coordinate.set_ticks(color="white", size=4, width=0.8)
        coordinate.set_ticklabel(color="black")
    axis.coords[0].set_axislabel("RA (J2000)" if show_ra else "")
    axis.coords[1].set_axislabel("Dec (J2000)" if show_dec else "")
    axis.coords[0].set_ticklabel_visible(show_ra)
    axis.coords[1].set_ticklabel_visible(show_dec)
    axis.coords.grid(True, color="white", alpha=0.20, linestyle=":", linewidth=0.5)
    for spine in axis.spines.values():
        spine.set_edgecolor("black")


def draw(panels, wcs, output, sources=(), pointings=8):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize
    from matplotlib.patches import Circle

    stage_names = (
        ("final exposure map", "Exposure"),
        ("before GTI filtering", "Before espfilt"),
        ("after GTI filtering", "After espfilt"),
        ("after pointing and quality cuts", "After Clean+Cut"),
        ("source-mask overlay", "Source-mask Overlay"),
    )
    reference = panels["soft"]["after GTI filtering"]
    positive = reference[reference > 0]
    event_max = max(float(np.percentile(positive, 99.9)), 0.22)
    exposure_max = max(
        float(panels[band]["final exposure map"].max()) / 1000 for band in panels
    )
    ny, nx = reference.shape
    figure_width, figure_height = 11, 994.349 / 72
    left, right, bottom, top = 0.07, 0.88, 0.06, 0.99
    cell_width = (right - left) * figure_width / 2
    visible_height = nx * ((top - bottom) * figure_height / 5) / cell_width
    y0 = 0.5 * (ny - visible_height)
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(figure_width, figure_height))
        grid = fig.add_gridspec(
            5, 2, left=left, right=right, bottom=bottom, top=top, wspace=0, hspace=0
        )
        exposure_artist = event_artist = None
        for row, (stage, title) in enumerate(stage_names):
            for column, band in enumerate(("soft", "hard")):
                axis = fig.add_subplot(grid[row, column], projection=wcs)
                lo, hi = CFG["bands"][band]["pi"]
                if row == 0:
                    cmap = plt.get_cmap("bone").copy()
                    cmap.set_bad(cmap(0))
                    cmap.set_under(cmap(0))
                    axis.set_facecolor(cmap(0))
                    data = panels[band][stage] / 1000
                    exposure_artist = axis.imshow(
                        data,
                        origin="lower",
                        cmap=cmap,
                        norm=Normalize(0, exposure_max),
                        aspect="equal",
                        interpolation="nearest",
                    )
                    tag = f"{title}  —  {band.title()} ({lo}-{hi} eV)  [{pointings} pointings]"
                else:
                    cmap = plt.get_cmap("gist_stern").copy()
                    cmap.set_bad(cmap(0))
                    cmap.set_under(cmap(0))
                    axis.set_facecolor(cmap(0))
                    data = (
                        panels[band]["after pointing and quality cuts"]
                        if row == 4
                        else panels[band][stage]
                    )
                    event_artist = axis.imshow(
                        data,
                        origin="lower",
                        cmap=cmap,
                        norm=LogNorm(0.11, event_max, clip=False),
                        aspect="equal",
                        interpolation="nearest",
                    )
                    tag = f"{title}  —  {band.title()} ({lo}-{hi} eV)"
                    if row == 4:
                        pixels_per_arcsec = 1 / (abs(wcs.wcs.cdelt[0]) * 3600)
                        for source in sources:
                            x, y = wcs.world_to_pixel_values(
                                float(source["ra"]), float(source["dec"])
                            )
                            axis.add_patch(
                                Circle(
                                    (x, y),
                                    float(source["radius_arcsec"]) * pixels_per_arcsec,
                                    facecolor=(1, 1, 0, 0.30),
                                    edgecolor=(1, 1, 1, 0.95),
                                    linewidth=0.4,
                                )
                            )
                _coordinates(axis, row == 4, column == 0)
                axis.set_ylim(y0, y0 + visible_height)
                letter = chr(ord("a") + row * 2 + column)
                axis.text(
                    0.02,
                    0.97,
                    rf"$\bf{{{letter}}}$   {tag}",
                    transform=axis.transAxes,
                    ha="left",
                    va="top",
                    fontsize=13,
                    color="white",
                    bbox={
                        "facecolor": "black",
                        "alpha": 0.6,
                        "edgecolor": "none",
                        "pad": 2.5,
                    },
                )

        cell_height = (top - bottom) / 5
        exposure_bottom = top - cell_height
        exposure_axis = fig.add_axes([right, exposure_bottom, 0.018, cell_height])
        event_axis = fig.add_axes([right, bottom, 0.018, exposure_bottom - bottom])
        bars = (
            fig.colorbar(exposure_artist, cax=exposure_axis, label="Exposure (ks)"),
            fig.colorbar(event_artist, cax=event_axis, label="Counts/Bin"),
        )
        for bar in bars:
            bar.ax.tick_params(
                which="both", color="black", labelcolor="black", direction="out"
            )
            bar.outline.set_edgecolor("black")
        fig.savefig(output / "pn-sky-frame-mosaic.pdf")
        fig.savefig(output / "pn-sky-frame-mosaic.png", dpi=220)
        plt.close(fig)


def run():
    output = directory(STAGE)
    frames = [
        r for r in read_tsv(directory("03_events") / "frames.tsv") if r["det"] == "pn"
    ]
    sources = read_tsv(directory("04_sources") / "sources.tsv")
    background = read_tsv(directory("05_background") / "background.tsv")
    before = [directory("02_reprocess") / "events" / "pn-s003.fits"]
    after = [
        p
        for p in (directory("02_reprocess") / "filter" / "pn-s003").glob(
            "*-allevc.fits"
        )
        if "oot" not in p.name
    ]
    quality = [Path(r["event"]) for r in frames]
    center = (
        float(np.median([float(r["ra"]) for r in frames])),
        float(np.median([float(r["dec"]) for r in frames])),
    )
    wcs, shape = event_bounds(before + after + quality, center)
    panels = {}
    for band in ("soft", "hard"):
        exposure = [
            Path(r["exposure_map"])
            for r in background
            if r["det"] == "pn" and r["band"] == band
        ]
        panels[band] = {
            "final exposure map": exposure_image(exposure, wcs, shape),
            "before GTI filtering": event_image(before, band, wcs, shape, False),
            "after GTI filtering": event_image(after, band, wcs, shape, False),
            "after pointing and quality cuts": event_image(
                quality, band, wcs, shape, True
            ),
        }
    draw(panels, wcs, output, sources, len(frames))
    finish(STAGE)
