import numpy as np
from scipy.ndimage import map_coordinates

CONTOURS = (0.002, 0.004, 0.007, 0.012, 0.020)
SUN = "#ffb24a"
MOTION = "#5fd2ff"
NORTH = "#22aa66"
EAST = "#aa22cc"
STRIP = "#7cfc00"
REFLECTED = "#e0218a"
MIN_PROFILE_COVERAGE = 0.5
STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Nimbus Sans", "Arial", "DejaVu Sans"],
    "figure.facecolor": "w",
    "axes.facecolor": "w",
    "savefig.facecolor": "w",
    "text.color": "k",
    "axes.labelcolor": "k",
    "axes.edgecolor": "k",
    "xtick.color": "w",
    "ytick.color": "w",
    "xtick.labelcolor": "k",
    "ytick.labelcolor": "k",
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "axes.labelsize": 15,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "hatch.linewidth": 0.25,
}


def _outlined(artist):
    from matplotlib import patheffects

    artist.set_path_effects(
        [patheffects.Stroke(linewidth=2.6, foreground="black"), patheffects.Normal()]
    )
    return artist


def _panel_label(ax, letter, title, colour):
    for text, offset, size in ((letter, -6, 21), (title, -36, 14)):
        artist = ax.annotate(
            text,
            (0, 1),
            xycoords="axes fraction",
            xytext=(9, offset),
            textcoords="offset points",
            ha="left",
            va="top",
            color=colour,
            fontsize=size,
            fontweight="bold",
        )
        if colour == "white":
            _outlined(artist)


def _arrow(ax, tip, colour, label, fontsize=13):
    from matplotlib.patches import FancyArrowPatch

    arrow = FancyArrowPatch(
        (0, 0), tip, arrowstyle="->", mutation_scale=14, lw=1.8, color=colour, zorder=30
    )
    ax.add_patch(_outlined(arrow))
    _outlined(
        ax.text(
            tip[0] * 1.16,
            tip[1] * 1.16,
            label,
            color=colour,
            fontsize=fontsize,
            fontweight="bold",
            ha="center",
            va="center",
            zorder=31,
        )
    )


def _extent(x, y):
    hx, hy = 0.5 * np.diff(x[:2])[0], 0.5 * np.diff(y[:2])[0]
    return [x[0] - hx, x[-1] + hx, y[0] - hy, y[-1] + hy]


def _trim(smooth, covered, dx, dy):
    rows, cols = np.where(covered)
    r = slice(rows.min(), rows.max() + 1)
    c = slice(cols.min(), cols.max() + 1)
    image = np.where(covered[r, c], smooth[r, c] * 1e3, np.nan)
    return image, dx[r, c], dy[r, c]


def rotated(image, mask, dx, dy, sun_pa, half_width=12, aspect=1, scale=1):
    step = 0.1
    x = np.arange(-half_width / scale, half_width / scale + step, step)
    y = np.arange(
        -half_width * aspect / scale, half_width * aspect / scale + step, step
    )
    xx, yy = np.meshgrid(x, y)
    angle = np.deg2rad(sun_pa)
    east = -xx * np.sin(angle) - yy * np.cos(angle)
    north = -xx * np.cos(angle) + yy * np.sin(angle)
    pixel = np.median(abs(np.diff(dx[0])))
    cx = np.interp(0, dx[0, ::-1], np.arange(dx.shape[1])[::-1])
    cy = np.interp(0, dy[:, 0], np.arange(dy.shape[0]))
    col, row = cx - east / pixel, cy + north / pixel
    weight = map_coordinates(
        mask.astype(float), [row, col], order=1, mode="constant", cval=0
    )
    values = map_coordinates(
        np.nan_to_num(image) * mask, [row, col], order=1, mode="constant", cval=0
    )
    support = weight > 0.5
    values[support] /= weight[support]
    values[~support] = np.nan
    return values, support, x * scale, y * scale


def _map_panel(
    fig, ax, cax, smooth, covered, dx, dy, sun_pa, motion_pa, scale, vmin, vmax
):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    image, x, y = _trim(smooth, covered, dx, dy)
    x, y = x[0] * scale, y[:, 0] * scale
    extent = _extent(x, y)
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad(cmap(0))
    cmap.set_under(cmap(0))
    ax.set_facecolor(cmap(0))
    handle = ax.imshow(
        image,
        origin="lower",
        extent=extent,
        cmap=cmap,
        aspect="equal",
        interpolation="nearest",
        norm=LogNorm(vmin=vmin, vmax=vmax, clip=False),
    )
    limit = max(abs(value) for value in extent)
    ax.set_xlim(limit, -limit)
    ax.set_ylim(-limit, limit)
    levels = np.asarray(CONTOURS) * 1e3
    ax.contour(
        x,
        y,
        image,
        levels=levels[:3],
        colors="white",
        linewidths=0.7,
        alpha=0.6,
        zorder=25,
    )
    ax.contour(
        x,
        y,
        image,
        levels=levels[3:],
        colors="black",
        linewidths=0.7,
        alpha=0.9,
        zorder=25,
    )
    for xs, ys in (
        ((-limit, -0.02 * limit), (0, 0)),
        ((limit, 0.02 * limit), (0, 0)),
        ((0, 0), (-limit, -0.02 * limit)),
        ((0, 0), (0.02 * limit, limit)),
    ):
        ax.plot(xs, ys, color="white", ls="--", lw=0.5, alpha=0.7, zorder=20)
    ax.plot(0, 0, "+", color="white", ms=9, mew=1.4, zorder=26)
    length = 0.30 * limit
    for pa, colour, label in ((sun_pa, SUN, "Sun"), (motion_pa, MOTION, "v")):
        angle = np.deg2rad(pa)
        _arrow(ax, (length * np.sin(angle), length * np.cos(angle)), colour, label)
    _arrow(ax, (0, 0.18 * limit), NORTH, "N")
    _arrow(ax, (0.18 * limit, 0), EAST, "E")
    ax.grid(True, color="white", alpha=0.35, ls=":", lw=0.45)
    ax.tick_params(which="both", top=True, right=True)
    label = r"Projected distance ($10^5$ km)"
    ax.set(xlabel=label, ylabel=label)
    _panel_label(ax, "a", "Smoothed image", "white")
    bar = fig.colorbar(handle, cax=cax)
    bar.ax.yaxis.set_ticks([])
    bar.ax.minorticks_off()
    bar.outline.set_edgecolor("black")
    for level in levels:
        colour = "black" if level >= 12 else "white"
        bar.ax.plot(
            [0, 1],
            [level, level],
            color=colour,
            lw=1.3,
            transform=bar.ax.get_yaxis_transform(),
            zorder=10,
        )
        bar.ax.annotate(
            f"{level:g}",
            (0.5, level),
            ha="center",
            va="bottom",
            xytext=(0, 2),
            textcoords="offset points",
            xycoords=bar.ax.get_yaxis_transform(),
            fontsize=12,
            color=colour,
            annotation_clip=False,
        )


def _rotated_panel(
    ax,
    smooth,
    support,
    excluded,
    dx,
    dy,
    sun_pa,
    motion_pa,
    aperture,
    half_width,
    aspect,
    scale,
    vmin,
    vmax,
):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Circle, FancyArrowPatch, Rectangle

    image, visible, x, y = rotated(
        smooth * 1e3, support, dx, dy, sun_pa, half_width, aspect, scale
    )
    masked = (
        rotated(
            excluded.astype(float), support, dx, dy, sun_pa, half_width, aspect, scale
        )[0]
        > 0.5
    )
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad(cmap(0))
    cmap.set_under(cmap(0))
    ax.set_facecolor(cmap(0))
    extent = _extent(x, y)
    ax.imshow(
        np.clip(image, vmin, None),
        origin="lower",
        extent=extent,
        cmap=cmap,
        aspect="auto",
        interpolation="nearest",
        norm=LogNorm(vmin=vmin, vmax=vmax, clip=False),
    )
    masked &= visible
    if masked.any():
        hatch = ax.contourf(
            x,
            y,
            masked.astype(float),
            levels=[0.5, 1.5],
            colors="none",
            hatches=["///"],
            zorder=6,
        )
        for collection in getattr(hatch, "collections", [hatch]):
            collection.set_edgecolor((1, 1, 1, 0.65))
            collection.set_linewidth(0.25)

    strip = 0.5 * scale
    ax.add_patch(
        Rectangle(
            (extent[0], -strip),
            extent[1] - extent[0],
            2 * strip,
            fill=False,
            edgecolor=STRIP,
            lw=0.9,
            ls="--",
            zorder=20,
        )
    )
    label_x, stem = -0.70 * half_width, 0.6 * strip
    for sign in (1, -1):
        arrow = FancyArrowPatch(
            (label_x, sign * (strip + stem)),
            (label_x, sign * strip),
            arrowstyle="->",
            mutation_scale=11,
            lw=1.5,
            color=STRIP,
            zorder=21,
        )
        ax.add_patch(_outlined(arrow))
    _outlined(
        ax.text(
            label_x,
            strip + stem + 0.02 * half_width,
            "1′",
            color=STRIP,
            fontsize=12,
            ha="center",
            va="bottom",
            zorder=22,
        )
    )
    radius = aperture * scale
    ax.add_patch(
        Circle(
            (0, 0),
            radius,
            fill=False,
            edgecolor="white",
            lw=1,
            ls="-.",
            alpha=0.9,
            zorder=24,
        )
    )
    corner = radius / np.sqrt(2)
    ax.plot(
        [0, -corner], [0, -corner], color="white", lw=1.1, ls="--", alpha=0.9, zorder=24
    )
    _outlined(
        ax.text(
            -0.55 * corner,
            -0.55 * corner + 0.03 * half_width,
            "370″",
            color="white",
            fontsize=12,
            ha="center",
            va="bottom",
            rotation=45,
            rotation_mode="anchor",
            zorder=25,
        )
    )
    ax.plot(0, 0, "+", color="white", ms=9, mew=1.4, zorder=25)
    length = 0.26 * half_width
    for pa, colour, label in (
        (sun_pa, SUN, "Sun"),
        (motion_pa, MOTION, "v"),
        (0, NORTH, "N"),
        (90, EAST, "E"),
    ):
        angle = np.deg2rad(pa)
        along = np.sin(angle) * np.sin(np.deg2rad(sun_pa))
        along += np.cos(angle) * np.cos(np.deg2rad(sun_pa))
        across = np.sin(angle) * np.cos(np.deg2rad(sun_pa))
        across -= np.cos(angle) * np.sin(np.deg2rad(sun_pa))
        _arrow(ax, (-length * along, -length * across), colour, label, 12)
    ax.set_xlim(-half_width, half_width)
    ax.set_ylim(-half_width * aspect, half_width * aspect)
    ax.set_ylabel(r"Perpendicular ($10^5$ km)")
    ax.yaxis.set_label_position("right")
    ax.tick_params(
        which="both",
        top=True,
        right=True,
        labelbottom=False,
        labelleft=False,
        labelright=True,
    )
    _panel_label(ax, "b", "Sun-aligned", "white")


def _profile_panel(ax, profile, centroid, aperture, half_width, scale):
    table = np.asarray(profile)
    x = -table[:, 0] * scale
    use = (table[:, 4] >= MIN_PROFILE_COVERAGE) & np.isfinite(table[:, 2:4]).all(axis=1)
    use &= abs(x) <= half_width
    ax.errorbar(
        x[use],
        table[use, 2] * 1e3,
        yerr=table[use, 3] * 1e3,
        fmt="o",
        ms=2.6,
        color="k",
        lw=0,
        elinewidth=0.9,
        capsize=1.8,
        capthick=0.7,
        zorder=5,
        label="Data",
    )
    order = np.argsort(-x)
    reflected = np.where(use, table[:, 2] * 1e3, np.nan)
    ax.step(
        -x[order],
        reflected[order],
        where="mid",
        color=REFLECTED,
        lw=1.3,
        zorder=4,
        label="Reflected",
    )
    visible = table[use, 2] * 1e3
    for position, colour, style, label, side in (
        (0, "#0072b2", "-", "nucleus", 1),
        (-(centroid / 60) * scale, "#d55e00", "--", "centroid", -1),
    ):
        for ymin, ymax in ((0, 0.02), (0.082, 1)):
            ax.axvline(
                position,
                ymin=ymin,
                ymax=ymax,
                color=colour,
                lw=1.1,
                ls=style,
                alpha=0.9,
                zorder=3,
            )
        ax.annotate(
            label,
            (position, 0.42),
            ha="left" if side > 0 else "right",
            va="center",
            rotation=90,
            fontsize=11,
            color=colour,
            xycoords=("data", "axes fraction"),
            xytext=(5 * side, 0),
            textcoords="offset points",
        )
    radius = aperture * scale
    for sign in (-1, 1):
        ax.axvline(sign * radius, color=".6", lw=0.8, ls="-.", zorder=1)
    ax.annotate(
        "370″ aperture",
        (0, 0.03),
        xycoords=("data", "axes fraction"),
        ha="center",
        va="bottom",
        fontsize=11,
        color=".35",
    )
    for sign in (-1, 1):
        ax.annotate(
            "",
            (sign * (radius - 0.2), 0.047),
            xytext=(sign * 1.9, 0.047),
            xycoords=("data", "axes fraction"),
            textcoords=("data", "axes fraction"),
            arrowprops={"arrowstyle": "->", "color": ".45", "lw": 1.1},
        )
    ax.annotate(
        "",
        (0.04, 0.74),
        xytext=(0.24, 0.74),
        xycoords="axes fraction",
        arrowprops={"arrowstyle": "->", "color": ".2", "lw": 1.5},
    )
    ax.text(
        0.255,
        0.74,
        "towards Sun",
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=12,
        color=".2",
    )
    ax.set_yscale("log")
    ax.set_ylim(0.7 * np.percentile(visible, 5), 2.5 * visible.max())
    ax.set_xlim(-half_width, half_width)
    ax.set_xlabel(r"Projected distance along the comet-Sun line ($10^5$ km)")
    ax.set_ylabel(r"S.B. ($10^{-3}$ count s$^{-1}$ arcmin$^{-2}$)")
    ax.yaxis.set_label_position("right")
    ax.tick_params(
        which="both",
        top=True,
        right=True,
        labelleft=False,
        labelright=True,
        color="black",
    )
    ax.grid(True, which="major", color="k", alpha=0.15, ls=":", lw=0.5)
    handles, labels = ax.get_legend_handles_labels()
    order = [labels.index("Data"), labels.index("Reflected")]
    ax.legend(
        [handles[i] for i in order],
        [labels[i] for i in order],
        loc="upper right",
        framealpha=0.92,
    )
    _panel_label(ax, "c", "Slice profile", "black")


def draw(
    smooth,
    mask,
    dx,
    dy,
    profile,
    sun_pa,
    motion_pa,
    aperture_arcmin,
    output,
    km_per_arcmin=81920.5,
    centroid_arcsec=0,
    coverage=None,
):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    scale = km_per_arcmin / 1e5
    support = np.isfinite(smooth) & (smooth != 0)
    covered = support if coverage is None else support & np.asarray(coverage, bool)
    excluded = support & ~mask
    trimmed = _trim(smooth, covered, dx, dy)[0]
    positive = trimmed[np.isfinite(trimmed) & (trimmed > 0)]
    vmin, vmax = np.percentile(positive, 5), positive.max()
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(15, 7.6))
        grid = fig.add_gridspec(
            1,
            2,
            width_ratios=(1, 1),
            wspace=0,
            left=0.04,
            right=0.945,
            top=0.97,
            bottom=0.09,
        )
        left = fig.add_subplot(grid[0])
        right = grid[1].subgridspec(2, 1, height_ratios=(1, 1.05), hspace=0)
        upper = fig.add_subplot(right[0])
        lower = fig.add_subplot(right[1], sharex=upper)
        left_edge = upper.get_position().x0
        height = 0.97 - 0.09
        width = height * fig.get_figheight() / fig.get_figwidth()
        left.set_position([left_edge - 0.026 - width, 0.09, width, height])
        cax = fig.add_axes([left_edge - 0.026, 0.09, 0.026, height])
        _map_panel(
            fig,
            left,
            cax,
            smooth,
            covered,
            dx,
            dy,
            sun_pa,
            motion_pa,
            scale,
            vmin,
            vmax,
        )
        box = upper.get_position()
        aspect = box.height * fig.get_figheight() / (box.width * fig.get_figwidth())
        _rotated_panel(
            upper,
            smooth,
            support,
            excluded,
            dx,
            dy,
            sun_pa,
            motion_pa,
            aperture_arcmin,
            12,
            aspect,
            scale,
            vmin,
            vmax,
        )
        _profile_panel(lower, profile, centroid_arcsec, aperture_arcmin, 12, scale)
        fig.align_ylabels([upper, lower])
        fig.savefig(output / "figure1.pdf", dpi=300, bbox_inches="tight")
        fig.savefig(output / "figure1.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
