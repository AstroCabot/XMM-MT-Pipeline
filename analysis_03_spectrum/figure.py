import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np

CAMERAS = {
    "pn": ("#202020", "o", "PN"),
    "mos1": ("#d55e00", "s", "MOS1"),
    "mos2": ("#56b4e9", "^", "MOS2"),
}
COMPONENTS = {
    "cx": ("#c51b7d", "-", "SWCX"),
    "lhb": ("#0072b2", "--", "LHB"),
    "halo": ("#3a923a", ":", "halo"),
    "cxb": ("#7b2cbf", "-.", "CXB"),
}
REGIONS = {"inner": "0–370″ aperture", "outer": "400–600″ annulus"}
STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def draw(rows, png, pdf, band):
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(7.2, 5.2))
        grid = fig.add_gridspec(
            2,
            2,
            height_ratios=(3.15, 1),
            left=0.105,
            right=0.985,
            bottom=0.12,
            top=0.86,
            hspace=0,
            wspace=0,
        )
        top_left = fig.add_subplot(grid[0, 0])
        top_right = fig.add_subplot(grid[0, 1], sharex=top_left, sharey=top_left)
        lower_left = fig.add_subplot(grid[1, 0], sharex=top_left)
        lower_right = fig.add_subplot(grid[1, 1], sharex=top_left, sharey=lower_left)
        axes = ((top_left, lower_left), (top_right, lower_right))

        rate = np.array([row["rate_ct_s_keV"] for row in rows])
        model = np.array([row["model_ct_s_keV"] for row in rows])
        positive = np.r_[
            rate[np.isfinite(rate) & (rate > 0)],
            model[np.isfinite(model) & (model > 0)],
        ]
        y_limits = (max(1e-4, positive.min() * 0.35), positive.max() * 2.2)

        for column, region in enumerate(("inner", "outer")):
            top, residual = axes[column]
            for detector, (colour, marker, _label) in CAMERAS.items():
                group = [
                    row
                    for row in rows
                    if row["region"] == region and row["det"] == detector
                ]
                energy = np.array([row["energy_keV"] for row in group])
                rate = np.array([row["rate_ct_s_keV"] for row in group])
                error = np.array([row["rate_error_ct_s_keV"] for row in group])
                model = np.array([row["model_ct_s_keV"] for row in group])
                delta = np.array([row["residual_sigma"] for row in group])
                finite = (
                    np.isfinite(energy)
                    & np.isfinite(rate)
                    & np.isfinite(error)
                    & (rate > 0)
                )
                top.step(
                    energy,
                    model,
                    where="mid",
                    color=colour,
                    lw=1.2,
                    alpha=0.95,
                    zorder=2,
                )
                top.errorbar(
                    energy[finite],
                    rate[finite],
                    yerr=error[finite],
                    fmt=marker,
                    ms=2.25,
                    mew=0.35,
                    lw=0.45,
                    color=colour,
                    ecolor=colour,
                    alpha=0.82,
                    capsize=0,
                    zorder=3,
                    rasterized=False,
                )
                residual.step(
                    energy, delta, where="mid", color=colour, lw=0.8, alpha=0.95
                )
                if detector == "pn":
                    for name, (line_colour, style, _text) in COMPONENTS.items():
                        values = np.array([row[f"{name}_ct_s_keV"] for row in group])
                        top.step(
                            energy,
                            values,
                            where="mid",
                            color=line_colour,
                            ls=style,
                            lw=1.15 if name == "cx" else 1,
                            alpha=0.95,
                            zorder=1,
                        )

            top.set_yscale("log")
            top.set_ylim(*y_limits)
            for x, text, align, size in (
                (0.035, "ab"[column], "left", 13),
                (0.96, REGIONS[region], "right", 9.5),
            ):
                top.text(
                    x,
                    0.955,
                    text,
                    transform=top.transAxes,
                    ha=align,
                    va="top",
                    fontsize=size,
                    fontweight="bold",
                    bbox={
                        "facecolor": "white",
                        "edgecolor": "none",
                        "alpha": 0.78,
                        "pad": 1.4,
                    },
                )
            top.grid(True, which="major", color=".2", alpha=0.12, lw=0.45)
            top.tick_params(
                which="both",
                direction="in",
                top=True,
                right=True,
                labelbottom=False,
                labelleft=column == 0,
            )

            residual.axhline(0, color=".25", lw=0.65, zorder=0)
            residual.set_ylim(-5, 5)
            residual.set_yticks((-5, -2.5, 0, 2.5, 5))
            residual.grid(True, which="major", color=".2", alpha=0.12, lw=0.45)
            residual.tick_params(
                which="both",
                direction="in",
                top=True,
                right=True,
                labelleft=column == 0,
            )
            residual.set_xlim(*band)
            residual.xaxis.set_major_locator(MultipleLocator(0.2))
            residual.xaxis.set_minor_locator(AutoMinorLocator(4))
            residual.set_xlabel("Energy (keV)")
        top_left.set_ylabel(r"Count rate (counts s$^{-1}$ keV$^{-1}$)")
        lower_left.set_ylabel(r"$\Delta\chi$")

        camera_handles = [
            Line2D(
                [0],
                [0],
                color=colour,
                marker=marker,
                lw=1.2,
                ms=4,
                label=f"{label} data + model",
            )
            for colour, marker, label in CAMERAS.values()
        ]
        component_handles = [
            Line2D([0], [0], color=colour, ls=style, lw=1.2, label=label)
            for colour, style, label in COMPONENTS.values()
        ]
        spacer = Line2D([], [], ls="none", label="")
        handles = [
            camera_handles[0],
            component_handles[0],
            camera_handles[1],
            component_handles[1],
            camera_handles[2],
            component_handles[2],
            spacer,
            component_handles[3],
        ]
        legend = fig.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(0.105, 0.974),
            ncol=4,
            frameon=True,
            fancybox=False,
            framealpha=1,
            facecolor="white",
            edgecolor=".25",
            handlelength=2.5,
            columnspacing=1.35,
            handletextpad=0.55,
            labelspacing=0.55,
            borderpad=0.55,
            borderaxespad=0,
        )
        legend.get_frame().set_linewidth(0.6)
        fig.savefig(
            pdf,
            metadata={
                "Creator": "Task 03 fixed-model paper renderer",
                "Title": "Joint two-region EPIC spectral fit",
                "CreationDate": None,
                "ModDate": None,
            },
        )
        fig.savefig(
            png, dpi=300, metadata={"Software": "Task 03 fixed-model paper renderer"}
        )
        plt.close(fig)
