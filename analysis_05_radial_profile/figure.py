import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def make_figure(output, epic, pps):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), layout="constrained")
    for ax, name, label, fit in zip(
        axes,
        ("profile.tsv", "pps_profile.tsv"),
        ("Reduced EPIC", "PPS pn"),
        (epic["broken"], pps),
    ):
        with (output / name).open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        column = lambda key: np.array([float(row[key]) for row in rows])
        radius = column("r_arcmin")
        sky = column("sky_ct_s_arcmin2")
        ax.errorbar(
            radius,
            column("net_ct_s_arcmin2"),
            column("error_ct_s_arcmin2"),
            fmt=".",
            color="black",
            ms=4,
            lw=0.8,
            label="Counts minus particle background",
        )
        ax.plot(
            radius,
            column("comet_model_ct_s_arcmin2") + sky,
            label="Broken profile + sky",
        )
        net = column("net_ct_s_arcmin2")
        lower, upper = net[net > 0].min() * 0.45, net.max() * 1.5
        ax.set_ylim(lower, upper)
        if np.any(sky >= lower):
            ax.plot(radius, sky, ls=":", color="gray", label="Sky")
        ax.axvline(fit["break_radius_arcmin"], color="gray", ls="--", lw=0.8)
        ax.set(
            xscale="log",
            yscale="log",
            xlabel="Radius (arcmin)",
            ylabel="Count rate per arcmin²",
            title=f"{label}: 0.35–1.10 keV",
        )
        ax.text(
            0.04,
            0.06,
            f"α inner = {fit['alpha_inner']:.3f}\nα outer = {fit['alpha_outer']:.3f}",
            transform=ax.transAxes,
            fontsize=9,
        )
    axes[0].legend(fontsize=7, frameon=False)
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"profiles.{suffix}", dpi=180)
    plt.close(fig)
