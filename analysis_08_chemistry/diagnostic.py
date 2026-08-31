from pathlib import Path

import numpy as np

import inference
from figure import STYLE

PARAMETERS = (
    ("Q_H2O", True, 1, r"$\log_{10}(Q_{\rm H_2O}/{\rm s}^{-1})$"),
    ("Q_CO2", True, 1, r"$\log_{10}(Q_{\rm CO_2}/{\rm s}^{-1})$"),
    ("Q_CO", True, 1, r"$\log_{10}(Q_{\rm CO}/{\rm s}^{-1})$"),
    ("Q_H2", True, 1, r"$\log_{10}(Q_{\rm H_2}/{\rm s}^{-1})$"),
    ("beta_H2", False, 1, r"$\beta_{\rm H_2}$"),
    (
        "phi_heavy",
        True,
        1,
        r"$\log_{10}(\Phi_{\rm heavy}/{\rm cm}^{-2}\,{\rm s}^{-1})$",
    ),
    ("v_H2_km_s", False, 1, r"$v_{\rm H_2}\ ({\rm km\,s}^{-1})$"),
    ("scale", False, 1, r"$s$"),
    ("offset", False, 1e5, r"$10^5 b\ ({\rm count\,s}^{-1}\,{\rm arcmin}^{-2})$"),
    ("jitter", False, 1, r"$\sigma_{\rm jit}$"),
)
MODES = (
    ("Low H$_2$", "#E68613", lambda ratio: ratio <= -1.5),
    ("High H$_2$", "#008C95", lambda ratio: ratio >= -0.5),
)


def values(physical):
    physical = np.asarray(physical, float)
    columns = []
    for name, logarithmic, factor, _ in PARAMETERS:
        column = factor * physical[:, inference.PHYSICAL_NAMES.index(name)]
        columns.append(np.log10(column) if logarithmic else column)
    return np.column_stack(columns)


def draw(physical, log_probability, path, limit=4000):
    import corner
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    matrix = values(physical)
    physical = np.asarray(physical, float)
    log_probability = np.asarray(log_probability, float)
    if len(log_probability) != len(matrix):
        raise ValueError("posterior coordinates and probabilities differ")

    ratio = physical[:, inference.PHYSICAL_NAMES.index("log10_H2_ratio")]
    families = []
    for label, color, criterion in MODES:
        indices = np.flatnonzero(criterion(ratio))
        if len(indices):
            families.append((label, color, indices))
    if not families:
        raise RuntimeError("no declared H2 mode in posterior draws")
    count = min(limit, *(len(indices) for _, _, indices in families))
    selected = []
    for label, color, indices in families:
        take = indices[np.linspace(0, len(indices) - 1, count, dtype=int)]
        best = matrix[indices[np.argmax(log_probability[indices])]]
        selected.append((label, color, matrix[take], best))

    combined = np.vstack([rows for _, _, rows, _ in selected])
    bounds = []
    for column in combined.T:
        lower, upper = np.quantile(column, (0.005, 0.995))
        pad = 0.05 * max(upper - lower, 1e-12)
        bounds.append((lower - pad, upper + pad))

    with plt.rc_context(STYLE):
        corner_figure = None
        labels = [parameter[-1] for parameter in PARAMETERS]
        for _, color, rows, _ in selected:
            corner_figure = corner.corner(
                rows,
                fig=corner_figure,
                labels=labels,
                range=bounds,
                bins=22,
                smooth=0.9,
                color=color,
                levels=(0.393, 0.865),
                plot_datapoints=False,
                plot_density=False,
                fill_contours=False,
                max_n_ticks=3,
                labelpad=0.18,
                label_kwargs={"fontsize": 8},
                hist_kwargs={"density": True, "linewidth": 1.1},
                contour_kwargs={"linewidths": 1.0},
            )
        corner_figure.set_size_inches(13, 13)
        axes = np.asarray(corner_figure.axes).reshape(len(labels), len(labels))
        for _, color, _, best in selected:
            for row in range(len(labels)):
                axes[row, row].axvline(best[row], color=color, lw=0.9)
                for column in range(row):
                    axes[row, column].plot(
                        best[column], best[row], marker="x", color=color, ms=4, mew=1
                    )
        corner_figure.legend(
            [
                Line2D([0], [0], color=color, lw=1.5, marker="x", markersize=4)
                for _, color, _, _ in selected
            ],
            [label for label, _, _, _ in selected],
            loc="upper right",
            frameon=False,
            fontsize=9,
        )
        corner_figure.text(
            0.98,
            0.955,
            "39/87% conditional contours; transition family omitted\n"
            "families normalized separately; occupancy is not probability",
            ha="right",
            va="top",
            fontsize=8,
            color="0.35",
        )
        corner_figure.savefig(
            Path(path), dpi=220, bbox_inches="tight", facecolor="white"
        )
        plt.close(corner_figure)
