import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
XLIM = (0.09, 18)
SPECIES = (
    ("H2O", r"H$_2$O", "#00a6c7", "o"),
    ("CO", "CO", "#228833", "s"),
    ("CO2", r"CO$_2$", "#cc3377", "^"),
)
HASER = ((1e5, "#9db8d2"), (1e6, "#5d87b3"), (1e7, "#2c5985"))
STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Nimbus Sans", "DejaVu Sans"],
    "mathtext.fontset": "cm",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "text.color": "black",
    "axes.labelcolor": "black",
    "axes.edgecolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "axes.labelsize": 10.2,
    "xtick.labelsize": 9.2,
    "ytick.labelsize": 9.2,
    "legend.fontsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def table(path):
    path = Path(path)
    with path.open(newline="") as stream:
        rows = list(
            csv.DictReader(stream, delimiter="," if path.suffix == ".csv" else "\t")
        )
    output = {
        name: np.array([float(row[name]) for row in rows])
        for name in rows[0]
        if name != "species"
    }
    if "species" in rows[0]:
        output["species"] = np.array([row["species"] for row in rows])
    return output


def haser(radius, scale):
    node, weight = np.polynomial.legendre.leggauss(64)
    angle = np.pi / 4 * (node + 1)
    integral = (
        np.pi
        / 4
        * np.sum(weight * np.exp(-radius[:, None] / scale / np.cos(angle)), axis=1)
    )
    return integral / radius


def _log_value(x, y, at):
    return float(np.exp(np.interp(np.log(at), np.log(x), np.log(y))))


def _pow3(value, _position):
    if value <= 0:
        return ""
    exponent = int(np.floor(np.log10(value)))
    mantissa = value / 10**exponent
    if np.isclose(mantissa, 1):
        return rf"$10^{{{exponent}}}$"
    return rf"${mantissa:g}\times10^{{{exponent}}}$"


def _broken(radius, model):
    radius = np.asarray(radius)
    anchor = model["normalization_at_break_ct_s_arcmin2"]
    point = model["break_radius_arcmin"]
    slope = np.where(radius < point, model["alpha_inner"], model["alpha_outer"])
    return anchor * (radius / point) ** -slope


def _radial_image(radius, brightness, half_width=12):
    axis = np.linspace(-half_width, half_width, 241)
    xx, yy = np.meshgrid(axis, axis)
    rr = np.hypot(xx, yy)
    positive = np.maximum(brightness, np.finfo(float).tiny)
    clipped = np.clip(rr, radius[0], radius[-1])
    image = np.exp(np.interp(np.log(clipped), np.log(radius), np.log(positive)))
    image[rr > radius[-1]] = np.nan
    return image


def draw(
    radial_path,
    chemical_path,
    intrinsic_path,
    background,
    tau_radius,
    scale,
    km_per_arcmin,
    broken,
    output,
    ecf,
):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.ticker import FuncFormatter, LogLocator

    radial = table(radial_path)
    chemical = table(chemical_path)
    intrinsic = table(intrinsic_path)
    spherex = table(HERE / "spherex_profiles.csv")
    conversion = ecf / 1e-13
    radius = radial["r_arcmin"]
    data = (radial["net_ct_s_arcmin2"] - background) * conversion
    error = radial["error_ct_s_arcmin2"] * conversion
    chemical_data = chemical["data_ct_s_arcmin2"] * conversion
    chemical_error = chemical["total_error_ct_s_arcmin2"] * conversion
    chemical_model = chemical["prediction_ct_s_arcmin2"] * conversion
    lower = np.r_[data - error, chemical_data - chemical_error]
    upper = np.r_[data + error, chemical_data + chemical_error]
    lower = lower[np.isfinite(lower) & (lower > 0)]
    upper = upper[np.isfinite(upper) & (upper > 0)]
    y_limits = (0.9 * lower.min(), 1.1 * upper.max())
    label = "Surface Brightness\n" r"($10^{-13}$ erg cm$^{-2}$ s$^{-1}$ arcmin$^{-2}$)"
    point = broken["break_radius_arcmin"]
    anchor = broken["normalization_at_break_ct_s_arcmin2"] * conversion

    with plt.rc_context(STYLE):
        width, height = 3.75, 6.94
        left, right, top, bottom = 0.19, 0.97, 0.960, 0.045
        square = (right - left) * width / height
        top_bottom = top - square
        middle_bottom = top_bottom - square
        fig = plt.figure(figsize=(width, height))
        comparison = fig.add_axes([left, top_bottom, right - left, square])
        profile = fig.add_axes(
            [left, middle_bottom, right - left, square], sharex=comparison
        )
        residual = fig.add_axes(
            [left, bottom, right - left, middle_bottom - bottom], sharex=comparison
        )

        inner, outer = radius < point, radius >= point
        comparison.errorbar(
            radius[outer],
            data[outer],
            yerr=error[outer],
            fmt="s",
            ms=3.2,
            color="black",
            mfc="black",
            mec="black",
            mew=0.45,
            ecolor="black",
            elinewidth=0.7,
            capsize=0,
            lw=0,
            zorder=8,
            label="EPIC SWCX",
        )
        comparison.errorbar(
            radius[inner],
            data[inner],
            yerr=error[inner],
            fmt="s",
            ms=3.2,
            color=".35",
            mfc=".35",
            mec=".35",
            mew=0.45,
            ecolor=".35",
            elinewidth=0.7,
            capsize=0,
            lw=0,
            alpha=0.35,
            zorder=7,
        )
        model_radius = np.logspace(np.log10(XLIM[0]), np.log10(XLIM[1]), 500)
        model_curve = _broken(model_radius, broken) * conversion
        comparison.plot(
            model_radius[model_radius < point],
            model_curve[model_radius < point],
            color="black",
            lw=1.45,
            alpha=0.25,
            zorder=6,
        )
        comparison.plot(
            model_radius[model_radius >= point],
            model_curve[model_radius >= point],
            color="black",
            lw=1.45,
            zorder=6,
            label="Broken power-law fit",
        )

        for key, text, colour, marker in SPECIES:
            use = spherex["species"] == key
            r = spherex["rho_km_at_1.83au"][use] / km_per_arcmin
            value = spherex["SB_norm"][use]
            uncertainty = spherex["err_norm"][use]
            factor = anchor / _log_value(r, value, point)
            value, uncertainty = value * factor, uncertainty * factor
            low = np.maximum(value - uncertainty, np.finfo(float).tiny)
            high = value + uncertainty
            inside, outside = r < point, r >= point
            at_low = _log_value(r, low, point)
            at_high = _log_value(r, high, point)
            comparison.fill_between(
                np.r_[r[inside], point],
                np.r_[low[inside], at_low],
                np.r_[high[inside], at_high],
                color=colour,
                alpha=0.035,
                lw=0,
                zorder=2,
            )
            comparison.fill_between(
                np.r_[point, r[outside]],
                np.r_[at_low, low[outside]],
                np.r_[at_high, high[outside]],
                color=colour,
                alpha=0.11,
                lw=0,
                zorder=2,
            )
            comparison.plot(
                r[inside],
                value[inside],
                ls="none",
                marker=marker,
                ms=3.2,
                color=colour,
                alpha=0.25,
                mec=colour,
                mfc=colour,
                mew=0.45,
                zorder=5,
            )
            comparison.plot(
                r[outside],
                value[outside],
                ls="none",
                marker=marker,
                ms=3.2,
                color=colour,
                mec=colour,
                mfc=colour,
                mew=0.45,
                zorder=5,
                label=f"SPHEREx {text}",
            )

        guide_radius = np.logspace(np.log10(point), np.log10(XLIM[1]), 300)
        haser_handles = []
        haser_labels = []
        for length, colour in HASER:
            shape = haser(guide_radius, length / km_per_arcmin)
            value = anchor * shape / haser(np.array([point]), length / km_per_arcmin)[0]
            (line,) = comparison.plot(
                guide_radius, value, color=colour, ls=(0, (5, 2.6)), lw=0.85, zorder=2
            )
            haser_handles.append(line)
            haser_labels.append(rf"$L = 10^{{{int(np.log10(length))}}}\,\mathrm{{km}}$")
        slopes = (
            (1, ":", ".25"),
            (1.5, (0, (5, 2, 1, 2)), ".45"),
            (2, (0, (1, 1.8)), ".58"),
        )
        full_radius = np.logspace(np.log10(XLIM[0]), np.log10(XLIM[1]), 300)
        slope_handles = []
        for slope, style, colour in slopes:
            value = anchor * (full_radius / point) ** -slope
            comparison.plot(
                full_radius[full_radius < point],
                value[full_radius < point],
                color=colour,
                ls=style,
                lw=0.75,
                alpha=0.22,
                zorder=1,
            )
            (line,) = comparison.plot(
                full_radius[full_radius >= point],
                value[full_radius >= point],
                color=colour,
                ls=style,
                lw=0.75,
                alpha=0.78,
                zorder=1,
                label=rf"$\alpha = {slope:g}$",
            )
            slope_handles.append(line)
        comparison.axvline(point, color=".35", lw=0.8, ls=(0, (3, 2)), zorder=1)
        comparison.annotate(
            "fit break",
            (point, 0.32),
            xycoords=("data", "axes fraction"),
            xytext=(-4, 0),
            textcoords="offset points",
            rotation=90,
            ha="right",
            va="bottom",
            fontsize=7.3,
            color=".35",
        )
        comparison.set(
            xscale="log", yscale="log", xlim=XLIM, ylim=y_limits, ylabel=label
        )
        comparison.xaxis.set_major_locator(LogLocator(10, subs=(1, 3), numticks=20))
        comparison.xaxis.set_major_formatter(
            FuncFormatter(lambda value, _: f"{value:g}")
        )
        comparison.xaxis.set_minor_locator(
            LogLocator(10, subs=np.arange(2, 10) * 0.1, numticks=100)
        )
        comparison.yaxis.set_minor_locator(
            LogLocator(10, subs=np.arange(2, 10) * 0.1, numticks=100)
        )
        comparison.tick_params(
            which="both", top=False, right=True, bottom=False, labelbottom=False
        )
        comparison.tick_params(which="major", length=4)
        comparison.tick_params(which="minor", length=2.2)
        comparison.grid(which="major", color="0", alpha=0.20, ls=":", lw=0.6)
        comparison.grid(which="minor", color="0", alpha=0.075, ls=":", lw=0.4)
        secondary = comparison.secondary_xaxis(
            "top", functions=(lambda r: r * km_per_arcmin, lambda d: d / km_per_arcmin)
        )
        secondary.set_xlabel("Projected distance (km)")
        secondary.xaxis.set_major_locator(LogLocator(10, subs=(1, 3), numticks=20))
        secondary.xaxis.set_major_formatter(FuncFormatter(_pow3))
        secondary.tick_params(which="major", direction="in", length=4)
        secondary.tick_params(which="minor", direction="in", length=2.2)
        handles, names = comparison.get_legend_handles_labels()
        by_name = dict(zip(names, handles))
        order = [
            "EPIC SWCX",
            "Broken power-law fit",
            r"SPHEREx H$_2$O",
            "SPHEREx CO",
            r"SPHEREx CO$_2$",
        ]
        first = comparison.legend(
            [by_name[name] for name in order],
            order,
            loc="upper right",
            frameon=False,
            handlelength=2.35,
            labelspacing=0.35,
        )
        comparison.add_artist(first)
        comparison.legend(
            haser_handles + slope_handles,
            haser_labels + [rf"$\alpha = {s:g}$" for s, _, _ in slopes],
            loc="lower left",
            bbox_to_anchor=(0.012, 0.012),
            frameon=False,
            handlelength=2.1,
            labelspacing=0.20,
            fontsize=7.5,
        )

        edges = np.r_[chemical["r_lo_arcsec"], chemical["r_hi_arcsec"][-1]] / 60
        center = (edges[:-1] + edges[1:]) / 2
        xerr = np.vstack((center - edges[:-1], edges[1:] - center))
        intrinsic_radius = intrinsic["radius_arcmin"]
        intrinsic_flux = scale * intrinsic["energy_flux_erg_cm2_s_arcmin2"] / 1e-13
        (intrinsic_handle,) = profile.plot(
            intrinsic_radius,
            intrinsic_flux,
            color="#7B3F61",
            lw=1.2,
            ls="--",
            label="Intrinsic MAP solution",
        )
        psf_handle = profile.stairs(
            chemical_model,
            edges,
            baseline=None,
            color="#176D71",
            lw=1.5,
            label="PSF-convolved MAP",
        )
        annulus_handle = profile.errorbar(
            center,
            chemical_data,
            xerr=xerr,
            yerr=chemical_error,
            fmt="o",
            markersize=3.5,
            markerfacecolor="white",
            markeredgecolor="#303030",
            markeredgewidth=0.8,
            ecolor="#303030",
            elinewidth=0.7,
            capsize=0,
            label="EPIC SWCX annuli",
            zorder=5,
        )
        if tau_radius > 0:
            profile.axvline(tau_radius, color=".35", lw=0.8, ls=":")
            profile.text(
                tau_radius * 1.20,
                0.95,
                r"$\tau = 1$",
                transform=profile.get_xaxis_transform(),
                ha="left",
                va="top",
                color=".3",
                fontsize=7.7,
            )
        inset = profile.inset_axes([0.035, 0.035, 0.47, 0.47], zorder=8)
        image = _radial_image(intrinsic_radius, intrinsic_flux)
        positive = image[np.isfinite(image) & (image > 0)]
        inset.imshow(
            image,
            origin="lower",
            extent=(-12, 12, -12, 12),
            cmap="magma",
            norm=LogNorm(positive.min(), positive.max()),
            interpolation="bilinear",
            rasterized=True,
        )
        inset.set(xticks=[], yticks=[], xlim=(-12, 12), ylim=(-12, 12))
        for spine in inset.spines.values():
            spine.set_color(".25")
            spine.set_linewidth(0.6)
        inset.plot(
            [-9.36, -4.36], [-9.36, -9.36], color="white", lw=1.1, solid_capstyle="butt"
        )
        inset.text(
            -6.86, -8.61, r"$5'$", color="white", ha="center", va="bottom", fontsize=6.5
        )
        inset.text(
            -10.56,
            10.32,
            "2D Reconstruction\nof MAP model",
            color="white",
            ha="left",
            va="top",
            fontsize=6.5,
        )
        profile.set(yscale="log", ylim=y_limits, ylabel=label)
        profile.legend(
            [intrinsic_handle, psf_handle, annulus_handle],
            ["Intrinsic MAP solution", "PSF-convolved MAP", "EPIC SWCX annuli"],
            loc="upper right",
            frameon=False,
            handlelength=2.2,
            labelspacing=0.25,
            fontsize=7.6,
        )
        profile.tick_params(
            which="both",
            direction="in",
            top=True,
            right=True,
            bottom=False,
            labelbottom=False,
        )
        profile.tick_params(which="major", length=4)
        profile.tick_params(which="minor", length=2.2)

        standardized = (
            chemical["residual_ct_s_arcmin2"] / chemical["total_error_ct_s_arcmin2"]
        )
        residual.axhspan(-2, 2, color=".5", alpha=0.08, lw=0)
        residual.axhline(0, color=".25", lw=0.7)
        residual.axhline(2, color=".55", lw=0.6, ls=":")
        residual.axhline(-2, color=".55", lw=0.6, ls=":")
        residual.errorbar(
            center,
            standardized,
            xerr=xerr,
            yerr=np.ones_like(center),
            fmt="o",
            markersize=3.5,
            markerfacecolor="white",
            markeredgecolor="#303030",
            markeredgewidth=0.8,
            ecolor="#303030",
            elinewidth=0.7,
            capsize=0,
        )
        residual.set(
            ylim=(-3.8, 3.8),
            yticks=(-2, 0, 2),
            ylabel="Residual\n" + r"($\sigma$)",
            xlabel="Comet-Centric Radius (arcmin)",
        )
        residual.tick_params(which="both", direction="in", top=True, right=True)
        residual.tick_params(which="major", length=4)
        residual.tick_params(which="minor", length=2.2)
        for text, axis in zip("ab", (comparison, profile)):
            axis.text(
                0.025,
                0.94,
                text,
                transform=axis.transAxes,
                ha="left",
                va="top",
                fontsize=10.8,
                fontweight="bold",
                clip_on=True,
                zorder=10,
            )
        for axis in (comparison, profile, residual):
            for spine in axis.spines.values():
                spine.set_linewidth(0.7)
        comparison.spines["bottom"].set_visible(False)
        profile.spines["bottom"].set_visible(False)
        fig.align_ylabels((comparison, profile, residual))
        fig.savefig(Path(output) / "figure3.pdf", bbox_inches="tight", pad_inches=0.03)
        fig.savefig(
            Path(output) / "figure3.png", dpi=350, bbox_inches="tight", pad_inches=0.03
        )
        plt.close(fig)
