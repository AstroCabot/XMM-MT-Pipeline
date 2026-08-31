import numpy as np
from astropy.time import Time

from solar_wind import AU_KM, GAP_HOURS, OMEGA_DEG_DAY, wrap180

WINDOW_SHADE = (0.72, 0.60, 0.30, 0.38)
PLOT_LIMITS_UTC = ("2025-12-01T23:00:00", "2025-12-06T01:00:00")
MAP_RADIAL_POINTS = 110
MAP_LONGITUDE_POINTS = 240
EMBER = (
    (0, "#000000"),
    (0.14, "#2b040a"),
    (0.34, "#7a1016"),
    (0.54, "#c9203a"),
    (0.72, "#f4610b"),
    (0.86, "#ffb627"),
    (0.95, "#ffe98a"),
    (1, "#ffffff"),
)
LINES = {
    "stereo_a": ("#e8590c", "-", 1.9, 5, "PLASTIC"),
    "solar_orbiter": ("#0b7285", "-", 1.4, 4, "PAS"),
    "wind": ("#1c3d8f", "-", 1.4, 3, "SWE (OMNI)"),
    "ace_scaled": ("#909090", "--", 1.3, 2, r"SWEPAM ($\times$1.7)"),
    "ace": ("#909090", ":", 1.3, 2, "SWEPAM (raw)"),
}
STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Nimbus Sans", "Arial", "DejaVu Sans"],
    "axes.labelsize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "figure.facecolor": "w",
    "savefig.facecolor": "w",
    "xtick.direction": "in",
    "ytick.direction": "in",
}


def nearest_source_epoch(mapped, requested):
    source = np.asarray(mapped["source"], float)
    epoch = float(source[np.argmin(abs(source - requested))])
    if abs(epoch - requested) > GAP_HOURS / 24:
        raise RuntimeError("Parker-map snapshot is outside the source support")
    return epoch


def parker_field(mapped, epoch):
    radius = np.linspace(0.1, 2, MAP_RADIAL_POINTS)
    longitude = np.linspace(0, 360, MAP_LONGITUDE_POINTS, endpoint=False)
    source, lon, r = mapped["source"], mapped["lon"], mapped["r"]
    if not source.min() <= epoch <= source.max():
        raise RuntimeError("Parker-map epoch is outside the source series")
    lon_ref = np.interp(epoch, source, lon)
    source_flux = np.convolve(
        np.pad(mapped["source_flux"], 1, mode="edge"), np.ones(3) / 3, mode="valid"
    )
    field = np.full((longitude.size, radius.size), np.nan)
    for i, phi in enumerate(longitude):
        delta = wrap180(phi - lon_ref) + lon_ref - lon
        for j, distance in enumerate(radius):
            arrival = (
                source
                + delta / OMEGA_DEG_DAY
                + (distance - r) * AU_KM / mapped["speed"] / 86400
            )
            order = np.argsort(arrival)
            if arrival[order][0] <= epoch <= arrival[order][-1]:
                field[i, j] = np.interp(
                    epoch,
                    arrival[order],
                    (source_flux * (r / distance) ** 2)[order],
                )
    return longitude, radius, field


def _hourly(grid, values):
    hours = np.floor(grid * 24).astype(np.int64)
    unique, index = np.unique(hours, return_inverse=True)
    output = np.full(unique.size, np.nan)
    for i in range(unique.size):
        sample = values[index == i]
        sample = sample[np.isfinite(sample)]
        if sample.size:
            output[i] = sample.mean()
    return (unique + 0.5) / 24, output


def _label(ax, letter, title):
    ax.annotate(
        letter,
        (0, 1),
        xycoords="axes fraction",
        xytext=(11, -10),
        textcoords="offset points",
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
    )
    ax.annotate(
        title,
        (1, 1),
        xycoords="axes fraction",
        xytext=(-9, -9),
        textcoords="offset points",
        ha="right",
        va="top",
        fontsize=12,
        fontweight="bold",
        multialignment="center",
        linespacing=1.3,
    )


def _map(fig, ax, cax, longitude, radius, field, epoch, tracks, mapped):
    from matplotlib.colors import LinearSegmentedColormap, LogNorm

    finite = field[np.isfinite(field) & (field > 0)]
    cmap = LinearSegmentedColormap.from_list("ember", EMBER)
    norm = LogNorm(
        max(float(np.percentile(finite, 1)), 1e5), float(np.percentile(finite, 99.7))
    )
    ax.set_facecolor("black")
    mesh = ax.pcolormesh(
        np.radians(longitude), radius, field.T, cmap=cmap, norm=norm, shading="auto"
    )
    comet_lon = np.interp(epoch, tracks["comet"]["mjd"], tracks["comet"]["lon"])
    comet_r = np.interp(epoch, tracks["comet"]["mjd"], tracks["comet"]["r"])
    speed = float(np.nanmedian(mapped["speed"]))
    radii = np.linspace(0.1, 2, 200)
    spiral = (
        np.radians(comet_lon)
        + np.radians(OMEGA_DEG_DAY) * (comet_r - radii) * AU_KM / speed / 86400
    )
    ax.plot(spiral, radii, color="white", lw=1, ls=":", alpha=0.85)
    for name, text, size in (("comet", "3I/ATLAS", 42), ("stereo_a", "STEREO-A", 18)):
        lon = np.interp(epoch, tracks[name]["mjd"], tracks[name]["lon"])
        distance = np.interp(epoch, tracks[name]["mjd"], tracks[name]["r"])
        ax.scatter(np.radians(lon), distance, s=size, color="white", zorder=6)
        ax.annotate(
            text,
            (np.radians(lon), distance),
            xytext=(12, -4),
            textcoords="offset points",
            color="white",
            fontsize=13,
            zorder=7,
        )
    ax.set_rmax(2)
    ax.set_rticks((0.5, 1, 1.5, 2))
    ax.set_yticklabels(("0.5", "1.0", "1.5", "2.0 AU"), fontsize=11)
    ax.set_rlabel_position(240)
    ax.grid(color="white", alpha=0.35, lw=0.5)
    ax.tick_params(axis="x", colors="black")
    ax.tick_params(axis="y", colors="white")
    bar = fig.colorbar(mesh, cax=cax)
    bar.set_label(r"Proton Flux (cm$^{-2}$ s$^{-1}$)")
    cax.yaxis.set_label_position("left")
    cax.yaxis.set_ticks_position("left")
    ax.annotate(
        "a",
        (0, 1),
        xycoords="axes fraction",
        xytext=(11, -10),
        textcoords="offset points",
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
    )


def _lightcurve_edges(midpoints):
    midpoints = np.asarray(midpoints)
    middle = 0.5 * (midpoints[1:] + midpoints[:-1])
    return np.r_[
        midpoints[0] - (middle[0] - midpoints[0]),
        middle,
        midpoints[-1] + (midpoints[-1] - middle[-1]),
    ]


def draw(wind, lightcurve, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    window, grid = wind["window"], wind["grid"]
    epoch = nearest_source_epoch(wind["mapped"]["stereo_a"], sum(window) / 2)
    longitude, radius, field = parker_field(wind["mapped"]["stereo_a"], epoch)
    shade = mdates.date2num(Time(window, format="mjd").to_datetime())
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(13.4, 5.8))
        polar = fig.add_axes([0.055, 0.100, 0.330, 0.795], projection="polar")
        cax = fig.add_axes([-0.004, 0.100, 0.016, 0.795])
        flux = fig.add_axes([0.435, 0.4975, 0.300, 0.3975])
        relative = fig.add_axes([0.435, 0.100, 0.300, 0.3975], sharex=flux)
        _map(
            fig,
            polar,
            cax,
            longitude,
            radius,
            field,
            epoch,
            wind["tracks"],
            wind["mapped"]["stereo_a"],
        )

        for name, (colour, style, width, zorder, text) in LINES.items():
            time, values = _hourly(grid, wind["curves"][name])
            dates = mdates.date2num(Time(time, format="mjd").to_datetime())
            flux.step(
                dates,
                values * 1e-7,
                where="mid",
                color=colour,
                ls=style,
                lw=width,
                zorder=zorder,
                label=text,
            )
            inside = (time >= window[0]) & (time <= window[1]) & np.isfinite(values)
            if name != "ace":
                relative.step(
                    dates,
                    values / np.nanmedian(values[inside]),
                    where="mid",
                    color=colour,
                    ls=style,
                    lw=width,
                    zorder=zorder,
                )
        for axis in (flux, relative):
            axis.axvspan(*shade, color=WINDOW_SHADE, lw=0)
            axis.yaxis.set_label_position("right")
            axis.yaxis.tick_right()
            axis.tick_params(top=True, right=True)
        map_date = mdates.date2num(Time(epoch, format="mjd").to_datetime())
        flux.axvline(map_date, color="#0b7285", ls="--", lw=1.2)
        flux.annotate(
            "Map",
            (map_date, 0.97),
            xytext=(-5, 0),
            textcoords="offset points",
            xycoords=("data", "axes fraction"),
            ha="right",
            va="top",
            color="#0b7285",
            fontsize=11,
        )
        flux.set_ylabel("Proton Flux\n" + r"(10$^7$ cm$^{-2}$ s$^{-1}$)")
        flux.tick_params(labelbottom=False)
        _label(flux, "b", "Ballistic\nMapping")

        midpoints = np.array(
            [Time(row["midpoint_utc"], format="isot").mjd for row in lightcurve]
        )
        rates = np.array([float(row["rate_ct_s"]) for row in lightcurve])
        edges = _lightcurve_edges(midpoints)
        xray = relative.stairs(
            rates / np.median(rates),
            mdates.date2num(Time(edges, format="mjd").to_datetime()),
            color="#111111",
            lw=2.4,
            zorder=6,
            baseline=None,
            label="EPIC-PN net rate",
        )
        relative.axhline(1, color=".5", lw=0.8, ls=":")
        relative.set_ylabel("Relative to\nWindow Median")
        relative.set_xlabel("Predicted Arrival Time at 3I/ATLAS (UTC)")
        relative.xaxis.set_major_locator(mdates.DayLocator())
        relative.xaxis.set_minor_locator(mdates.HourLocator(byhour=[12]))
        relative.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        limits = Time(PLOT_LIMITS_UTC, format="isot", scale="utc").to_datetime()
        relative.set_xlim(*mdates.date2num(limits))
        _label(relative, "c", "Relative\nVariability")
        handles, names = flux.get_legend_handles_labels()
        fig.legend(
            handles + [xray],
            names + ["EPIC-PN net rate"],
            ncol=3,
            loc="upper center",
            bbox_to_anchor=(0.585, 1.01),
            frameon=False,
        )
        for suffix in ("pdf", "png"):
            fig.savefig(
                output / f"figure4.{suffix}",
                dpi=300,
                bbox_inches="tight",
                pad_inches=0.05,
            )
        plt.close(fig)
