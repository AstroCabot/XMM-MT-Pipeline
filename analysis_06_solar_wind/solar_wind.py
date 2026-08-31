import csv
import datetime as dt
import os
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from astropy.time import Time

HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"
OUTPUT = HERE / "outputs"

AU_KM = 1.495978707e8
OMEGA_DEG_DAY = 14.1844
GRID_MINUTES = 10
SMOOTH_MINUTES = 120
GAP_HOURS = 1.5
EVAL_PAD_DAYS = 3.0
OBSERVATION_UTC = ("2025-12-03T12:04:07", "2025-12-04T10:11:54")
ACE_ANALYSIS_SCALE = 1.7
HAPI = "https://cdaweb.gsfc.nasa.gov/hapi/data"
HAPI_WINDOW = ("2025-11-08T00:00:00Z", "2025-12-30T00:00:00Z")
TRACK_WINDOW = ("2025-11-01", "2026-02-01")
PENTICTON = (
    "https://www.spaceweather.gc.ca/solar_flux_data/daily_flux_values/fluxtable.txt"
)
F107_DATES = ("2025-12-03", "2025-12-04")
F107_HALF_WINDOW = 40
MONITORS = {
    "stereo_a": (
        "STA_COHO1HR_MERGED_MAG_PLASMA",
        ("plasmaSpeed", "plasmaDensity"),
        1,
        0,
    ),
    "solar_orbiter": (
        "SOLO_COHO1HR_MERGED_MAG_PLASMA",
        ("ProtonSpeed", "protonDensity"),
        1,
        0,
    ),
    "ace": ("AC_K1_SWE", ("Np", "Vp"), 0, 1),
    "wind": ("OMNI2_H0_MRG1HR", ("PLS1800", "N1800", "V1800"), 1, 2),
}
SERIES = (*MONITORS, "ace_scaled")
TRACK_IDS = {
    "comet": "C/2025 N1",
    "stereo_a": "-234",
    "solar_orbiter": "-144",
    "ace": "-92",
    "wind": "399",
}


def write_tsv(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def fetch(url, path):
    if path.exists() and path.stat().st_size:
        return path.read_text()
    path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        url, headers={"User-Agent": "3I-ATLAS-reproduction"}
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        text = response.read().decode()
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text)
    temporary.replace(path)
    return text


def plasma(name):
    dataset, fields, density_col, speed_col = MONITORS[name]
    query = urllib.parse.urlencode(
        {
            "id": dataset,
            "parameters": ",".join(fields),
            "time.min": HAPI_WINDOW[0],
            "time.max": HAPI_WINDOW[1],
            "format": "csv",
        }
    )
    rows = []
    for row in csv.reader(fetch(f"{HAPI}?{query}", CACHE / f"{name}.csv").splitlines()):
        try:
            values = [float(value) for value in row[1:]]
            density, speed = values[density_col], values[speed_col]
            valid = 0.01 < density < 500 and 100 < speed < 1500
            if name == "wind":
                # OMNI code 52 is Wind definitive plasma.
                valid &= int(values[0]) == 52
            if valid:
                rows.append(
                    (Time(row[0], format="isot", scale="utc").mjd, density, speed)
                )
        except (ValueError, IndexError):
            continue
    data = np.asarray(rows)
    if len(data) < 10 or np.any(np.diff(data[:, 0]) <= 0):
        raise RuntimeError(f"invalid {dataset} series")
    segment = np.r_[0, np.cumsum(np.diff(data[:, 0]) > GAP_HOURS / 24)]
    return {
        "mjd": data[:, 0],
        "density": data[:, 1],
        "speed": data[:, 2],
        "segment": segment.astype(int),
    }


def track(name):
    path = CACHE / f"track_{name}.tsv"
    if not path.exists():
        os.environ["XDG_CACHE_HOME"] = str(CACHE / ".xdg")
        from astroquery.jplhorizons import Horizons

        table = Horizons(
            id=TRACK_IDS[name],
            location="@sun",
            id_type=None,
            epochs={"start": TRACK_WINDOW[0], "stop": TRACK_WINDOW[1], "step": "1h"},
        ).vectors(refplane="ecliptic", cache=False)
        x, y, z = (np.asarray(table[key], float) for key in ("x", "y", "z"))
        mjd = Time(
            np.asarray(table["datetime_jd"], float), format="jd", scale="tdb"
        ).utc.mjd
        write_tsv(
            path,
            ("mjd_utc", "longitude_deg", "radius_au"),
            (
                {
                    "mjd_utc": f"{t:.9f}",
                    "longitude_deg": f"{p:.9f}",
                    "radius_au": f"{r:.9f}",
                }
                for t, p, r in zip(
                    mjd,
                    np.degrees(np.unwrap(np.arctan2(y, x))),
                    np.sqrt(x * x + y * y + z * z),
                )
            ),
        )
    rows = read_tsv(path)
    return {
        "mjd": np.array([float(row["mjd_utc"]) for row in rows]),
        "lon": np.array([float(row["longitude_deg"]) for row in rows]),
        "r": np.array([float(row["radius_au"]) for row in rows]),
    }


def at(track_data, times):
    times = np.asarray(times)
    if times.min() < track_data["mjd"].min() or times.max() > track_data["mjd"].max():
        raise RuntimeError("Horizons track does not cover mapped times")
    return (
        np.interp(times, track_data["mjd"], track_data["lon"]),
        np.interp(times, track_data["mjd"], track_data["r"]),
    )


def wrap180(angle):
    return (np.asarray(angle) + 180) % 360 - 180


def map_samples(series, monitor, comet):
    source = series["mjd"]
    lon_m, r_m = at(monitor, source)
    arrival = source.copy()
    for _ in range(20):
        lon_c, r_c = at(comet, arrival)
        updated = (
            source
            + wrap180(lon_c - lon_m) / OMEGA_DEG_DAY
            + (r_c - r_m) * AU_KM / series["speed"] / 86400
        )
        converged = np.max(abs(updated - arrival)) < 1e-8
        arrival = updated
        if converged:
            break
    else:
        raise RuntimeError("ballistic mapping did not converge")
    source_flux = series["density"] * series["speed"] * 1e5
    return {
        "arrival": arrival,
        "flux": source_flux * (r_m / at(comet, arrival)[1]) ** 2,
        "speed": series["speed"],
        "source": source,
        "segment": series["segment"],
        "lon": lon_m,
        "r": r_m,
        "source_flux": source_flux,
    }


def smooth_segments(mapped, grid):
    """Smooth each segment over SMOOTH_MINUTES without bridging long gaps."""
    points = SMOOTH_MINUTES // GRID_MINUTES + 1
    half = points // 2
    kernel = np.r_[0.5, np.ones(points - 2), 0.5] / (points - 1)
    raw = np.full(grid.size, np.nan)
    smooth = np.full(grid.size, np.nan)
    for segment in np.unique(mapped["segment"]):
        use = mapped["segment"] == segment
        order = np.argsort(mapped["arrival"][use])
        times, unique = np.unique(mapped["arrival"][use][order], return_index=True)
        values = mapped["flux"][use][order][unique]
        inside = np.flatnonzero((grid >= times[0]) & (grid <= times[-1]))
        if not len(inside):
            continue
        raw[inside] = np.interp(grid[inside], times, values)
        if len(inside) >= points:
            smooth[inside[half:-half]] = np.convolve(raw[inside], kernel, mode="valid")
    return smooth


def f107():
    """Penticton 20:00 UT flux adjusted to 1 au, with the centred 81-day mean
    F10.7A and the Moore & Mapaye proxy F10.7P = (F10.7 + F10.7A)/2."""
    path = CACHE / "penticton_f107.tsv"
    if not path.exists():
        wanted = set()
        for date in F107_DATES:
            day = dt.date.fromisoformat(date)
            wanted |= {
                day + dt.timedelta(days=k)
                for k in range(-F107_HALF_WINDOW, F107_HALF_WINDOW + 1)
            }
        table = {}
        for line in fetch(PENTICTON, CACHE / "fluxtable.txt").splitlines()[2:]:
            parts = line.split()
            if len(parts) < 6 or parts[1] != "200000":
                continue
            try:
                day = dt.date(int(parts[0][:4]), int(parts[0][4:6]), int(parts[0][6:8]))
                if day in wanted:
                    table[day] = float(parts[5])
            except ValueError:
                continue
        write_tsv(
            path,
            ("date", "adjusted_flux_sfu"),
            (
                {"date": day.isoformat(), "adjusted_flux_sfu": f"{table[day]:.1f}"}
                for day in sorted(table)
            ),
        )
        (CACHE / "fluxtable.txt").unlink()
    daily = {
        dt.date.fromisoformat(row["date"]): float(row["adjusted_flux_sfu"])
        for row in read_tsv(path)
    }
    days = {}
    for date in F107_DATES:
        day = dt.date.fromisoformat(date)
        window = [
            daily.get(day + dt.timedelta(days=k))
            for k in range(-F107_HALF_WINDOW, F107_HALF_WINDOW + 1)
        ]
        if any(value is None for value in window):
            raise RuntimeError(f"incomplete 81-day F10.7 window for {date}")
        average = float(np.mean(window))
        days[date] = {
            "f107_sfu": daily[day],
            "f107a_sfu": round(average, 5),
            "f107p_sfu": round(0.5 * (daily[day] + average), 5),
        }
    return {
        "f107p_sfu": round(float(np.mean([d["f107p_sfu"] for d in days.values()])), 5),
        "units": "sfu",
        "daily": days,
        "convention": "Penticton 20:00 UT flux adjusted to 1 au; "
        "F10.7A is the centred 81-day mean; F10.7P = (F10.7 + F10.7A)/2",
    }


def observation_window():
    return tuple(
        Time(value, format="isot", scale="utc").mjd for value in OBSERVATION_UTC
    )


def ace_quality(series, mapped, window):
    use = (mapped["ace"]["arrival"] >= window[0]) & (
        mapped["ace"]["arrival"] <= window[1]
    )
    source = series["ace"]["mjd"][use]
    wind_density = np.interp(source, series["wind"]["mjd"], series["wind"]["density"])
    wind_speed = np.interp(source, series["wind"]["mjd"], series["wind"]["speed"])
    return {
        "dataset": MONITORS["ace"][0],
        "quality": "preliminary browse",
        "analysis_scale": ACE_ANALYSIS_SCALE,
        "source_interval_utc": [
            Time(source.min(), format="mjd").isot,
            Time(source.max(), format="mjd").isot,
        ],
        "source_samples": int(use.sum()),
        "density_ace_over_wind": float(
            np.mean(series["ace"]["density"][use] / wind_density)
        ),
        "speed_ace_over_wind": float(np.mean(series["ace"]["speed"][use] / wind_speed)),
    }


def products(output=OUTPUT):
    window = observation_window()
    bins = 1440 // GRID_MINUTES
    start = np.floor((window[0] - EVAL_PAD_DAYS) * bins) / bins
    stop = np.ceil((window[1] + EVAL_PAD_DAYS) * bins) / bins
    grid = np.arange(round((stop - start) * bins) + 1) / bins + start
    tracks = {name: track(name) for name in TRACK_IDS}
    series = {name: plasma(name) for name in MONITORS}
    mapped, curves = {}, {}
    for name in MONITORS:
        mapped[name] = map_samples(series[name], tracks[name], tracks["comet"])
        curves[name] = smooth_segments(mapped[name], grid)
    curves["ace_scaled"] = curves["ace"] * ACE_ANALYSIS_SCALE
    write_tsv(
        output / "mapped_flux.tsv",
        ("mjd_utc", *SERIES),
        (
            {
                "mjd_utc": f"{grid[i]:.9f}",
                **{
                    name: (
                        ""
                        if not np.isfinite(curves[name][i])
                        else f"{curves[name][i]:.8g}"
                    )
                    for name in SERIES
                },
            }
            for i in range(grid.size)
        ),
    )
    stats = {}
    for name, values in curves.items():
        sample = values[(grid >= window[0]) & (grid <= window[1]) & np.isfinite(values)]
        if not sample.size or sample.min() <= 0:
            raise RuntimeError(f"no positive mapped flux inside the window for {name}")
        stats[name] = {
            "median": float(np.median(sample)),
            "dynamic_range": float(sample.max() / sample.min()),
        }
    return {
        "window": window,
        "grid": grid,
        "tracks": tracks,
        "mapped": mapped,
        "curves": curves,
        "stats": stats,
        "ace_quality": ace_quality(series, mapped, window),
        "f107": f107(),
    }
