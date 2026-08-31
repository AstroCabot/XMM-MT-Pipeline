import math

import numpy as np
from scipy.special import ellipeinc, ellipkinc

AU_M = 149_597_870_700.0
DAY_S = 86_400.0
K_B = 1.380_649e-23
U_KG = 1.660_539_066_60e-27
MASS_U = {"H2": 2.01588, "H2O": 18.01528, "CO2": 44.0095, "CO": 28.0101}
SPECIES = ("H2O", "CO2", "CO", "H2")


def q3(values):
    values = np.asarray(values, float)
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("invalid posterior draws")
    lo, mid, hi = np.quantile(values, (0.16, 0.5, 0.84))
    return {"p16": float(lo), "median": float(mid), "p84": float(hi)}


def fraction_envelope(momentum, forces):
    summaries = [q3(values / force) for values in momentum.values() for force in forces]
    medians = [item["median"] for item in summaries]
    return {
        "p16_min": min(item["p16"] for item in summaries),
        "median_range": [min(medians), max(medians)],
        "p84_max": max(item["p84"] for item in summaries),
    }


def thermal_speed(mass_kg, temperature_K):
    return math.sqrt(8 * K_B * temperature_K / (math.pi * mass_kg))


def ellipsoid_area(a, b, c):
    a, b, c = sorted((a, b, c), reverse=True)
    phi = math.acos(c / a)
    if phi == 0:
        return 4 * math.pi * a**2
    m = a**2 * (b**2 - c**2) / (b**2 * (a**2 - c**2))
    return 2 * math.pi * c**2 + 2 * math.pi * a * b / math.sin(phi) * (
        ellipeinc(phi, m) * math.sin(phi) ** 2 + ellipkinc(phi, m) * math.cos(phi) ** 2
    )


def shell_volume(radius_m, depth_m):
    if not 0 <= depth_m <= radius_m:
        raise ValueError("invalid shell depth")
    return 4 * math.pi / 3 * (radius_m**3 - (radius_m - depth_m) ** 3)


def delayed_radius(current_r_au, q_au, eccentricity, lookback_days):
    a = q_au / (1 - eccentricity)
    n = 0.01720209894846 * (-a) ** -1.5
    H = math.acosh((1 - current_r_au / a) / eccentricity)
    mean = eccentricity * math.sinh(H) - H - n * lookback_days
    H = math.asinh(mean / eccentricity)
    for _ in range(64):
        step = (eccentricity * math.sinh(H) - H - mean) / (
            eccentricity * math.cosh(H) - 1
        )
        H -= step
        if abs(step) < 1e-14:
            break
    return a * (1 - eccentricity * math.cosh(H))


def acceleration(A1_au_d2, activity_r_au):
    return A1_au_d2 * AU_M / DAY_S**2 / activity_r_au**2


def production_draws(draws, names):
    if draws.ndim != 2 or draws.shape[1] != len(names):
        raise ValueError("physical draw shape does not match names")
    missing = [name for name in SPECIES if f"Q_{name}" not in names]
    if missing:
        raise ValueError(f"missing production draws: {missing}")
    q = {name: draws[:, names.index(f"Q_{name}")] for name in SPECIES}
    valid = np.all(np.column_stack(list(q.values())) > 0, axis=1)
    if not np.all(valid):
        raise ValueError("invalid production draws")
    return {name: values[valid] for name, values in q.items()}, valid


def nucleus(parameters):
    p = parameters["nucleus"]
    radius_km = math.sqrt(p["area_product_km2"] / (math.pi * p["albedo"]))
    radius_sigma = 0.5 * p["area_product_sigma_km2"] / p["area_product_km2"] * radius_km
    radius_m = 1e3 * radius_km
    mass = 4 * math.pi / 3 * p["density_kg_m3"] * radius_m**3
    mass_sigma = 1.5 * p["area_product_sigma_km2"] / p["area_product_km2"] * mass
    return (
        radius_m,
        mass,
        {
            "area_product_km2": p["area_product_km2"],
            "albedo": p["albedo"],
            "radius_km": {"value": radius_km, "sigma": radius_sigma},
            "mass_kg": {"value": mass, "sigma": mass_sigma},
            "surface_area_m2": 4 * math.pi * radius_m**2,
            "albedo_exponents": {"radius": -0.5, "mass": -1.5},
        },
    )


def capacity_rates(radius_m, interval_days, depths, densities, ratios):
    water_mass = MASS_U["H2O"] * U_KG
    interval_s = interval_days * DAY_S
    return [
        shell_volume(radius_m, depth) * density / water_mass * ratio / interval_s
        for depth, density, ratio in zip(depths, densities, ratios, strict=True)
    ]


def calculate(draws, names, task08, parameters):
    q, selected = production_draws(draws, names)
    velocity_names = {name: f"v_{name}_km_s" for name in SPECIES[:3]}
    missing = [field for field in velocity_names.values() if field not in names]
    if missing:
        raise ValueError(f"missing velocity draws: {missing}")
    velocity = {
        name: draws[selected, names.index(field)]
        for name, field in velocity_names.items()
    }
    if any(
        np.any(~np.isfinite(values) | (values <= 0)) for values in velocity.values()
    ):
        raise ValueError("invalid velocity draws")
    masses = {name: MASS_U[name] * U_KG for name in SPECIES}
    radius_m, mass, nucleus_result = nucleus(parameters)
    r_xmm = float(task08["xmm_heliocentric_distance_au"])
    orbit = task08["orbit"]
    if float(task08["xmm_observation_jd_tdb"]) <= float(orbit["perihelion_jd_tdb"]):
        raise ValueError("Task 09 requires the post-perihelion XMM epoch")
    spada = parameters["spada"]
    r_delayed = delayed_radius(
        r_xmm, orbit["perihelion_q_au"], orbit["eccentricity"], spada["delay_days"]
    )
    a = {
        "symmetric": acceleration(spada["symmetric_A1_au_d2"], r_xmm),
        "delayed": acceleration(spada["delayed_A1_au_d2"], r_delayed),
    }
    force = {key: mass * value for key, value in a.items()}
    force_sigma = {
        key: nucleus_result["mass_kg"]["sigma"] * value for key, value in a.items()
    }
    parent_momentum = sum(
        q[name] * masses[name] * velocity[name] * 1e3 for name in SPECIES[:3]
    )
    h2_speed = {
        f"{int(T)}K": thermal_speed(masses["H2"], T)
        for T in parameters["gas"]["h2_temperature_K"]
    }
    h2_momentum = {
        key: q["H2"] * masses["H2"] * speed for key, speed in h2_speed.items()
    }
    total_momentum = {
        key: parent_momentum + value for key, value in h2_momentum.items()
    }
    q_summary = {name: q3(values) for name, values in q.items()}
    velocity_summary = {name: q3(values) for name, values in velocity.items()}

    oum = parameters["oumuamua"]
    semi = [axis / 2 for axis in oum["full_axes_m"]]
    oum_mass = 4 * math.pi / 3 * np.prod(semi) * oum["density_kg_m3"]
    oum_area = ellipsoid_area(*semi)
    atlas_area = 4 * math.pi * radius_m**2
    oum_speed = thermal_speed(masses["H2"], oum["temperature_K"])
    published_mass_loss = oum["published_H2_mass_loss_kg_s"]

    def oum_q(r_au):
        recoil_force = oum_mass * oum["acceleration_1au_m_s2"] / r_au**2
        return recoil_force / (oum["collimation"] * masses["H2"] * oum_speed)

    area_scaled = q["H2"] / (atlas_area / oum_area)

    reservoir = parameters["reservoir"]

    def capacity_case(depths):
        return {
            "depth_m": depths,
            "H2_molecules_s": capacity_rates(
                reservoir["radius_m"],
                reservoir["activity_interval_days"],
                depths,
                reservoir["water_equivalent_density_kg_m3"],
                reservoir["h2_per_h2o"],
            ),
        }

    return {
        "nucleus": nucleus_result,
        "radial_recoil": {
            "headline": "symmetric",
            "symmetric": {
                "A1_au_d2": spada["symmetric_A1_au_d2"],
                "delay_days": 0.0,
                "activity_r_au": r_xmm,
                "acceleration_m_s2": a["symmetric"],
                "required_force_N": {
                    "value": force["symmetric"],
                    "sigma": force_sigma["symmetric"],
                },
            },
            "delayed": {
                "A1_au_d2": spada["delayed_A1_au_d2"],
                "delay_days": spada["delay_days"],
                "activity_r_au": r_delayed,
                "acceleration_m_s2": a["delayed"],
                "required_force_N": {
                    "value": force["delayed"],
                    "sigma": force_sigma["delayed"],
                },
            },
        },
        "gas": {
            "Q_molecules_s": q_summary,
            "parent_speed_km_s": velocity_summary,
            "H2_thermal_speed_m_s": h2_speed,
            "H2_mass_loss_kg_s": q3(q["H2"] * masses["H2"]),
            "momentum_N": {
                "known_parents": q3(parent_momentum),
                **{f"H2_{key}": q3(value) for key, value in h2_momentum.items()},
                **{f"all_{key}": q3(value) for key, value in total_momentum.items()},
            },
            "supported_force_fraction_zeta1_central_mass": {
                "H2": fraction_envelope(h2_momentum, force.values()),
                "all_gas": fraction_envelope(total_momentum, force.values()),
            },
        },
        "oumuamua": {
            "mass_kg": oum_mass,
            "surface_area_m2": oum_area,
            "axes_albedo": oum["axes_albedo"],
            "temperature_K": oum["temperature_K"],
            "H2_thermal_speed_m_s": oum_speed,
            "collimation": oum["collimation"],
            "H2_requirement_molecules_s": {
                "at_1.25_au": oum_q(1.25),
                "at_xmm_distance": oum_q(r_xmm),
            },
            "published_H2_range_at_1.25_au": {
                "mass_loss_kg_s": published_mass_loss,
                "molecules_s": [value / masses["H2"] for value in published_mass_loss],
            },
            "atlas_to_oumuamua_area_ratio": atlas_area / oum_area,
            "atlas_H2_area_scaled_molecules_s": q3(area_scaled),
            "same_distance_rate_ratio": q3(area_scaled / oum_q(r_xmm)),
        },
        "radiolytic_capacity": {
            "radius_m": reservoir["radius_m"],
            "activity_interval_days": reservoir["activity_interval_days"],
            "water_equivalent_density_kg_m3": reservoir[
                "water_equivalent_density_kg_m3"
            ],
            "H2_per_H2O": reservoir["h2_per_h2o"],
            "gas_erosion": capacity_case(reservoir["gas_erosion_m"]),
            "dust_inclusive_erosion": capacity_case(reservoir["dust_erosion_m"]),
        },
    }
