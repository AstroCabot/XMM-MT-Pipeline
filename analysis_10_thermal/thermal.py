import numpy as np
from scipy.linalg import solve_banded

AU_M = 1.496e11
MU_SUN = 1.32712440018e20
DAY_S = 86_400.0
SIGMA_SB = 5.670374e-8


def distance(q_au, eccentricity, days):
    a = q_au / (eccentricity - 1) * AU_M
    mean = np.sqrt(MU_SUN / a**3) * np.asarray(days, float) * DAY_S
    H = np.arcsinh(mean / eccentricity)
    for _ in range(50):
        H -= (eccentricity * np.sinh(H) - H - mean) / (eccentricity * np.cosh(H) - 1)
    return a * (eccentricity * np.cosh(H) - 1) / AU_M


def deepest(z, T, thresholds):
    return [z[T >= threshold].max(initial=0.0) for threshold in thresholds]


def simulate(conductivity, model, q_au, eccentricity, end_days, thresholds, depths=()):
    dz = model["step_m"]
    dt = model["step_days"]
    z = np.arange(0.0, model["depth_m"] + dz / 2, dz)
    absorbed = (1 - model["albedo"]) * model["solar_constant_W_m2"] / np.pi
    emissivity = model["emissivity"]
    T = np.full(
        z.size,
        (absorbed / (4 * emissivity * SIGMA_SB * model["initial_au"] ** 2)) ** 0.25,
    )
    T_max = T.copy()
    times = np.arange(model["start_days"], end_days + dt / 2, dt)
    r = distance(q_au, eccentricity, times)
    c = model["density_kg_m3"] * model["heat_capacity_J_kg_K"] * dz**2
    c /= conductivity * dt * DAY_S
    matrix = np.zeros((3, z.size))
    matrix[0, 2:] = -1.0
    matrix[1, 1:] = c + 2.0
    matrix[2, :-2] = -1.0
    matrix[2, -2] = -2.0
    link = conductivity / dz
    matrix[0, 1] = -link
    rows = [round(depth / dz) for depth in depths]
    days = np.arange(model["start_days"], end_days + 0.5)
    out = {"days": days, "T_K": [], "deepest_m": [], "reached_m": []}
    day = 0
    for k in range(1, times.size):
        rhs = c * T
        surface = T[0]
        for _ in range(3):
            emitted = emissivity * SIGMA_SB * surface**4
            slope = 4 * emitted / surface
            matrix[1, 0] = link + slope
            rhs[0] = absorbed / r[k] ** 2 - emitted + slope * surface
            new = solve_banded((1, 1), matrix, rhs)
            if abs(new[0] - surface) < 1e-6:
                break
            surface = new[0]
        T = new
        np.maximum(T_max, T, out=T_max)
        while day < days.size and times[k] >= days[day] - 1e-9:
            out["T_K"].append(T[rows])
            out["deepest_m"].append(deepest(z, T, thresholds))
            out["reached_m"].append(deepest(z, T_max, thresholds))
            day += 1
    return {key: np.array(value) for key, value in out.items()}
