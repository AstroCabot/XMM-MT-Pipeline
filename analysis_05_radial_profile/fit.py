import math

import numpy as np
from scipy.optimize import brentq, minimize

# Half a 4-arcsec pixel regularizes the unresolved center.
CORE_RADIUS_ARCMIN = 1 / 30


def deviance(counts, model):
    model = np.maximum(model, 1e-30)
    positive = counts > 0
    terms = model - counts
    terms[positive] += counts[positive] * np.log(counts[positive] / model[positive])
    return 2 * float(terms.sum())


def radial_sum(sample, shape):
    size = len(sample.get("counts", sample.get("soft")))
    return np.bincount(sample["ring"], sample["weight"] * shape, minlength=size)


def power_shape(radius, alpha):
    return np.maximum(radius, CORE_RADIUS_ARCMIN) ** -alpha


def broken_shape(radius, alpha_in, alpha_out, radius_break):
    ratio = np.maximum(radius, CORE_RADIUS_ARCMIN) / radius_break
    return np.where(radius <= radius_break, ratio**-alpha_in, ratio**-alpha_out)


def best(objective, starts, bounds):
    fits = [
        minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 2000},
        )
        for start in starts
    ]
    finite = [fit for fit in fits if fit.success and np.isfinite(fit.fun)]
    if not finite:
        raise RuntimeError("profile fit failed")
    return min(finite, key=lambda fit: fit.fun)


def profile_upper(objective, point, bounds, index):
    free = [i for i in range(len(point)) if i != index]
    target = objective(point) + 1

    def delta(value):
        def reduced(values):
            trial = point.copy()
            trial[index], trial[free] = value, values
            return objective(trial)

        fit = minimize(
            reduced, point[free], method="L-BFGS-B", bounds=[bounds[i] for i in free]
        )
        if not fit.success:
            raise RuntimeError("background profile failed")
        return float(fit.fun - target)

    lower, upper = float(point[index]), bounds[index][1]
    trial = min(lower + 0.1, upper)
    while trial < upper and delta(trial) < 0:
        trial = min(lower + 2 * (trial - lower), upper)
    if delta(trial) < 0:
        raise RuntimeError("background interval is unbounded")
    crossing = brentq(delta, lower, trial)
    return math.exp(crossing) - math.exp(lower)


def fit_epic(sample, broken=False, background_interval=False):
    rate = (sample["counts"] - sample["fixed"]) / sample["area"]
    positive = rate[np.isfinite(rate) & (rate > 0)]
    level = float(np.median(positive))
    sky = max(float(np.percentile(positive, 10)), 1e-8)
    amplitude = max(level - sky, 1e-5)
    if broken:
        bounds = [(-30, 5), (-1, 1.9), (0, 3), (math.log(0.1), math.log(2)), (-30, 5)]

        def expected(theta):
            shape = broken_shape(
                sample["radius"], theta[1], theta[2], math.exp(theta[3])
            )
            return (
                math.exp(theta[0]) * radial_sum(sample, shape)
                + math.exp(theta[4]) * sample["area"]
                + sample["fixed"]
            )

        starts = [
            np.array([math.log(amplitude), 0.3, 1, math.log(radius), math.log(sky)])
            for radius in (0.3, 0.6, 1)
        ]
        background_index = 4
    else:
        bounds = [(-30, 5), (0, 3), (-30, 5)]

        def expected(theta):
            shape = power_shape(sample["radius"], theta[1])
            return (
                math.exp(theta[0]) * radial_sum(sample, shape)
                + math.exp(theta[2]) * sample["area"]
                + sample["fixed"]
            )

        starts = [
            np.array([math.log(amplitude), slope, math.log(sky)])
            for slope in (0.8, 1, 1.2)
        ]
        background_index = 2

    objective = lambda theta: deviance(sample["counts"], expected(theta))
    fit = best(objective, starts, bounds)
    error = (
        profile_upper(objective, fit.x, bounds, background_index)
        if background_interval
        else None
    )
    return {"theta": fit.x, "stat": float(fit.fun), "background_error": error}


def hard_expectation(soft, hard, vignetted, ratio):
    if ratio < 1e-8:
        return np.maximum(hard, 1e-30)
    b = vignetted * (1 + ratio) - ratio * (soft + hard)
    root = np.sqrt(b * b + 4 * ratio * (1 + ratio) * hard * vignetted)
    output = np.zeros_like(vignetted)
    positive = b >= 0
    np.divide(
        2 * hard * vignetted, b + root, out=output, where=positive & (b + root > 0)
    )
    output[~positive] = (-b[~positive] + root[~positive]) / (2 * ratio * (1 + ratio))
    return np.maximum(output, 1e-30)


def fit_pps(sample):
    rate = sample["soft"] / sample["area"]
    level = float(np.median(rate[rate > 0]))

    def components(theta):
        shape = power_shape(sample["radius"], theta[1])
        vignetted = (
            math.exp(theta[0]) * radial_sum(sample, shape)
            + math.exp(theta[2]) * sample["area"]
        )
        ratio = math.exp(theta[3])
        particle = hard_expectation(sample["soft"], sample["hard"], vignetted, ratio)
        return vignetted, particle, ratio

    def objective(theta):
        vignetted, particle, ratio = components(theta)
        return deviance(sample["soft"], vignetted + ratio * particle) + deviance(
            sample["hard"], particle
        )

    starts = [
        np.array(
            [math.log(level), slope, math.log(max(level / 20, 1e-8)), math.log(0.1)]
        )
        for slope in (0.8, 1, 1.2)
    ]
    fit = best(objective, starts, [(-30, 5), (0, 3), (-30, 5), (-20, 5)])
    _, particle, ratio = components(fit.x)
    return {
        "theta": fit.x,
        "stat": float(fit.fun),
        "hard_expectation": particle,
        "soft_to_hard": ratio,
    }


def integrated_broken(amplitude, alpha_in, alpha_out, radius_break, radius):
    def integral(lo, hi, alpha):
        if abs(alpha - 2) < 1e-8:
            return math.log(hi / lo)
        return (hi ** (2 - alpha) - lo ** (2 - alpha)) / (2 - alpha)

    if radius <= CORE_RADIUS_ARCMIN:
        return (
            math.pi
            * radius**2
            * amplitude
            * (CORE_RADIUS_ARCMIN / radius_break) ** -alpha_in
        )
    inner = (
        0.5 * CORE_RADIUS_ARCMIN**2 * (CORE_RADIUS_ARCMIN / radius_break) ** -alpha_in
    )
    inner += radius_break**alpha_in * integral(
        CORE_RADIUS_ARCMIN, min(radius, radius_break), alpha_in
    )
    outer = (
        0
        if radius <= radius_break
        else radius_break**alpha_out * integral(radius_break, radius, alpha_out)
    )
    return 2 * math.pi * amplitude * (inner + outer)
