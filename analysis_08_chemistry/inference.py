import math
import multiprocessing as mp

import numpy as np
from scipy.special import ndtr

COORDS = (
    "ln_Qref_H2O",
    "ln_Qref_CO2",
    "ln_Qref_CO",
    "beta_H2O",
    "beta_CO2",
    "beta_CO",
    "beta_H2",
    "log10_H2_ratio",
    "phi_heavy",
    "v_H2O",
    "v_CO2",
    "v_CO",
    "v_H2",
    "v_OH",
    "v_O_H2O",
    "v_O_OH",
    "v_H",
    "ln_scale",
    "offset",
    "jitter",
)
PHYSICAL_NAMES = (
    "Q_H2O",
    "Q_CO2",
    "Q_CO",
    "Q_H2",
    "beta_H2O",
    "beta_CO2",
    "beta_CO",
    "beta_H2",
    "log10_H2_ratio",
    "phi_heavy",
    "v_H2O_km_s",
    "v_CO2_km_s",
    "v_CO_km_s",
    "v_H2_km_s",
    "v_OH_km_s",
    "v_O_H2O_km_s",
    "v_O_OH_km_s",
    "v_H_km_s",
    "scale",
    "offset",
    "jitter",
)


def normal_logpdf(value, mean, sigma):
    return -0.5 * ((value - mean) / sigma) ** 2 - math.log(
        sigma * math.sqrt(2 * math.pi)
    )


def positive_normal_logpdf(value, mean, sigma):
    if value <= 0:
        return -math.inf
    return normal_logpdf(value, mean, sigma) - math.log(ndtr(mean / sigma))


def uniform_logpdf(value, lower, upper):
    return -math.log(upper - lower) if lower <= value <= upper else -math.inf


def halfnormal_logpdf(value, sigma):
    if value < 0:
        return -math.inf
    return math.log(2) + normal_logpdf(value, 0, sigma)


def unpack(theta, config, orbit):
    value = dict(zip(COORDS, np.asarray(theta, float), strict=True))
    beta = {name: value[f"beta_{name}"] for name in ("H2O", "CO2", "CO", "H2")}
    reference = config["production_reference_r_au"]
    q = {
        name: math.exp(value[f"ln_Qref_{name}"])
        * (orbit.current_r_au / reference) ** (-beta[name])
        for name in ("H2O", "CO2", "CO")
    }
    q["H2"] = 10 ** value["log10_H2_ratio"] * sum(q.values())
    speed = {name: value[f"v_{name}"] for name in ("H2O", "CO2", "CO", "H2")}
    speed.update(
        OH=value["v_OH"],
        O_H2O=value["v_O_H2O"],
        O_OH=value["v_O_OH"],
        H=value["v_H"],
        heavy=config["fixed_velocity_km_s"]["heavy_daughter"],
        H2_photo=config["fixed_velocity_km_s"]["photolytic_H2"],
    )
    physical = [q[name] for name in ("H2O", "CO2", "CO", "H2")]
    physical += [beta[name] for name in ("H2O", "CO2", "CO", "H2")]
    physical += [value[name] for name in COORDS[7:17]]
    physical += [math.exp(value["ln_scale"]), value["offset"], value["jitter"]]
    return value, q, beta, speed, np.asarray(physical)


class Posterior:
    def __init__(
        self,
        data,
        error,
        ecf,
        transport,
        chains,
        psf_matrix,
        config,
        photon_energy_ev,
        scale_sigma,
        offset_sigma,
        heavy_location,
        heavy_scale,
    ):
        self.data, self.error = np.asarray(data), np.asarray(error)
        self.ecf, self.transport, self.chains = ecf, transport, chains
        self.psf = np.asarray(psf_matrix)
        self.config = {
            **config,
            "fixed_velocity_km_s": config["fixed_velocity_km_s"].copy(),
        }
        self.photon_energy_ev = photon_energy_ev
        self.scale_sigma, self.offset_sigma = scale_sigma, offset_sigma
        self.heavy_location, self.heavy_scale = heavy_location, heavy_scale
        if self.data.shape != self.error.shape or self.psf.shape[0] != len(self.data):
            raise ValueError("profile and PSF dimensions differ")

    def log_prior(self, theta):
        v = dict(zip(COORDS, np.asarray(theta, float), strict=True))
        p = self.config["priors"]
        result = sum(
            normal_logpdf(
                v[f"ln_Qref_{name}"],
                math.log(p[f"Qref_{name}"][0]),
                p[f"Qref_{name}"][1],
            )
            for name in ("H2O", "CO2", "CO")
        )
        for name in ("H2O", "CO2", "CO"):
            result += normal_logpdf(v[f"beta_{name}"], *p[f"beta_{name}"])
        result += uniform_logpdf(v["beta_H2"], *p["beta_H2"])
        result += uniform_logpdf(v["log10_H2_ratio"], *p["log10_H2_ratio"])
        result += positive_normal_logpdf(
            v["phi_heavy"], self.heavy_location, self.heavy_scale
        )
        result += sum(
            positive_normal_logpdf(v[f"v_{name}"], *p["parent_velocity"])
            for name in ("H2O", "CO2", "CO", "H2")
        )
        result += uniform_logpdf(v["v_OH"], *p["OH_velocity"])
        result += uniform_logpdf(v["v_O_H2O"], *p["O_velocity"])
        result += uniform_logpdf(v["v_O_OH"], *p["O_velocity"])
        result += uniform_logpdf(v["v_H"], *p["H_velocity"])
        result += normal_logpdf(v["ln_scale"], 0, self.scale_sigma)
        result += normal_logpdf(v["offset"], 0, self.offset_sigma)
        result += halfnormal_logpdf(v["jitter"], p["jitter_sigma"])
        return result

    def predict(self, theta):
        value, q, beta, speed, physical = unpack(
            theta, self.config, self.transport.orbit
        )
        source, tau = self.transport.source_profile(
            q,
            beta,
            speed,
            self.chains,
            value["phi_heavy"],
            self.photon_energy_ev,
            self.config["photons_per_captured_ion"],
        )
        model = self.psf @ source / self.ecf
        scale = math.exp(value["ln_scale"])
        prediction = scale * model + value["offset"]
        variance = self.error**2 + (value["jitter"] * scale * model) ** 2
        return prediction, variance, model, source, tau, physical

    def evaluate(self, theta):
        prior = self.log_prior(theta)
        if not math.isfinite(prior):
            return -math.inf, -math.inf
        try:
            prediction, variance, *_ = self.predict(theta)
        except (ValueError, FloatingPointError, OverflowError):
            return -math.inf, -math.inf
        if np.any(~np.isfinite(prediction)) or np.any(variance <= 0):
            return -math.inf, -math.inf
        residual = self.data - prediction
        likelihood = -0.5 * np.sum(
            residual**2 / variance + np.log(2 * math.pi * variance)
        )
        return float(prior + likelihood), float(likelihood)

    def log_probability(self, theta):
        return self.evaluate(theta)


def draw_walkers(posterior, count, rng):
    p, rows = posterior.config["priors"], []
    families = ((0.4, (-4, -1.5)), (0.2, (-1.5, -0.5)), (0.4, (-0.5, 1)))
    family = np.concatenate(
        [np.tile(bounds, (max(1, round(frac * count)), 1)) for frac, bounds in families]
    )[:count]
    if len(family) < count:
        family = np.vstack((family, np.tile((-0.5, 1), (count - len(family), 1))))
    rng.shuffle(family)

    def positive(mean, sigma):
        return next(x for x in rng.normal(mean, sigma, 1000) if x > 0)

    for lower, upper in family:
        row = [
            rng.normal(math.log(p[f"Qref_{x}"][0]), p[f"Qref_{x}"][1])
            for x in ("H2O", "CO2", "CO")
        ]
        row += [rng.normal(*p[f"beta_{x}"]) for x in ("H2O", "CO2", "CO")]
        row += [
            rng.uniform(*p["beta_H2"]),
            rng.uniform(lower, upper),
            positive(posterior.heavy_location, posterior.heavy_scale),
        ]
        row += [positive(*p["parent_velocity"]) for _ in range(4)]
        row += [
            rng.uniform(*p["OH_velocity"]),
            rng.uniform(*p["O_velocity"]),
            rng.uniform(*p["O_velocity"]),
            rng.uniform(*p["H_velocity"]),
            rng.normal(0, posterior.scale_sigma),
            rng.normal(0, posterior.offset_sigma),
            abs(rng.normal(0, p["jitter_sigma"])),
        ]
        rows.append(row)
    return np.asarray(rows)


def run_sampler(posterior, start, burn, retained, processes=1, seed=None):
    import emcee

    method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
    pool = mp.get_context(method).Pool(processes) if processes > 1 else None
    try:
        sampler = emcee.EnsembleSampler(
            len(start),
            len(COORDS),
            posterior.log_probability,
            pool=pool,
            blobs_dtype=float,
        )
        if seed is not None:
            sampler.random_state = np.random.RandomState(seed).get_state()
        state = sampler.run_mcmc(start, burn, progress=False)
        sampler.reset()
        sampler.run_mcmc(state, retained, progress=False)
    finally:
        if pool:
            pool.close()
            pool.join()
    return (
        sampler.get_chain(),
        sampler.get_log_prob(),
        sampler.get_blobs(),
        float(np.mean(sampler.acceptance_fraction)),
    )
