import math
import multiprocessing as mp
import pickle
from pathlib import Path

import numpy as np
from scipy.special import ndtr, ndtri

CHECKPOINT_STEPS = 500

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
        heavy_lower,
        heavy_upper,
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
        self.heavy_lower, self.heavy_upper = heavy_lower, heavy_upper
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
        result += uniform_logpdf(v["phi_heavy"], self.heavy_lower, self.heavy_upper)
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
            rng.uniform(posterior.heavy_lower, posterior.heavy_upper),
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


def load_checkpoint(path):
    if path is None or not Path(path).exists():
        return None
    with Path(path).open("rb") as stream:
        return pickle.load(stream)


def save_checkpoint(backend, path):
    temporary = Path(path).with_suffix(".tmp")
    with temporary.open("wb") as stream:
        pickle.dump(backend, stream)
    temporary.replace(path)


def run_sampler(
    posterior, start, burn, retained, processes=1, seed=None, checkpoint=None
):
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
            moves=emcee.moves.DEMove(),
            backend=load_checkpoint(checkpoint),
        )
        if seed is not None and not sampler.iteration:
            sampler.random_state = np.random.RandomState(seed).get_state()
        while sampler.iteration < burn + retained:
            steps = min(CHECKPOINT_STEPS, burn + retained - sampler.iteration)
            state = None if sampler.iteration else start
            sampler.run_mcmc(state, steps, progress=False)
            if checkpoint is not None:
                save_checkpoint(sampler.backend, checkpoint)
    finally:
        if pool:
            pool.close()
            pool.join()
    chain = sampler.get_chain(discard=burn)
    return (
        chain,
        sampler.get_log_prob(discard=burn),
        sampler.get_blobs(discard=burn),
        float(np.mean(np.any(chain[1:] != chain[:-1], axis=-1))),
    )


class Gaussianized:
    def __init__(self, posterior):
        self.posterior = posterior
        p = posterior.config["priors"]
        self.normal = {
            COORDS.index("ln_scale"): (0.0, posterior.scale_sigma),
            COORDS.index("offset"): (0.0, posterior.offset_sigma),
        }
        for name in ("H2O", "CO2", "CO"):
            median, width = p[f"Qref_{name}"]
            self.normal[COORDS.index(f"ln_Qref_{name}")] = (math.log(median), width)
            self.normal[COORDS.index(f"beta_{name}")] = tuple(p[f"beta_{name}"])
        self.uniform = {
            COORDS.index("beta_H2"): tuple(p["beta_H2"]),
            COORDS.index("log10_H2_ratio"): tuple(p["log10_H2_ratio"]),
            COORDS.index("phi_heavy"): (posterior.heavy_lower, posterior.heavy_upper),
            COORDS.index("v_OH"): tuple(p["OH_velocity"]),
            COORDS.index("v_O_H2O"): tuple(p["O_velocity"]),
            COORDS.index("v_O_OH"): tuple(p["O_velocity"]),
            COORDS.index("v_H"): tuple(p["H_velocity"]),
        }
        self.speeds = [COORDS.index(f"v_{name}") for name in ("H2O", "CO2", "CO", "H2")]
        self.speed_mean, self.speed_sigma = p["parent_velocity"]
        self.speed_floor = float(ndtr(-self.speed_mean / self.speed_sigma))
        self.jitter = COORDS.index("jitter")
        self.jitter_sigma = p["jitter_sigma"]

    def physical(self, latent):
        latent = np.asarray(latent, float)
        theta = np.empty_like(latent)
        for index, (mean, sigma) in self.normal.items():
            theta[..., index] = mean + sigma * latent[..., index]
        for index, (lower, upper) in self.uniform.items():
            theta[..., index] = lower + (upper - lower) * ndtr(latent[..., index])
        for index in self.speeds:
            level = self.speed_floor + (1 - self.speed_floor) * ndtr(latent[..., index])
            theta[..., index] = self.speed_mean + self.speed_sigma * ndtri(level)
        level = 0.5 * (1 + ndtr(latent[..., self.jitter]))
        theta[..., self.jitter] = self.jitter_sigma * ndtri(level)
        return theta

    def latent(self, theta):
        theta = np.asarray(theta, float)
        latent = np.empty_like(theta)
        for index, (mean, sigma) in self.normal.items():
            latent[..., index] = (theta[..., index] - mean) / sigma
        for index, (lower, upper) in self.uniform.items():
            latent[..., index] = ndtri((theta[..., index] - lower) / (upper - lower))
        for index in self.speeds:
            level = ndtr((theta[..., index] - self.speed_mean) / self.speed_sigma)
            latent[..., index] = ndtri(
                (level - self.speed_floor) / (1 - self.speed_floor)
            )
        level = 2 * ndtr(theta[..., self.jitter] / self.jitter_sigma) - 1
        latent[..., self.jitter] = ndtri(level)
        return latent

    def log_probability(self, latent):
        _, likelihood = self.posterior.evaluate(self.physical(latent))
        if not math.isfinite(likelihood):
            return -math.inf, -math.inf
        return likelihood - 0.5 * float(np.sum(np.square(latent))), likelihood


def sample(posterior, start, burn, retained, processes=1, seed=None, checkpoint=None):
    model = Gaussianized(posterior)
    latent, _, likelihood, acceptance = run_sampler(
        model, model.latent(start), burn, retained, processes, seed, checkpoint
    )
    chain = model.physical(latent)
    prior = np.array([[posterior.log_prior(row) for row in step] for step in chain])
    return chain, prior + likelihood, likelihood, acceptance
