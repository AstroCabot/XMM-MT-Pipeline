import hashlib
import json

import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import logsumexp
from scipy.stats import chi2

from config import CONFIG, SHARED
from inference import Gaussianized, coordinates

FITS = {
    "parents_shallow": ("no_direct_h2", False, "beta"),
    "radiolytic_shallow": ("", True, "beta"),
    "parents_steep": ("no_direct_h2/water_frincke", False, "beta_steep"),
    "radiolytic_steep": ("water_frincke", True, "beta_steep"),
}
MASSES = (0.5, 0.3)
DRAWS = 4000
SEED = 20261127


def radius2(points, mean, factor):
    return np.sum(solve_triangular(factor, (points - mean).T, lower=True) ** 2, 0)


def estimate(log_likelihood, samples, log_joint, mass, draws, seed):
    ndim = samples.shape[-1]
    core = samples[log_joint >= np.median(log_joint)]
    mean = core.mean(axis=0)
    factor = np.linalg.cholesky(np.cov(core, rowvar=False))
    cut = chi2.ppf(mass, ndim)
    fraction = np.mean(radius2(samples, mean, factor) <= cut)
    rng = np.random.default_rng(seed)
    latent = np.empty((0, ndim))
    while len(latent) < draws:
        proposal = mean + rng.normal(size=(4 * draws, ndim)) @ factor.T
        latent = np.vstack([latent, proposal[radius2(proposal, mean, factor) <= cut]])
    latent = latent[:draws]
    log_reference = -0.5 * radius2(latent, mean, factor) - np.log(mass)
    log_reference -= np.log(np.diag(factor)).sum() + 0.5 * ndim * np.log(2 * np.pi)
    log_prior = -0.5 * (np.sum(latent**2, axis=1) + ndim * np.log(2 * np.pi))
    weights = log_likelihood(latent) + log_prior - log_reference
    total = logsumexp(weights)
    return {
        "gaussian_mass": mass,
        "seed": seed,
        "draws": draws,
        "posterior_fraction_inside": float(fraction),
        "effective_draws": float(np.exp(2 * total - logsumexp(2 * weights))),
        "log_evidence": float(total - np.log(draws) - np.log(fraction)),
    }


def main():
    import run

    inputs = run.load_inputs()
    operator = run.load_psf(run.OUTPUT / "psf.npz", inputs, run.psf_dependencies())
    results = {}
    for name, (folder, direct, key) in FITS.items():
        directory = run.OUTPUT / folder
        beta = SHARED["spherex"]["parent_priors"]["H2O"][key]
        CONFIG["direct_h2"] = direct
        CONFIG["priors"]["beta_H2O"] = [beta["mean"], beta["sigma"]]
        result = json.loads((directory / "result.json").read_text())
        if (
            result.get("direct_h2", True) != direct
            or list(result["beta_H2O_prior"]) != CONFIG["priors"]["beta_H2O"]
        ):
            raise ValueError(f"{directory}: chain priors differ from {name}")
        posterior, _ = run.build_model(inputs, operator)
        transform = Gaussianized(posterior)
        with np.load(directory / "chain.npz") as chain:
            if tuple(chain["coordinate_names"]) != coordinates(CONFIG):
                raise ValueError(f"{directory}: chain coordinates differ from priors")
            samples = transform.latent(chain["samples"].astype(float))
            log_joint = chain["log_likelihood"].astype(float)
        log_joint -= 0.5 * (
            np.sum(samples**2, axis=-1) + samples.shape[-1] * np.log(2 * np.pi)
        )

        def log_likelihood(latent):
            return np.array(
                [posterior.evaluate(row)[1] for row in transform.physical(latent)]
            )

        reference, *checks = [
            estimate(log_likelihood, samples, log_joint, mass, DRAWS, SEED + index)
            for index, mass in enumerate(MASSES)
        ]
        results[name] = {
            "chain_sha256": hashlib.sha256(
                (directory / "chain.npz").read_bytes()
            ).hexdigest(),
            "samples": len(samples),
            "dimensions": samples.shape[-1],
            **reference,
            "checks": checks,
        }
        print(f"{name}: ln Z = {reference['log_evidence']:.3f}", flush=True)
    for history in ("shallow", "steep"):
        parents = results[f"parents_{history}"]
        h2 = results[f"radiolytic_{history}"]
        h2["delta_log_evidence"] = h2["log_evidence"] - parents["log_evidence"]
    for name, (folder, *_) in FITS.items():
        (run.OUTPUT / folder / "evidence.json").write_text(
            json.dumps(results[name], indent=2, allow_nan=False) + "\n"
        )


if __name__ == "__main__":
    main()
