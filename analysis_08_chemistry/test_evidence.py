import numpy as np
from evidence import estimate

rng = np.random.default_rng(20261130)
for ndim in (17, 20):
    center = rng.normal(scale=0.5, size=ndim)
    width = np.geomspace(0.05, 0.5, ndim)
    precision = 1 / width**2 + 1
    mean = center / width**2 / precision
    samples = mean + rng.normal(size=(48000, ndim)) / np.sqrt(precision)

    def log_likelihood(latent):
        return 7.3 - 0.5 * np.sum(
            ((latent - center) / width) ** 2 + np.log(2 * np.pi * width**2), axis=-1
        )

    exact = 7.3 - 0.5 * np.sum(
        center**2 / (width**2 + 1) + np.log(2 * np.pi * (width**2 + 1))
    )
    joint = log_likelihood(samples) - 0.5 * (
        np.sum(samples**2, axis=-1) + ndim * np.log(2 * np.pi)
    )
    for mass in (0.5, 0.3):
        value = estimate(log_likelihood, samples, joint, mass, 4000, 1)["log_evidence"]
        assert abs(value - exact) < 0.05, (ndim, mass, value, exact)
