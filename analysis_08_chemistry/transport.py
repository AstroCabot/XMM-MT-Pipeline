import math

import numpy as np

DAY_S = 86400.0
GAUSS_K = 0.01720209894846
EV_ERG = 1.602176634e-12
ARCMIN2_SR = (math.pi / 10800) ** 2


class Orbit:
    def __init__(
        self, q_au, eccentricity, perihelion_jd, observation_jd, maximum_age_s
    ):
        self.q, self.e = q_au, eccentricity
        self.tp, self.now = perihelion_jd, observation_jd
        self.a = self.q / (1 - self.e)
        self.n = GAUSS_K * (-self.a) ** -1.5
        self.perihelion_age_s = (self.now - self.tp) * DAY_S
        nodes = np.geomspace(1e-3, maximum_age_s, 8190)
        self.age = np.unique(np.r_[0.0, nodes, self.perihelion_age_s])
        self.radius_grid = self.radius(self.age)
        inverse_r2 = self.radius_grid**-2
        self.solar_exposure = np.r_[
            0.0, np.cumsum(0.5 * (inverse_r2[1:] + inverse_r2[:-1]) * np.diff(self.age))
        ]
        self.current_r_au = float(self.radius(0.0))

    def radius(self, age_s):
        age = np.asarray(age_s, float)
        mean = self.n * (self.now - self.tp - age / DAY_S)
        anomaly = np.arcsinh(mean / self.e)
        for _ in range(12):
            step = (self.e * np.sinh(anomaly) - anomaly - mean) / (
                self.e * np.cosh(anomaly) - 1
            )
            anomaly -= step
            if np.max(np.abs(step)) < 1e-13:
                break
        else:
            raise RuntimeError("hyperbolic orbit solve failed")
        radius = self.a * (1 - self.e * np.cosh(anomaly))
        return float(radius) if age.ndim == 0 else radius

    def exposure(self, age_s):
        age = np.asarray(age_s, float)
        if np.any(age < 0):
            raise ValueError("reaction time precedes observation")
        return np.interp(np.minimum(age, self.age[-1]), self.age, self.solar_exposure)

    def inverse_r2(self, age_s):
        age = np.asarray(age_s, float)
        value = (
            np.interp(np.minimum(age, self.age[-1]), self.age, self.radius_grid) ** -2
        )
        return np.where(age <= self.age[-1], value, 0)

    def activity(self, beta, age_s):
        age = np.asarray(age_s, float)
        radius = np.interp(
            np.minimum(age, self.perihelion_age_s), self.age, self.radius_grid
        )
        return (radius / self.current_r_au) ** (-beta)


class Transport:
    def __init__(
        self,
        orbit,
        network,
        shell_count,
        inner_km,
        outer_km,
        source_radius_arcmin,
        km_per_arcmin,
        granddaughter_stride=4,
        granddaughter_radial_stride=4,
    ):
        self.orbit, self.network = orbit, network
        self.edges = np.geomspace(inner_km, outer_km, shell_count + 1)
        self.radius = np.sqrt(self.edges[:-1] * self.edges[1:])
        self.source_radius = np.asarray(source_radius_arcmin, float)
        rho = self.source_radius[None, :] * km_per_arcmin
        inner, outer = self.edges[:-1, None], self.edges[1:, None]
        self.chord_cm = 2e5 * (
            np.sqrt(np.maximum(outer**2 - rho**2, 0))
            - np.sqrt(np.maximum(inner**2 - rho**2, 0))
        )
        self.stride = granddaughter_stride
        self.radial_stride = granddaughter_radial_stride

    def density(self, chain, beta, speed):
        segments, weight = chain
        family = [item[0] for item in segments]
        velocity = [speed[item[1]] for item in segments]
        rate = [self.network[name]["rate_1au_s-1"] for name in family]
        r0, radius = self.edges[0], self.radius
        norm = weight / (4 * math.pi * (radius * 1e5) ** 2 * velocity[-1] * 1e5)
        if len(segments) == 1:
            age = (radius - r0) / velocity[0]
            survival = np.exp(-rate[0] * self.orbit.exposure(age))
            return norm * self.orbit.activity(beta, age) * survival
        if len(segments) == 2:
            lo, hi = self.edges[:-1], self.edges[1:]
            upper = np.minimum(hi[:, None], radius[None, :])
            valid = upper > lo[:, None]
            birth = 0.5 * (lo[:, None] + upper)
            final_age = np.maximum(radius[None, :] - birth, 0) / velocity[1]
            emission_age = final_age + (birth - r0) / velocity[0]
            probability = (
                (upper - lo[:, None])
                * rate[0]
                / velocity[0]
                * self.orbit.inverse_r2(final_age)
                * np.exp(
                    -rate[0]
                    * (
                        self.orbit.exposure(emission_age)
                        - self.orbit.exposure(final_age)
                    )
                )
            )
            survival = np.exp(-rate[1] * self.orbit.exposure(final_age))
            terms = np.where(
                valid,
                probability * survival * self.orbit.activity(beta, emission_age),
                0,
            )
            return norm * terms.sum(axis=0)
        return norm * self._granddaughter(rate, velocity, beta)

    def _granddaughter(self, rate, velocity, beta):
        coarse = self.edges[:: self.stride]
        if coarse[-1] != self.edges[-1]:
            coarse = np.r_[coarse, self.edges[-1]]
        lo1, hi1 = coarse[:-1][None, :, None], coarse[1:][None, :, None]
        lo2, hi2 = coarse[:-1][None, None, :], coarse[1:][None, None, :]
        index = np.unique(
            np.r_[
                np.arange(0, len(self.radius), self.radial_stride), len(self.radius) - 1
            ]
        )
        final_radius = self.radius[index, None, None]
        b1 = 0.5 * (lo1 + hi1)
        upper2 = np.minimum(hi2, final_radius)
        lower2 = np.maximum(lo2, b1)
        valid = upper2 > lower2
        lower2, upper2 = np.where(valid, lower2, b1), np.where(valid, upper2, b1)
        b2 = 0.5 * (lower2 + upper2)
        final_age = np.maximum(final_radius - b2, 0) / velocity[2]
        middle_age = (b2 - b1) / velocity[1]
        parent_start = final_age + middle_age
        emission_age = parent_start + (b1 - self.edges[0]) / velocity[0]
        parent = (
            (hi1 - lo1)
            * rate[0]
            / velocity[0]
            * self.orbit.inverse_r2(parent_start)
            * np.exp(
                -rate[0]
                * (
                    self.orbit.exposure(emission_age)
                    - self.orbit.exposure(parent_start)
                )
            )
        )
        daughter = (
            (upper2 - lower2)
            * rate[1]
            / velocity[1]
            * self.orbit.inverse_r2(final_age)
            * np.exp(
                -rate[1]
                * (self.orbit.exposure(parent_start) - self.orbit.exposure(final_age))
            )
        )
        survival = np.exp(-rate[2] * self.orbit.exposure(final_age))
        terms = parent * daughter * survival * self.orbit.activity(beta, emission_age)
        values = np.where(valid, terms, 0).sum(axis=(1, 2))
        return np.exp(
            np.interp(
                np.log(self.radius),
                np.log(self.radius[index]),
                np.log(np.maximum(values, np.finfo(float).tiny)),
            )
        )

    def source_profile(
        self, primaries, beta, speed, chains, phi, photon_energy_ev, photons_per_ion
    ):
        opacity = np.zeros_like(self.radius)
        for primary, q_s in primaries.items():
            for chain in chains[primary]:
                density = q_s * self.density(chain, beta[primary], speed)
                opacity += density * self.network[chain[0][-1][0]]["sigma_cm2"]
        tau = opacity @ self.chord_cm
        saturation = (
            phi
            * photons_per_ion
            * photon_energy_ev
            * EV_ERG
            * ARCMIN2_SR
            / (4 * math.pi)
        )
        return saturation * -np.expm1(-tau), tau


def tau_one_radius(radius, tau):
    thick = np.flatnonzero(tau >= 1)
    if not len(thick):
        return 0.0, "below_grid"
    index = int(thick[-1])
    if index == len(tau) - 1:
        return float(radius[-1]), "above_grid"
    value = np.interp(
        0.0,
        np.log(tau[index : index + 2])[::-1],
        np.log(radius[index : index + 2])[::-1],
    )
    return float(np.exp(value)), "resolved"
