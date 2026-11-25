import math
from functools import lru_cache

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import exprel

DAY_S = 86400.0
GAUSS_K = 0.01720209894846
EV_ERG = 1.602176634e-12
ARCMIN2_SR = (math.pi / 10800) ** 2
TINY = np.finfo(float).tiny
REMNANT_DEPTH = 30.0


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
        self.current_r_au = float(self.radius(0.0))
        self.log_radius_grid = np.log(self.radius_grid / self.current_r_au)
        inverse_r2 = self.radius_grid**-2
        self.solar_exposure = np.r_[
            0.0, np.cumsum(0.5 * (inverse_r2[1:] + inverse_r2[:-1]) * np.diff(self.age))
        ]

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

    def activity(self, beta, age_s):
        age = np.asarray(age_s, float)
        return np.exp(
            -beta
            * np.interp(
                np.minimum(age, self.perihelion_age_s), self.age, self.log_radius_grid
            )
        )


def intervals(value, log_value, exposure, step):
    change = log_value[:, step:] - log_value[:, :-step]
    flat = np.abs(change) < 1e-6
    mean = np.where(
        flat,
        0.5 * (value[:, step:] + value[:, :-step]),
        (value[:, step:] - value[:, :-step]) / np.where(flat, 1, change),
    )
    return (exposure[:, :-step] - exposure[:, step:]) * mean


def sequential(old, new, span, exposure, growth, decay):
    x, y = growth * span, decay * span
    x = np.where(np.abs(x) < 1e-6, 1e-6, x)
    shape = (exprel(x + y) - exprel(y)) / np.expm1(x)
    return np.exp(decay * exposure) * span * (new * exprel(y) + (old - new) * shape)


def cumulate(single, paired):
    total = np.zeros((single.shape[0], single.shape[1] + 1))
    np.cumsum(single, axis=1, out=total[:, 1:])
    coarse = total.copy()
    coarse[:, 2::2] = np.cumsum(paired[:, ::2], axis=1)
    coarse[:, 3::2] = total[:, 1:2] + np.cumsum(paired[:, 1::2], axis=1)
    return total + np.clip((total - coarse) / 3, -0.1 * total, 0.1 * total)


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
        node_stride=2,
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
        self.node = self.edges[::node_stride]
        self.step = math.log(self.node[1] / self.node[0])
        self.separation = np.maximum(self.node[:, None] - self.node[None, :], 0)

    @lru_cache(maxsize=16)
    def _flight(self, speed):
        return self.orbit.exposure(self.separation / speed)

    @lru_cache(maxsize=16)
    def _release(self, parent_speed, speed, beta):
        age = self.separation / speed + (self.node - self.node[0]) / parent_speed
        return np.log(self.orbit.activity(beta, age)), self.orbit.exposure(age)

    @lru_cache(maxsize=4)
    def _daughter(self, parent_rate, rate, parent_speed, speed, beta):
        flight = self._flight(speed)
        log_activity, release = self._release(parent_speed, speed, beta)
        log_birth = log_activity - parent_rate * (release - flight) - rate * flight
        birth = np.exp(log_birth)
        flux = (parent_rate * speed / parent_speed) * cumulate(
            intervals(birth, log_birth, flight, 1),
            intervals(birth, log_birth, flight, 2),
        )
        log_flux = np.log(np.maximum(flux, TINY)) + rate * flight
        with np.errstate(divide="ignore", invalid="ignore"):
            width = np.diff(flight, axis=0)
            secant = np.diff(log_flux, axis=0) / width
            lower, upper = secant[:-1], secant[1:]
            first, second = 2 * width[1:] + width[:-1], width[1:] + 2 * width[:-1]
            slope = np.empty_like(flux)
            slope[1:-1] = np.where(
                lower * upper > 0,
                (first + second) / (first / lower + second / upper),
                0,
            )
        slope[0], slope[-1] = secant[0], secant[-1]
        diagonal = np.arange(len(self.node) - 1)
        slope[diagonal, diagonal] = secant[diagonal, diagonal]
        return flux, log_flux, np.nan_to_num(slope)

    def _granddaughter(self, daughter, rate, velocity):
        flux, log_flux, slope = daughter
        reach, flight = self._flight(velocity[1]), self._flight(velocity[2])
        offset, last = self.node - self.node[0], len(self.node) - 1
        column = np.broadcast_to(np.arange(last + 1), flight.shape)
        arrival = self.node + velocity[1] / velocity[2] * self.separation
        start = np.floor(np.log(arrival / self.node[0]) / self.step).astype(int)
        start = np.minimum(np.maximum(start, column), last - 1)
        width = reach[start + 1, column] - reach[start, column]
        with np.errstate(divide="ignore", invalid="ignore"):
            u = np.where(width > 0, (flight - reach[start, column]) / width, 1.0)
        u = np.clip(u, 0, 1)
        log_value = (
            (1 + 2 * u) * (1 - u) ** 2 * log_flux[start, column]
            + u * (1 - u) ** 2 * width * slope[start, column]
            + u**2 * (3 - 2 * u) * log_flux[start + 1, column]
            + u**2 * (u - 1) * width * slope[start + 1, column]
            - rate[1] * flight
        )
        growth, decay = rate[0] - rate[1], rate[1] - rate[2]
        burst = rate[1] * width * offset > REMNANT_DEPTH * (
            self.node[start + 1] - self.node[start]
        )
        lower, upper = (
            flux[start[burst], column[burst]],
            flux[start[burst] + 1, column[burst]],
        )
        x = growth * width[burst]
        x = np.where(np.abs(x) < 1e-6, 1e-6, x)
        frozen = lower + (upper - lower) * np.expm1(x * u[burst]) / np.expm1(x)
        log_value[burst] = np.log(np.maximum(frozen, TINY))
        log_birth = log_value + decay * flight
        birth = np.exp(log_birth)
        parts = []
        for step in (1, 2):
            part = intervals(birth, log_birth, flight, step)
            span = flight[:, :-step] - flight[:, step:]
            burst = rate[1] * velocity[2] * span * offset[:-step] > (
                REMNANT_DEPTH * velocity[1] * (self.node[step:] - self.node[:-step])
            )
            part[burst] = sequential(
                np.exp(log_value[:, :-step][burst]),
                np.exp(log_value[:, step:][burst]),
                span[burst],
                flight[:, step:][burst],
                growth,
                decay,
            )
            parts.append(part)
        return rate[1] * velocity[2] / velocity[1] * np.diagonal(cumulate(*parts))

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
        daughter = self._daughter(rate[0], rate[1], velocity[0], velocity[1], beta)
        if len(segments) == 2:
            flux = np.diagonal(daughter[0])
        else:
            flux = self._granddaughter(daughter, rate, velocity)
        return norm * np.exp(
            PchipInterpolator(
                np.log(self.node[1:] - r0), np.log(np.maximum(flux[1:], TINY))
            )(np.log(radius - r0))
        )

    def source_profile(
        self, primaries, beta, speed, chains, phi, photon_energy_ev, photons_per_ion
    ):
        opacity = np.zeros_like(self.radius)
        for primary, q_s in primaries.items():
            if q_s == 0:
                continue
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
