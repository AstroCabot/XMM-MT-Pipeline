import math

from scipy.integrate import quad
from scipy.optimize import brentq

EV_TO_ERG = 1.602176634e-12
AU_KM = 1.495978707e8
AVOGADRO = 6.02214076e23
BOLTZMANN = 1.380649e-23
ARCMIN2_PER_SR = (180 * 60 / math.pi) ** 2

SIGMA_CM2 = 3e-15
LAMBDA_CM = 5e10
KAPPA = 8.1
PROTON_DENSITY_1AU_CM3 = 7.0
SOLAR_WIND_SPEED_CM_S = 4e7
HEAVY_ION_FRACTION = 1e-3
PHOTONS_PER_ION = 1.0
PHOTON_ENERGY_EV = 550.0
WEGMANN_C_S2_CM2 = 1e-38

REFERENCE_Q_S = 1.5e29
REFERENCE_RADIUS_KM = 2.4
ICE_PROPERTIES = {
    "H2O": (54.46, 8 / 6, 155),
    "CO2": (28.84, 8 / 6, 82),
    "CO": (7.60, 7 / 5, 58),
    "N2": (7.34, 7 / 5, 25),
    "H2": (1.00, 7 / 5, 6),
}


def molecular_energy_kj(delta_h_kj_mol, gamma, temperature_k):
    return delta_h_kj_mol / AVOGADRO + gamma * BOLTZMANN * temperature_k / 1000


WATER_ENERGY_KJ = molecular_energy_kj(*ICE_PROPERTIES["H2O"])


def proton_density_cm3(r_au):
    return PROTON_DENSITY_1AU_CM3 * r_au**-2


def ion_number_flux_cm2_s(r_au):
    return proton_density_cm3(r_au) * SOLAR_WIND_SPEED_CM_S * HEAVY_ION_FRACTION


def ion_energy_flux(phi_cm2_s, energy_ev=PHOTON_ENERGY_EV):
    return phi_cm2_s * PHOTONS_PER_ION * energy_ev * EV_TO_ERG


def cabot_energy_flux(r_au, energy_ev=PHOTON_ENERGY_EV):
    return ion_energy_flux(ion_number_flux_cm2_s(r_au), energy_ev)


def production_rate(radius_km, r_au, energy_kj):
    return (
        REFERENCE_Q_S
        * (radius_km / REFERENCE_RADIUS_KM) ** 2
        * r_au**-2
        * WATER_ENERGY_KJ
        / energy_kj
    )


def thick_radius_cm(q_s, speed_km_s=1.0, sigma_cm2=SIGMA_CM2, length_cm=LAMBDA_CM):
    speed = speed_km_s * 1e5
    linear = sigma_cm2 * q_s / (4 * speed)
    coefficient = sigma_cm2 * q_s / (4 * math.pi * speed)

    def residual(radius):
        integral = quad(
            lambda angle: math.exp(-radius / (length_cm * math.cos(angle))),
            0,
            math.pi / 2,
            epsabs=1e-11,
            epsrel=1e-11,
        )[0]
        return 2 * coefficient * integral / radius - 1

    return brentq(residual, linear * 1e-12, linear * (1 + 1e-10), rtol=1e-12)


def luminosity(q_s, r_au, speed_km_s=1.0, energy_flux=None):
    radius = KAPPA * thick_radius_cm(q_s, speed_km_s)
    energy_flux = cabot_energy_flux(r_au) if energy_flux is None else energy_flux
    return math.pi * radius**2 * energy_flux


def invert_luminosity(luminosity_erg_s, r_au, speed_km_s=1.0, energy_flux=None):
    def residual(log_q):
        return math.log(
            luminosity(10**log_q, r_au, speed_km_s, energy_flux)
        ) - math.log(luminosity_erg_s)

    return 10 ** brentq(residual, 24, 32, rtol=1e-12)


def empirical_rate(luminosity_erg_s, energy_flux):
    return math.sqrt(luminosity_erg_s / (WEGMANN_C_S2_CM2 * energy_flux))


def thin_rate(surface_flux_arcmin2, radius_cm, speed_km_s, energy_flux, sigma_cm2):
    intensity_sr = surface_flux_arcmin2 * ARCMIN2_PER_SR
    return (
        16
        * math.pi
        * speed_km_s
        * 1e5
        * radius_cm
        * intensity_sr
        / (energy_flux * sigma_cm2)
    )


def thin_tau(q_s, radius_cm, speed_km_s, sigma_cm2):
    return sigma_cm2 * q_s / (4 * speed_km_s * 1e5 * radius_cm)
