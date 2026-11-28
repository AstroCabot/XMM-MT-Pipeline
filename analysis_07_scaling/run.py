import csv
import json
import math
from pathlib import Path

import numpy as np

from analysis_05_radial_profile.fit import integrated_broken

from .model import (
    AU_KM,
    EV_TO_ERG,
    HEAVY_ION_FRACTION,
    ICE_PROPERTIES,
    KAPPA,
    LAMBDA_CM,
    PHOTON_ENERGY_EV,
    PHOTONS_PER_ION,
    PROTON_DENSITY_1AU_CM3,
    SIGMA_CM2,
    SOLAR_WIND_SPEED_CM_S,
    WEGMANN_C_S2_CM2,
    WEGMANN_FAST_H_ERG_CM2_S,
    cabot_energy_flux,
    coefficient_scatter_dex,
    empirical_rate,
    ion_energy_flux,
    ion_number_flux_cm2_s,
    invert_luminosity,
    luminosity,
    molecular_energy_kj,
    production_rate,
    rate_range,
    thin_rate,
    thin_tau,
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "outputs"
FINAL_PRODUCTS = ("result.json", "scaling.tsv")
CALIBRATION = HERE / "wegmann2004_fig19.tsv"
TASK5 = ROOT / "analysis_05_radial_profile/outputs/result.json"
TASK3 = ROOT / "analysis_03_spectrum/outputs/result.json"
TASK6 = ROOT / "analysis_06_solar_wind/outputs/result.json"
MORPHOLOGY = ROOT / "analysis_04_morphology/outputs"
SHARED = json.loads((ROOT / "parameters.json").read_text())

SPHEREX = SHARED["spherex"]
REFERENCE_R_AU = float(SPHEREX["reference_distance_au"])
REFERENCE_EPOCH_UTC = SPHEREX["reference_epoch_utc"]
END_MEMBER_RADIUS_KM = 1.3
END_MEMBER_SPEED_KM_S = 0.37
H2_SIGMA_CM2 = float(SHARED["h2_single_capture_cross_section_cm2"])
SCIENCE_APERTURE = float(SHARED["science_aperture_radius_arcsec"])
SCIENCE_NAME = f"{SCIENCE_APERTURE:g}_arcsec"
WEGMANN_APERTURE_KM = 1.1e5
WIND_MONITORS = ("stereo_a", "solar_orbiter", "wind")
ION_INSTRUMENTS = ("ace_swics", "solar_orbiter_his")
THIN_RADIUS_ARCMIN = 6.0
PARENTS = {
    name: (values["reference_rate_s-1"]["median"], values["beta"]["mean"])
    for name, values in SPHEREX["parent_priors"].items()
}


def physical_radius_km(arcsec, distance_au):
    return distance_au * AU_KM * math.tan(math.radians(arcsec / 3600))


def parent_rates(r_au):
    xmm = {
        name: rate * (r_au / REFERENCE_R_AU) ** -beta
        for name, (rate, beta) in PARENTS.items()
    }
    return xmm, sum(rate for rate, _ in PARENTS.values()), sum(xmm.values())


def missing(inferred, measured):
    value = np.maximum(0.0, inferred - measured)
    return value, value / inferred


def read_calibration(path=CALIBRATION):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    hg2 = np.array([float(row["h_g2_1e54_erg_s-3_cm-2"]) for row in rows])
    lx = np.array([float(row["luminosity_1e16_erg_s"]) for row in rows])
    if len(rows) != 48 or np.any(hg2 <= 0) or np.any(lx <= 0):
        raise ValueError("unexpected Wegmann Fig. 19 points")
    return lx * 1e16 / (WEGMANN_C_S2_CM2 * hg2 * 1e54)


def aperture_luminosity(task5, task3, radius_km):
    broken = task5["epic"]["broken"]
    distance = float(task3["distance_au"])
    arcsec = math.degrees(math.atan(radius_km / (distance * AU_KM))) * 3600
    rate = integrated_broken(
        float(broken["normalization_at_break_ct_s_arcmin2"]),
        float(broken["alpha_inner"]),
        float(broken["alpha_outer"]),
        float(broken["break_radius_arcmin"]),
        arcsec / 60,
    )
    flux = rate * float(task3["mos2_ecf_erg_cm2_per_count"])
    return arcsec, flux * 4 * math.pi * (distance * AU_KM * 1e5) ** 2


def measured_energy_flux(task6, photon_energy_ev):
    fluxes = np.log10([task6["proton_flux"][name]["median"] for name in WIND_MONITORS])
    ions = np.log10(
        [task6["ion_fraction"][name]["f_16_50_84"] for name in ION_INSTRUMENTS]
    )
    flux, fraction = 10 ** fluxes.mean(), 10 ** ions[:, 1].mean()
    flux_sigma = float(fluxes.std(ddof=1))
    fraction_sigma = float(np.mean(ions[:, 2] - ions[:, 0]) / 2)
    inputs = {
        "proton_flux_cm-2_s-1": flux,
        "proton_flux_log10_sigma": flux_sigma,
        "ion_fraction": fraction,
        "ion_fraction_log10_sigma": fraction_sigma,
        "photon_energy_eV": photon_energy_ev,
    }
    energy_flux = fraction * flux * photon_energy_ev * EV_TO_ERG
    return inputs, energy_flux, math.hypot(flux_sigma, fraction_sigma)


def model_surface_flux(task5, task3, radius_arcmin):
    profile = task5["epic"]
    rate = float(profile["normalization_ct_s_arcmin2"]) * (
        radius_arcmin ** -float(profile["slope"])
    )
    return rate * float(task3["mos2_ecf_erg_cm2_per_count"])


def calculate(task5, task3, task6, ratios=None):
    fit_band = task3.get("fit_band_keV", [])
    if (
        task5.get("band_eV") != SHARED["soft_band_ev"]
        or len(fit_band) != 2
        or any(
            not math.isclose(value * 1000, expected, abs_tol=1e-9)
            for value, expected in zip(fit_band, SHARED["soft_band_ev"])
        )
    ):
        raise RuntimeError("upstream products differ from the shared soft band")
    r_au = float(task3["heliocentric_distance_au"])
    distance = float(task3["distance_au"])
    mean_energy = float(task3["mean_photon_energy_eV"])
    task6_prior = task6["heavy_ion_prior"]
    task6_bounds = [float(task6_prior["lower"]), float(task6_prior["upper"])]
    task6_phi = sum(task6_bounds) / 2
    task6_h = ion_energy_flux(task6_phi, mean_energy)
    cabot_phi = ion_number_flux_cm2_s(r_au)
    cabot_h = cabot_energy_flux(r_au)
    xmm_rates, reference_sum, xmm_sum = parent_rates(r_au)
    ratios = read_calibration() if ratios is None else np.asarray(ratios, float)
    wegmann_arcsec, wegmann_lx = aperture_luminosity(task5, task3, WEGMANN_APERTURE_KM)
    lx_sigma = float(task3["mos2_ecf_sigma_ln"]) / math.log(10)
    c_sigma = coefficient_scatter_dex(ratios)
    measured, measured_h, measured_sigma = measured_energy_flux(task6, mean_energy)
    energy_fluxes = {
        "w04_fast_wind": (WEGMANN_FAST_H_ERG_CM2_S, 0.0),
        "w04_fast_wind_scaled": (WEGMANN_FAST_H_ERG_CM2_S / r_au**2, 0.0),
        "measured": (measured_h, measured_sigma),
    }

    rows, estimates = [], {}
    for name, (energy_flux, h_sigma) in energy_fluxes.items():
        q = empirical_rate(wegmann_lx, energy_flux)
        q_range = rate_range(q, lx_sigma, c_sigma, h_sigma)
        deficit, fraction = missing(q, xmm_sum)
        deficit_range, fraction_range = missing(q_range, xmm_sum)
        estimates[name] = {
            "energy_flux_erg_cm-2_s-1": energy_flux,
            "energy_flux_log10_sigma": h_sigma,
            "inferred_production_rate_s-1": q,
            "inferred_production_rate_range_s-1": q_range.tolist(),
            "missing_rate_s-1": deficit,
            "missing_rate_range_s-1": deficit_range.tolist(),
            "missing_fraction": fraction,
            "missing_fraction_range": fraction_range.tolist(),
        }
        rows.append(
            {
                "kind": "wegmann",
                "name": name,
                "method": "Wegmann Eq. 30",
                "radius_arcsec": wegmann_arcsec,
                "radius_km": WEGMANN_APERTURE_KM,
                "measured_luminosity_erg_s": wegmann_lx,
                "inferred_production_rate_s-1": q,
                "xmm_parent_rate_s-1": xmm_sum,
                "missing_rate_s-1": deficit,
                "missing_fraction": fraction,
            }
        )

    observed = float(task3["regions"]["inner"]["luminosity_erg_s"])
    q = invert_luminosity(observed, r_au)
    q_task6 = invert_luminosity(observed, r_au, energy_flux=task6_h)
    deficit, fraction = missing(q, xmm_sum)
    radius_km = physical_radius_km(SCIENCE_APERTURE, distance)
    row = {
        "kind": "aperture",
        "name": SCIENCE_NAME,
        "method": "Cabot B2-B7/B11-B12; v=1 km s-1",
        "radius_arcsec": SCIENCE_APERTURE,
        "radius_km": radius_km,
        "measured_luminosity_erg_s": observed,
        "inferred_production_rate_s-1": q,
        "task06_prior_inferred_rate_s-1": q_task6,
        "xmm_parent_rate_s-1": xmm_sum,
        "missing_rate_s-1": deficit,
        "missing_fraction": fraction,
        "required_heavy_ion_flux_multiple_for_xmm_parents": observed
        / luminosity(xmm_sum, r_au),
    }
    rows.append(row)
    apertures = {SCIENCE_NAME: row}
    slow_q = invert_luminosity(observed, r_au, END_MEMBER_SPEED_KM_S)
    slow_missing, slow_fraction = missing(slow_q, xmm_sum)
    speed_check = {
        "speed_km_s": END_MEMBER_SPEED_KM_S,
        "inferred_production_rate_s-1": slow_q,
        "missing_rate_s-1": slow_missing,
        "missing_fraction": slow_fraction,
        "required_heavy_ion_flux_multiple_for_xmm_parents": observed
        / luminosity(xmm_sum, r_au, END_MEMBER_SPEED_KM_S),
    }

    end_members = {}
    for name, (enthalpy, gamma, temperature) in ICE_PROPERTIES.items():
        energy = molecular_energy_kj(enthalpy, gamma, temperature)
        predicted_q = production_rate(END_MEMBER_RADIUS_KM, r_au, energy)
        predicted_lx = luminosity(
            predicted_q, r_au, END_MEMBER_SPEED_KM_S, energy_flux=measured_h
        )
        comparison = {
            "b8_production_rate_s-1": predicted_q,
            "b12_luminosity_erg_s": predicted_lx,
            "observed_over_prediction": wegmann_lx / predicted_lx,
            "prediction_over_observed": predicted_lx / wegmann_lx,
        }
        if name in xmm_rates:
            comparison["b8_over_xmm_measured_rate"] = predicted_q / xmm_rates[name]
        end_members[name] = comparison
        rows.append(
            {
                "kind": "end_member",
                "name": name,
                "method": "Cabot B2-B9/B11-B12; R=1.3 km; v=0.37 km s-1",
                "radius_arcsec": wegmann_arcsec,
                "radius_km": WEGMANN_APERTURE_KM,
                "measured_luminosity_erg_s": wegmann_lx,
                "b8_production_rate_s-1": predicted_q,
                "predicted_luminosity_erg_s": predicted_lx,
                **{
                    key: comparison[key]
                    for key in ("observed_over_prediction", "prediction_over_observed")
                },
            }
        )

    radius_km = physical_radius_km(THIN_RADIUS_ARCMIN * 60, distance)
    surface_flux = model_surface_flux(task5, task3, THIN_RADIUS_ARCMIN)
    h2_h = ion_energy_flux(cabot_phi, mean_energy)
    thin = {}
    for speed in (1.0, 1.2):
        q_thin = thin_rate(surface_flux, radius_km * 1e5, speed, h2_h, H2_SIGMA_CM2)
        thin[f"{speed:.1f}_km_s"] = q_thin
    optical_depth = thin_tau(thin["1.0_km_s"], radius_km * 1e5, 1.0, H2_SIGMA_CM2)

    result = {
        "heliocentric_distance_au": r_au,
        "observer_distance_au": distance,
        "parent_rates": {
            "reference_distance_au": REFERENCE_R_AU,
            "reference_epoch_utc": REFERENCE_EPOCH_UTC,
            "reference_sum_s-1": reference_sum,
            "xmm_s-1": xmm_rates,
            "xmm_sum_s-1": xmm_sum,
        },
        "cabot_baseline": {
            "proton_density_1au_cm-3": PROTON_DENSITY_1AU_CM3,
            "solar_wind_speed_cm_s": SOLAR_WIND_SPEED_CM_S,
            "heavy_ion_fraction": HEAVY_ION_FRACTION,
            "photons_per_ion": PHOTONS_PER_ION,
            "photon_energy_eV": PHOTON_ENERGY_EV,
            "cross_section_cm2": SIGMA_CM2,
            "neutral_scale_length_cm": LAMBDA_CM,
            "kappa": KAPPA,
        },
        "heavy_ion": {
            "cabot_number_flux_cm-2_s-1": cabot_phi,
            "cabot_energy_flux_erg_cm-2_s-1": cabot_h,
            "task06_prior_midpoint_cm-2_s-1": task6_phi,
            "task06_prior_bounds_cm-2_s-1": task6_bounds,
            "task06_contemporaneous": bool(task6_prior["contemporaneous"]),
            "task06_energy_eV": mean_energy,
        },
        "wegmann": {
            "aperture_km": WEGMANN_APERTURE_KM,
            "aperture_arcsec": wegmann_arcsec,
            "luminosity_erg_s": wegmann_lx,
            "luminosity_log10_sigma": lx_sigma,
            "coefficient_s2_cm2": WEGMANN_C_S2_CM2,
            "coefficient_log10_sigma": c_sigma,
            "calibration_points": len(ratios),
            "measured_inputs": measured,
            "estimates": estimates,
        },
        "end_member_assumptions": {
            "nucleus_radius_km": END_MEMBER_RADIUS_KM,
            "outflow_speed_km_s": END_MEMBER_SPEED_KM_S,
            "energy_flux_erg_cm-2_s-1": measured_h,
            "aperture_km": WEGMANN_APERTURE_KM,
            "comparison_scope": "total B12 model versus aperture measurement",
            "luminosity_ratio_log10_sigma": math.hypot(
                lx_sigma, c_sigma, measured_sigma
            ),
        },
        "apertures": {
            name: {
                key: value for key, value in row.items() if key not in ("kind", "name")
            }
            for name, row in apertures.items()
        },
        "cabot_370_slow_speed_check": speed_check,
        "end_members": end_members,
        "h2_thin_check": {
            "radius_arcmin": THIN_RADIUS_ARCMIN,
            "radius_km": radius_km,
            "surface_flux_erg_cm-2_s-1_arcmin-2": surface_flux,
            "cross_section_cm2": H2_SIGMA_CM2,
            "cabot_ion_number_flux_cm-2_s-1": cabot_phi,
            "task03_mean_photon_energy_eV": mean_energy,
            "ion_energy_flux_erg_cm-2_s-1": h2_h,
            "optical_depth": optical_depth,
            "production_rate_s-1": thin,
        },
    }
    return result, rows


def write_outputs(result, rows, output=OUTPUT):
    columns = (
        "kind",
        "name",
        "method",
        "radius_arcsec",
        "radius_km",
        "measured_luminosity_erg_s",
        "inferred_production_rate_s-1",
        "task06_prior_inferred_rate_s-1",
        "xmm_parent_rate_s-1",
        "missing_rate_s-1",
        "missing_fraction",
        "b8_production_rate_s-1",
        "predicted_luminosity_erg_s",
        "observed_over_prediction",
        "prediction_over_observed",
    )
    output.mkdir(exist_ok=True)
    (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    with (output / "scaling.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def require_current_profile(task5=TASK5, task3=TASK3, morphology=MORPHOLOGY):
    inputs = (task3, morphology / "result.json", morphology / "image.fits")
    missing = [path for path in (task5, *inputs) if not path.is_file()]
    if missing:
        raise RuntimeError(f"missing prerequisite: {missing[0]}")
    if max(path.stat().st_mtime_ns for path in inputs) > task5.stat().st_mtime_ns:
        raise RuntimeError("Task 05 profile predates its calibration inputs")


def main():
    for name in FINAL_PRODUCTS:
        (OUTPUT / name).unlink(missing_ok=True)
    require_current_profile()
    task5, task3, task6 = (
        json.loads(path.read_text()) for path in (TASK5, TASK3, TASK6)
    )
    result, rows = calculate(task5, task3, task6)
    write_outputs(result, rows)


if __name__ == "__main__":
    main()
