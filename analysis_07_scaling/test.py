import csv
import json
import math
import os
import tempfile
from pathlib import Path

import numpy as np

from analysis_05_radial_profile.run import aperture_rows

from .model import (
    ARCMIN2_PER_SR,
    EV_TO_ERG,
    ICE_PROPERTIES,
    KAPPA,
    SIGMA_CM2,
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
    thick_radius_cm,
    thin_rate,
    thin_tau,
)
from .run import (
    FINAL_PRODUCTS,
    HERE,
    H2_SIGMA_CM2,
    SCIENCE_APERTURE,
    SCIENCE_NAME,
    SHARED,
    WEGMANN_APERTURE_KM,
    aperture_luminosity,
    calculate,
    parent_rates,
    physical_radius_km,
    read_calibration,
    require_current_profile,
    write_outputs,
)


def close(left, right, tolerance=1e-10):
    assert math.isclose(left, right, rel_tol=tolerance)


def fixtures():
    task5 = {
        "band_eV": SHARED["soft_band_ev"],
        "epic": {
            "slope": 0.9956183898388086,
            "normalization_ct_s_arcmin2": 0.009399913547404288,
            "broken": {
                "normalization_at_break_ct_s_arcmin2": 0.014709138249849813,
                "alpha_inner": 0.3463125079853897,
                "alpha_outer": 0.988526521467759,
                "break_radius_arcmin": 0.6282113792932272,
            },
        },
    }
    task3 = {
        "heliocentric_distance_au": 1.8688,
        "fit_band_keV": [value / 1000 for value in SHARED["soft_band_ev"]],
        "distance_au": 1.88252630197118,
        "mean_photon_energy_eV": 502.2,
        "mos2_ecf_erg_cm2_per_count": 8.672299198621948e-12,
        "mos2_ecf_sigma_ln": 0.09594655763290505,
        "regions": {"inner": {"luminosity_erg_s": 3.0e16}},
    }
    task6 = {
        "heavy_ion_prior": {**SHARED["heavy_ion_prior"], "lower": 2e3, "upper": 1.6e5},
        "proton_flux": {
            "stereo_a": {"median": 8.34e7},
            "solar_orbiter": {"median": 1.40e8},
            "wind": {"median": 6.02e7},
        },
        "ion_fraction": {
            "ace_swics": {"f_16_50_84": [1.386e-4, 2.492e-4, 4.396e-4]},
            "solar_orbiter_his": {"f_16_50_84": [1.131e-4, 2.020e-4, 3.554e-4]},
        },
    }
    return task5, task3, task6


def test_cabot_equations():
    water = molecular_energy_kj(*ICE_PROPERTIES["H2O"])
    close(water, 9.328629886495437e-23)
    q = production_rate(2.4, 1.0, water)
    close(q, 1.5e29)
    close(ion_number_flux_cm2_s(1.0), 2.8e5)
    close(cabot_energy_flux(1.0), 2.46735201636e-4)

    linear = SIGMA_CM2 * q / (4 * 1e5)
    radius = thick_radius_cm(q)
    assert 0 < radius < linear
    near_infinite = thick_radius_cm(q, length_cm=1e30)
    close(near_infinite, linear, 1e-8)
    lx = luminosity(q, 1.0)
    close(lx, math.pi * (KAPPA * radius) ** 2 * cabot_energy_flux(1.0))
    close(invert_luminosity(lx, 1.0), q, 1e-9)


def test_wegmann_and_thin_limits():
    h = cabot_energy_flux(1.8688)
    q = empirical_rate(6.128e15, h)
    close(q, 9.313365487643337e28)
    close(1e-38 * h * q * q, 6.128e15)

    radius = 4.9e10
    energy_flux = ion_energy_flux(8e4, 502.2)
    expected_q = 1.4e29
    surface_sr = energy_flux * H2_SIGMA_CM2 * expected_q / (16 * math.pi * 1e5 * radius)
    surface_arcmin = surface_sr / ARCMIN2_PER_SR
    recovered = thin_rate(surface_arcmin, radius, 1.0, energy_flux, H2_SIGMA_CM2)
    close(recovered, expected_q)
    close(
        thin_tau(recovered, radius, 1.0, H2_SIGMA_CM2),
        H2_SIGMA_CM2 * expected_q / (4e5 * radius),
    )


def test_calibration_scatter():
    ratios = read_calibration()
    assert len(ratios) == 48
    close(ratios.min(), 0.33351311373456954, 1e-9)
    close(ratios.max(), 4.651661408201318, 1e-9)
    assert 0.35 < coefficient_scatter_dex(ratios) < 0.36
    close(coefficient_scatter_dex([10.0, 0.1]), 1.0)
    assert coefficient_scatter_dex(np.ones(5)) == 0


def test_wegmann_estimates():
    task5, task3, task6 = fixtures()
    broken = task5["epic"]["broken"]
    theta = [
        math.log(broken["normalization_at_break_ct_s_arcmin2"]),
        broken["alpha_inner"],
        broken["alpha_outer"],
        math.log(broken["break_radius_arcmin"]),
    ]
    reference = aperture_rows(
        theta, task3["distance_au"], task3["mos2_ecf_erg_cm2_per_count"]
    )[0]
    close(aperture_luminosity(task5, task3, reference[2])[1], reference[5], 1e-9)

    result, rows = calculate(task5, task3, task6)
    wegmann = result["wegmann"]
    assert wegmann["aperture_km"] == WEGMANN_APERTURE_KM == 1.1e5
    assert wegmann["calibration_points"] == 48
    assert 0.35 < wegmann["coefficient_log10_sigma"] < 0.36
    close(wegmann["luminosity_log10_sigma"], 0.09594655763290505 / math.log(10))
    estimates = wegmann["estimates"]
    fast, scaled, measured = (
        estimates[name]
        for name in ("w04_fast_wind", "w04_fast_wind_scaled", "measured")
    )
    close(fast["energy_flux_erg_cm-2_s-1"], WEGMANN_FAST_H_ERG_CM2_S)
    close(
        scaled["inferred_production_rate_s-1"] / fast["inferred_production_rate_s-1"],
        task3["heliocentric_distance_au"],
    )
    inputs = wegmann["measured_inputs"]
    assert 0.18 < inputs["proton_flux_log10_sigma"] < 0.19
    assert 0.24 < inputs["ion_fraction_log10_sigma"] < 0.26
    close(
        measured["energy_flux_erg_cm-2_s-1"],
        inputs["ion_fraction"]
        * inputs["proton_flux_cm-2_s-1"]
        * task3["mean_photon_energy_eV"]
        * EV_TO_ERG,
    )
    for estimate in (fast, scaled, measured):
        low, high = estimate["inferred_production_rate_range_s-1"]
        q = estimate["inferred_production_rate_s-1"]
        close(low * high, q * q)
        spread = math.hypot(
            wegmann["luminosity_log10_sigma"],
            wegmann["coefficient_log10_sigma"],
            estimate["energy_flux_log10_sigma"],
        )
        close(high / low, 10**spread)
    assert [row["name"] for row in rows[:3]] == list(estimates)
    assert all(row["method"] == "Wegmann Eq. 30" for row in rows[:3])


def test_paper_calculation():
    task5, task3, task6 = fixtures()
    result, rows = calculate(task5, task3, task6)
    _, reference, xmm = parent_rates(task3["heliocentric_distance_au"])
    close(reference, 2.46e28)
    close(xmm, 3.156043338888213e28)
    assert (
        result["parent_rates"]["reference_epoch_utc"]
        == SHARED["spherex"]["reference_epoch_utc"]
    )
    assert (
        2.0e29
        < result["apertures"][SCIENCE_NAME]["inferred_production_rate_s-1"]
        < 2.2e29
    )
    assert result["apertures"][SCIENCE_NAME]["missing_fraction"] > 0.8
    speed_check = result["cabot_370_slow_speed_check"]
    close(speed_check["speed_km_s"], 0.37)
    close(
        speed_check["inferred_production_rate_s-1"],
        0.37 * result["apertures"][SCIENCE_NAME]["inferred_production_rate_s-1"],
    )
    assert speed_check["missing_fraction"] > 0.5
    assert speed_check["required_heavy_ion_flux_multiple_for_xmm_parents"] > 5
    assert (
        38
        < result["apertures"][SCIENCE_NAME][
            "required_heavy_ion_flux_multiple_for_xmm_parents"
        ]
        < 39
    )
    assert not result["heavy_ion"]["task06_contemporaneous"]
    assert result["cabot_baseline"] == {
        "proton_density_1au_cm-3": 7.0,
        "solar_wind_speed_cm_s": 4e7,
        "heavy_ion_fraction": 1e-3,
        "photons_per_ion": 1.0,
        "photon_energy_eV": 550.0,
        "cross_section_cm2": 3e-15,
        "neutral_scale_length_cm": 5e10,
        "kappa": 8.1,
    }

    end = result["end_members"]
    assert 9 < end["CO"]["b8_over_xmm_measured_rate"] < 11
    thin = result["h2_thin_check"]
    rates = thin["production_rate_s-1"]
    close(rates["1.0_km_s"], 1.4618168793339194e29)
    close(thin["optical_depth"], 0.031524961453890665)
    close(rates["1.2_km_s"] / rates["1.0_km_s"], 1.2)
    close(
        thin["cabot_ion_number_flux_cm-2_s-1"],
        result["heavy_ion"]["cabot_number_flux_cm-2_s-1"],
    )
    close(thin["task03_mean_photon_energy_eV"], task3["mean_photon_energy_eV"])
    changed = json.loads(json.dumps(task5))
    changed["epic"]["broken"]["normalization_at_break_ct_s_arcmin2"] *= 2
    repeated, _ = calculate(changed, task3, task6)
    close(
        repeated["h2_thin_check"]["production_rate_s-1"]["1.0_km_s"], rates["1.0_km_s"]
    )
    close(
        result["h2_thin_check"]["radius_km"],
        physical_radius_km(360, task3["distance_au"]),
    )
    assert rows[3]["name"] == SCIENCE_NAME
    assert SCIENCE_APERTURE == 370
    assert rows[3]["method"] == "Cabot B2-B7/B11-B12; v=1 km s-1"
    assert all(
        row["method"] == "Cabot B2-B9/B11-B12; R=1.3 km; v=0.37 km s-1"
        for row in rows[4:]
    )

    with tempfile.TemporaryDirectory(dir=HERE) as folder:
        output = Path(folder)
        write_outputs(result, rows, output)
        assert json.loads((output / "result.json").read_text()) == result
        with (output / "scaling.tsv").open(newline="") as stream:
            table = list(csv.DictReader(stream, delimiter="\t"))
        assert len(table) == len(rows)
        assert tuple(table[0]) == (
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
        assert FINAL_PRODUCTS == ("result.json", "scaling.tsv")


def test_end_member_inputs():
    task5, task3, task6 = fixtures()
    result, rows = calculate(task5, task3, task6)
    end = result["end_members"]
    assumptions = result["end_member_assumptions"]
    wegmann = result["wegmann"]
    measured = wegmann["estimates"]["measured"]
    assert assumptions["aperture_km"] == 110000
    close(assumptions["energy_flux_erg_cm-2_s-1"], measured["energy_flux_erg_cm-2_s-1"])
    low, high = measured["inferred_production_rate_range_s-1"]
    close(assumptions["luminosity_ratio_log10_sigma"], math.log10(high / low))
    for row in rows:
        if row["kind"] == "end_member":
            assert row["radius_km"] == 110000
            close(row["radius_arcsec"], wegmann["aperture_arcsec"])
            close(row["measured_luminosity_erg_s"], wegmann["luminosity_erg_s"])

    task3["regions"]["inner"]["luminosity_erg_s"] *= 2
    repeated, _ = calculate(task5, task3, task6)
    assert repeated["end_members"] == end
    assert repeated["end_member_assumptions"] == assumptions

    for monitor in task6["proton_flux"].values():
        monitor["median"] *= 2
    repeated, _ = calculate(task5, task3, task6)
    for name, original in end.items():
        updated = repeated["end_members"][name]
        close(updated["b8_production_rate_s-1"], original["b8_production_rate_s-1"])
        close(updated["b12_luminosity_erg_s"], 2 * original["b12_luminosity_erg_s"])
        close(
            updated["observed_over_prediction"],
            original["observed_over_prediction"] / 2,
        )
        close(
            updated["prediction_over_observed"],
            2 * original["prediction_over_observed"],
        )

    task5["epic"]["broken"]["normalization_at_break_ct_s_arcmin2"] *= 2
    repeated, _ = calculate(task5, task3, task6)
    for name, original in end.items():
        updated = repeated["end_members"][name]
        close(updated["observed_over_prediction"], original["observed_over_prediction"])
        close(updated["prediction_over_observed"], original["prediction_over_observed"])


def test_profile_freshness():
    with tempfile.TemporaryDirectory(dir=HERE) as folder:
        task5, task3 = (Path(folder) / name for name in ("task5", "task3"))
        morphology = Path(folder) / "morphology"
        morphology.mkdir()
        morphology_inputs = (morphology / "result.json", morphology / "image.fits")
        task5.touch()
        task3.touch()
        for path in morphology_inputs:
            path.touch()
        newer = (
            max(path.stat().st_mtime_ns for path in (task3, *morphology_inputs))
            + 2_000_000_000
        )
        os.utime(task5, ns=(newer, newer))
        require_current_profile(task5, task3, morphology)
        newer += 2_000_000_000
        os.utime(morphology_inputs[0], ns=(newer, newer))
        try:
            require_current_profile(task5, task3, morphology)
        except RuntimeError as error:
            assert "predates" in str(error)
        else:
            raise AssertionError("stale Task 05 profile accepted")
        task3.unlink()
        try:
            require_current_profile(task5, task3, morphology)
        except RuntimeError as error:
            assert "missing prerequisite" in str(error)
        else:
            raise AssertionError("missing Task 03 result accepted")


def main():
    test_cabot_equations()
    test_wegmann_and_thin_limits()
    test_calibration_scatter()
    test_wegmann_estimates()
    test_paper_calculation()
    test_end_member_inputs()
    test_profile_freshness()


if __name__ == "__main__":
    main()
