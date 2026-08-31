import csv
import json
import math
import os
import tempfile
from pathlib import Path

from .model import (
    ARCMIN2_PER_SR,
    ICE_PROPERTIES,
    KAPPA,
    SIGMA_CM2,
    cabot_energy_flux,
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
    MATCHED_NAMES,
    N2_APERTURE,
    SCIENCE_APERTURE,
    SCIENCE_NAME,
    SHARED,
    calculate,
    parent_rates,
    physical_radius_km,
    require_current_profile,
    write_outputs,
)


def close(left, right, tolerance=1e-10):
    assert math.isclose(left, right, rel_tol=tolerance)


def fixtures():
    small, large = MATCHED_NAMES
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
        "apertures": {
            small: {
                "radius_arcsec": 87.8901,
                "radius_km": 120000.0,
                "luminosity_erg_s": 6.128e15,
            },
            large: {
                "radius_arcsec": 98.8764,
                "radius_km": 135000.0,
                "luminosity_erg_s": 7.0565e15,
            },
        },
    }
    task3 = {
        "heliocentric_distance_au": 1.8688,
        "fit_band_keV": [value / 1000 for value in SHARED["soft_band_ev"]],
        "distance_au": 1.88252630197118,
        "mean_photon_energy_eV": 502.2,
        "mos2_ecf_erg_cm2_per_count": 8.672299198621948e-12,
        "regions": {"inner": {"luminosity_erg_s": 3.0e16}},
    }
    task6 = {"heavy_ion_prior": SHARED["heavy_ion_prior"]}
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


def test_paper_calculation():
    task5, task3, task6 = fixtures()
    result, rows = calculate(task5, task3, task6)
    _, reference, xmm = parent_rates(task3["heliocentric_distance_au"])
    close(reference, 2.46e28)
    close(xmm, 3.1418250350221415e28)
    close(
        result["apertures"][MATCHED_NAMES[0]]["inferred_production_rate_s-1"],
        9.313365487643337e28,
    )
    close(
        result["apertures"][MATCHED_NAMES[1]]["inferred_production_rate_s-1"],
        9.994059738837967e28,
    )
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
    assert 32 < end["H2O"]["observed_over_prediction"] < 34
    assert 9 < end["CO2"]["observed_over_prediction"] < 10
    assert 38 < end["H2"]["prediction_over_observed"] < 40
    assert 9 < end["CO"]["b8_over_xmm_measured_rate"] < 11
    close(result["n2_alternative"]["Q_N2_over_Q_CO"], 7.89434563374471)
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
    assert [row["name"] for row in rows[:3]] == [*MATCHED_NAMES, SCIENCE_NAME]
    assert result["n2_alternative"]["aperture"] == N2_APERTURE
    assert SCIENCE_APERTURE == 370
    assert all(
        "Wegmann" in row["method"] and "approximate" in row["method"]
        for row in rows[:2]
    )
    assert rows[2]["method"] == "Cabot B2-B7/B11-B12; v=1 km s-1"
    assert all(
        row["method"] == "Cabot B2-B9/B11-B12; R=1.3 km; v=0.37 km s-1"
        for row in rows[3:]
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
    test_paper_calculation()
    test_profile_freshness()


if __name__ == "__main__":
    main()
