import json
import math
import tempfile
from pathlib import Path

import numpy as np
from physics import (
    MASS_U,
    U_KG,
    acceleration,
    calculate,
    capacity_rates,
    delayed_radius,
    ellipsoid_area,
    nucleus,
    production_draws,
    shell_volume,
    thermal_speed,
)
from run import load_draws, run

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PARAMETERS = json.loads((HERE / "parameters.json").read_text())


def close(value, expected, relative=1e-3):
    assert math.isclose(value, expected, rel_tol=relative)


def test_physics():
    m_h2 = MASS_U["H2"] * U_KG
    close(thermal_speed(m_h2, 100), 1024.84, 2e-5)
    close(ellipsoid_area(57.5, 55.5, 9.5), 21472.0, 2e-5)
    close(ellipsoid_area(2, 2, 2), 16 * math.pi, 1e-12)
    close(shell_volume(100, 1), 4 * math.pi / 3 * (100**3 - 99**3), 1e-12)
    radius, mass, result = nucleus(PARAMETERS)
    close(radius / 1000, 1.323, 5e-4)
    close(result["radius_km"]["sigma"], 0.210, 5e-3)
    close(mass, 4.85e12, 3e-3)


def test_spada():
    current = 1.868783666266
    delayed = delayed_radius(current, 1.3564, 6.1396, 7.49)
    close(delayed, 1.69377, 2e-5)
    close(acceleration(4.941e-8, current), 2.835e-7, 2e-3)
    close(acceleration(5.381e-8, delayed), 3.759e-7, 2e-3)


def synthetic_draws():
    names = [
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
    ]
    rows = [
        (2.0e28, 3.0e27, 8.0e27, 7.0e28, 0.32, 0.42, 0.52),
        (2.2e28, 3.2e27, 8.5e27, 6.8e28, 0.37, 0.47, 0.57),
        (1.8e28, 2.8e27, 7.5e27, 7.2e28, 0.42, 0.52, 0.62),
        (2.0e28, 3.0e27, 8.0e27, 1.0e27, 0.37, 0.47, 0.57),
    ]
    draws = np.ones((len(rows), len(names)))
    for i, row in enumerate(rows):
        for name, value in zip(
            ("Q_H2O", "Q_CO2", "Q_CO", "Q_H2", "v_H2O_km_s", "v_CO2_km_s", "v_CO_km_s"),
            row,
            strict=True,
        ):
            draws[i, names.index(name)] = value
    return draws, names


def test_production_draws_and_result():
    draws, names = synthetic_draws()
    q, keep = production_draws(draws, names)
    assert len(q["H2"]) == 4 and keep.all()
    invalid = draws.copy()
    invalid[0, names.index("Q_H2")] = 0
    try:
        production_draws(invalid, names)
    except ValueError as error:
        assert "invalid production draws" in str(error)
    else:
        raise AssertionError("invalid production draw accepted")
    task08 = {
        "xmm_heliocentric_distance_au": 1.868783666266,
        "xmm_observation_jd_tdb": 2461013.291968092,
        "orbit": {
            "perihelion_q_au": 1.3564,
            "eccentricity": 6.1396,
            "perihelion_jd_tdb": 2460977.995262848,
        },
    }
    result = calculate(draws, names, task08, PARAMETERS)
    recoil = result["radial_recoil"]
    assert recoil["headline"] == "symmetric"
    assert recoil["symmetric"]["delay_days"] == 0
    close(recoil["symmetric"]["required_force_N"]["value"], 1.374e6, 5e-3)
    close(recoil["delayed"]["required_force_N"]["value"], 1.822e6, 5e-3)
    close(
        recoil["symmetric"]["required_force_N"]["sigma"]
        / recoil["symmetric"]["required_force_N"]["value"],
        1.5 * 0.07 / 0.22,
    )
    oum = result["oumuamua"]
    assert oum["axes_albedo"] == PARAMETERS["oumuamua"]["axes_albedo"] == 0.1
    close(oum["H2_requirement_molecules_s"]["at_1.25_au"], 7.8e25, 0.02)
    close(oum["H2_requirement_molecules_s"]["at_xmm_distance"], 3.5e25, 0.02)
    published = oum["published_H2_range_at_1.25_au"]
    assert published["mass_loss_kg_s"] == [0.01, 0.8]
    close(published["molecules_s"][0], 3e24, 0.01)
    close(published["molecules_s"][1], 2.4e26, 0.01)
    ratio = oum["same_distance_rate_ratio"]
    assert ratio["p16"] < ratio["median"] < ratio["p84"]
    assert 1.7 < ratio["median"] < 2.2
    fractions = result["gas"]["supported_force_fraction_zeta1_central_mass"]
    assert (
        0 < fractions["H2"]["median_range"][0] < fractions["H2"]["median_range"][1] < 1
    )
    assert (
        fractions["H2"]["median_range"][1] < fractions["all_gas"]["median_range"][1] < 1
    )
    assert fractions["H2"]["p16_min"] < fractions["H2"]["median_range"][0]
    assert fractions["all_gas"]["median_range"][1] < fractions["all_gas"]["p84_max"]
    expected = np.zeros(np.count_nonzero(keep))
    for species in ("H2O", "CO2", "CO"):
        rates = draws[keep, names.index(f"Q_{species}")]
        speed = draws[keep, names.index(f"v_{species}_km_s")]
        expected += rates * MASS_U[species] * U_KG * speed * 1e3
    close(
        result["gas"]["momentum_N"]["known_parents"]["median"],
        np.median(expected),
        1e-12,
    )
    capacity = result["radiolytic_capacity"]
    assert (
        capacity["activity_interval_days"]
        == PARAMETERS["reservoir"]["activity_interval_days"]
    )
    for case in ("gas_erosion", "dust_inclusive_erosion"):
        values = capacity[case]["H2_molecules_s"]
        assert 0 < values[0] < values[1]
    task08["xmm_observation_jd_tdb"] = task08["orbit"]["perihelion_jd_tdb"] - 1
    try:
        calculate(draws, names, task08, PARAMETERS)
    except ValueError as error:
        assert "post-perihelion" in str(error)
    else:
        raise AssertionError("pre-perihelion recoil calculation accepted")


def test_capacity_regression():
    p = PARAMETERS["reservoir"]
    args = (
        p["radius_m"],
        p["activity_interval_days"],
        p["water_equivalent_density_kg_m3"],
        p["h2_per_h2o"],
    )
    cases = [
        (
            p["gas_erosion_m"],
            [6.309416252e27, 7.245227478e28],
        ),
        (
            p["dust_erosion_m"],
            [1.265555418e28, 3.596520196e29],
        ),
    ]
    for depths, expected in cases:
        result = capacity_rates(args[0], args[1], depths, args[2], args[3])
        for value, target in zip(result, expected, strict=True):
            close(value, target, 2e-9)


def test_io():
    draws, names = synthetic_draws()
    task08 = {
        "xmm_heliocentric_distance_au": 1.868783666266,
        "xmm_observation_jd_tdb": 2461013.291968092,
        "orbit": {
            "perihelion_q_au": 1.3564,
            "eccentricity": 6.1396,
            "perihelion_jd_tdb": 2460977.995262848,
        },
    }
    with tempfile.TemporaryDirectory(dir=HERE) as folder:
        folder = Path(folder)
        chain = folder / "chain.npz"
        np.savez_compressed(chain, physical_draws=draws, physical_names=np.array(names))
        loaded, loaded_names = load_draws(chain)
        assert loaded_names == names and np.array_equal(loaded, draws)
        source = folder / "result.json"
        source.write_text(json.dumps(task08))
        output = folder / "out/result.json"
        result = run(source, chain, HERE / "parameters.json", output)
        assert json.loads(output.read_text()) == result
        output.write_text("stale")
        try:
            run(folder / "missing.json", chain, HERE / "parameters.json", output)
        except FileNotFoundError:
            assert not output.exists()
        else:
            raise AssertionError("missing input accepted")
        before = source.read_text()
        try:
            run(source, chain, HERE / "parameters.json", source)
        except ValueError as error:
            assert "differ" in str(error) and source.read_text() == before
        else:
            raise AssertionError("input/output alias accepted")
        case_alias = source.with_name(source.name.upper())
        if case_alias.exists():
            try:
                run(source, chain, HERE / "parameters.json", case_alias)
            except ValueError as error:
                assert "differ" in str(error) and source.read_text() == before
            else:
                raise AssertionError("case-variant input/output alias accepted")
        try:
            run(source, chain, HERE / "parameters.json", ROOT.parent / "forbidden.json")
        except ValueError as error:
            assert "inside analysis_repro" in str(error)
        else:
            raise AssertionError("outside output accepted")


def main():
    test_physics()
    test_spada()
    test_production_draws_and_result()
    test_capacity_regression()
    test_io()


if __name__ == "__main__":
    main()
