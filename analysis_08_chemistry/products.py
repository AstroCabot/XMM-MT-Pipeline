import csv
import json
from pathlib import Path

import chemistry
import diagnostic
import figure
import inference
import numpy as np
from config import CONFIG
from transport import tau_one_radius

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "outputs"
PROFILE = HERE.parent / "analysis_05_radial_profile" / "outputs"
FINAL_PRODUCTS = (
    "result.json",
    "chain.npz",
    "posterior.tsv",
    "profile.tsv",
    "intrinsic.tsv",
    "chemistry.tsv",
    "figure3.png",
    "figure3.pdf",
    "corner.png",
)


def write_tsv(path, columns, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(columns)
        writer.writerows(rows)


def quantiles(values, weights=None):
    values = np.asarray(values, float)
    if weights is None:
        result = np.quantile(values, (0.16, 0.5, 0.84))
    else:
        order = np.argsort(values)
        cumulative = np.cumsum(np.asarray(weights)[order])
        result = np.interp(
            (0.16, 0.5, 0.84), cumulative / cumulative[-1], values[order]
        )
    return [float(value) for value in result]


def jitter_checks(values):
    base = CONFIG["priors"]["jitter_sigma"]
    jitter = values[:, inference.PHYSICAL_NAMES.index("jitter")]
    q_h2 = values[:, inference.PHYSICAL_NAMES.index("Q_H2")]
    fraction = q_h2 / values[:, :4].sum(axis=1)
    output = {}
    for sigma in (0.05, 0.20):
        log_weight = np.log(base / sigma) - 0.5 * jitter**2 * (sigma**-2 - base**-2)
        weight = np.exp(log_weight - log_weight.max())
        output[f"sigma_{sigma:.2f}"] = {
            "Q_H2_median_s-1": quantiles(q_h2, weight)[1],
            "H2_fraction_median": quantiles(fraction, weight)[1],
            "effective_sample_size": float(weight.sum() ** 2 / np.sum(weight**2)),
        }
    return output


def sensitivity_checks(posterior, theta, reference):
    def change(trial):
        scale = np.maximum(np.abs(reference), np.max(np.abs(reference)) * 1e-12)
        return float(np.max(np.abs(trial - reference) / scale))

    fixed = posterior.config["fixed_velocity_km_s"]
    original_speed = fixed.copy()
    velocity_change = []
    try:
        for name, bounds in (
            ("heavy_daughter", (0.30, 3.07)),
            ("photolytic_H2", (0.30, 22.37)),
        ):
            for value in bounds:
                fixed[name] = value
                velocity_change.append(change(posterior.predict(theta)[2]))
            fixed[name] = original_speed[name]
    finally:
        fixed.update(original_speed)

    parents = {"H2O", "CO2", "CO", "H2"}
    proxy = {
        name: row["sigma_cm2"]
        for name, row in posterior.transport.network.items()
        if name not in parents
    }
    cross_section_change = []
    try:
        for value in (4.40e-15, 5.41e-15):
            for name in proxy:
                posterior.transport.network[name]["sigma_cm2"] = value
            cross_section_change.append(change(posterior.predict(theta)[2]))
    finally:
        for name, value in proxy.items():
            posterior.transport.network[name]["sigma_cm2"] = value
    return {
        "fixed_minor_velocity_ranges_km_s": {
            "heavy_daughter": [0.30, 3.07],
            "photolytic_H2": [0.30, 22.37],
        },
        "fixed_minor_velocity_maximum_profile_fractional_change": max(velocity_change),
        "daughter_proxy_cross_section_range_cm2": [4.40e-15, 5.41e-15],
        "daughter_proxy_maximum_profile_fractional_change": max(cross_section_change),
    }


def point_record(posterior, theta, log_probability, log_likelihood):
    physical = posterior.predict(theta)[-1]
    total = sum(physical[:4])
    return {
        "log_probability": float(log_probability),
        "log_likelihood": float(log_likelihood),
        "Q_total_s-1": float(total),
        "Q_H2_s-1": float(physical[3]),
        "H2_fraction": float(physical[3] / total),
    }


def orbit_output():
    return {
        name: CONFIG["orbit"][name]
        for name in ("perihelion_q_au", "eccentricity", "perihelion_jd_tdb")
    }


def mode_selection(theta, log_likelihood):
    ratio = theta[:, inference.COORDS.index("log10_H2_ratio")]
    high, low = ratio >= -0.5, ratio <= -1.5
    if not high.any():
        raise RuntimeError("chain did not sample the required high-H2 mode")

    def best(mask):
        return int(np.flatnonzero(mask)[np.argmax(log_likelihood[mask])])

    return high, best(high), best(low) if low.any() else None


def write(inputs, posterior, reactions, chain, logp, logl, acceptance_fraction):
    flat = chain.reshape(-1, chain.shape[-1])
    flat_lp, flat_ll = logp.ravel(), logl.ravel()
    ratio_index = inference.COORDS.index("log10_H2_ratio")
    high, high_ml, low_ml = mode_selection(flat, flat_ll)
    low = flat[:, ratio_index] <= -1.5
    overall_ml = int(np.argmax(flat_ll))
    high_map = int(np.flatnonzero(high)[np.argmax(flat_lp[high])])
    low_map = (
        None if low_ml is None else int(np.flatnonzero(low)[np.argmax(flat_lp[low])])
    )
    best_theta = flat[high_map]

    stride = max(1, chain.shape[0] // 500)
    theta = chain[::stride].reshape(-1, chain.shape[-1])
    lp, ll = logp[::stride].ravel(), logl[::stride].ravel()
    physical = np.array(
        [inference.unpack(row, CONFIG, posterior.transport.orbit)[-1] for row in theta]
    )
    q = physical[:, :4]
    values = np.c_[physical, q.sum(axis=1), q[:, 3] / q.sum(axis=1)]
    map_physical = inference.unpack(best_theta, CONFIG, posterior.transport.orbit)[-1]
    map_q = map_physical[:4]
    map_values = np.r_[map_physical, map_q.sum(), map_q[3] / map_q.sum()]
    names = inference.PHYSICAL_NAMES + ("Q_total", "H2_fraction")
    units = {name: "1" for name in names}
    units.update(
        {name: "s^-1" for name in ("Q_H2O", "Q_CO2", "Q_CO", "Q_H2", "Q_total")}
    )
    units.update({name: "km s^-1" for name in names if name.startswith("v_")})
    units.update(phi_heavy="cm^-2 s^-1", offset="count s^-1 arcmin^-2")
    write_tsv(
        OUTPUT / "posterior.tsv",
        ("parameter", "unit", "map", "q16", "median", "q84"),
        [
            [name, units.get(name, "1"), map_values[i], *quantiles(values[:, i])]
            for i, name in enumerate(names)
        ],
    )

    prediction, variance, model, source, tau, _ = posterior.predict(best_theta)
    write_tsv(
        OUTPUT / "profile.tsv",
        (
            "r_lo_arcsec",
            "r_hi_arcsec",
            "data_ct_s_arcmin2",
            "total_error_ct_s_arcmin2",
            "physical_model_ct_s_arcmin2",
            "prediction_ct_s_arcmin2",
            "residual_ct_s_arcmin2",
        ),
        [
            [
                *inputs["edges"][i : i + 2],
                inputs["data"][i],
                np.sqrt(variance[i]),
                model[i],
                prediction[i],
                inputs["data"][i] - prediction[i],
            ]
            for i in range(10)
        ],
    )
    write_tsv(
        OUTPUT / "intrinsic.tsv",
        (
            "radius_arcmin",
            "energy_flux_erg_cm2_s_arcmin2",
            "count_rate_ct_s_arcmin2",
            "tau",
        ),
        zip(posterior.transport.source_radius, source, source / posterior.ecf, tau),
    )
    totals = {
        name: sum(row["rate_s-1"] for row in reactions if row["parent"] == name)
        for name in chemistry.FAMILIES
    }
    columns = ("parent", "products", "rate_s-1", "branch") + chemistry.FAMILIES
    write_tsv(
        OUTPUT / "chemistry.tsv",
        columns,
        [
            [
                row["parent"],
                row["products"],
                row["rate_s-1"],
                row["rate_s-1"] / totals[row["parent"]],
            ]
            + [row[name] for name in chemistry.FAMILIES]
            for row in reactions
        ],
    )

    keep = physical
    best = [("overall_ml", overall_ml), ("high_ml", high_ml)]
    if low_ml is not None:
        best.append(("low_ml", low_ml))
    best.append(("high_map", high_map))
    if low_map is not None:
        best.append(("low_map", low_map))
    best_index = np.array([index for _, index in best])
    np.savez_compressed(
        OUTPUT / "chain.npz",
        samples=theta.astype("f4"),
        coordinate_names=np.asarray(inference.COORDS),
        log_probability=lp.astype("f4"),
        log_likelihood=ll.astype("f4"),
        physical_draws=keep.astype("f4"),
        physical_names=np.asarray(inference.PHYSICAL_NAMES),
        best_labels=np.asarray([label for label, _ in best]),
        best_coordinates=flat[best_index],
        best_log_probability=flat_lp[best_index],
        best_log_likelihood=flat_ll[best_index],
    )
    summary = {
        name: dict(zip(("q16", "median", "q84"), quantiles(values[:, i])))
        for i, name in enumerate(names)
        if name.startswith("Q_") or name == "H2_fraction"
    }
    radius, status = tau_one_radius(posterior.transport.source_radius, tau)
    result = {
        "xmm_heliocentric_distance_au": posterior.transport.orbit.current_r_au,
        "xmm_observation_jd_tdb": CONFIG["orbit"]["observation_jd_tdb"],
        "orbit": orbit_output(),
        "f107p_sfu": inputs["f107p"],
        "photons_per_captured_ion": CONFIG["photons_per_captured_ion"],
        "mean_acceptance_fraction": float(acceptance_fraction),
        "background": {"model": "continuous broken-law B", **inputs["background"]},
        "posterior": summary,
        "maximum_likelihood": point_record(
            posterior, flat[overall_ml], flat_lp[overall_ml], flat_ll[overall_ml]
        ),
        "maximum_likelihood_high": point_record(
            posterior, flat[high_ml], flat_lp[high_ml], flat_ll[high_ml]
        ),
        "maximum_a_posteriori_high": point_record(
            posterior, best_theta, flat_lp[high_map], flat_ll[high_map]
        ),
        "maximum_a_posteriori_low": (
            None
            if low_map is None
            else point_record(
                posterior, flat[low_map], flat_lp[low_map], flat_ll[low_map]
            )
        ),
        "sampled_delta_log_likelihood_high_minus_low": (
            None if low_ml is None else float(flat_ll[high_ml] - flat_ll[low_ml])
        ),
        "sampled_delta_log_probability_high_minus_low": (
            None if low_map is None else float(flat_lp[high_map] - flat_lp[low_map])
        ),
        "tau_one_radius_arcmin": radius,
        "tau_one_status": status,
        "jitter_prior_checks": jitter_checks(physical),
        "fixed_input_checks": sensitivity_checks(posterior, best_theta, model),
    }
    (OUTPUT / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    scale = map_physical[inference.PHYSICAL_NAMES.index("scale")]
    figure.draw(
        PROFILE / "profile.tsv",
        OUTPUT / "profile.tsv",
        OUTPUT / "intrinsic.tsv",
        background=inputs["background"]["mean_ct_s_arcmin2"],
        tau_radius=radius,
        scale=scale,
        km_per_arcmin=inputs["geometry"]["km_per_arcmin"],
        broken=inputs["broken"],
        output=OUTPUT,
        ecf=posterior.ecf,
    )
    diagnostic.draw(physical, lp, OUTPUT / "corner.png")
