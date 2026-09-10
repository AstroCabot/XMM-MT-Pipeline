import csv
import os
import sys

import numpy as np

from _common import DATA, OUT, REPRO, read_json, restore_env, stage_dir, write_json

sys.path.insert(0, str(REPRO / "analysis_08_chemistry"))

import run as run8
import figure
import inference
import products
import psf_operator

restore_env()

DEMO_MCMC = {
    "walkers": 48,
    "burn_steps": 200,
    "retained_steps": 600,
    "processes": max(2, min(8, os.cpu_count() or 2)),
    "seed": 20260910,
}
DEMO_MODEL = {
    "shells": [80, 1.3, 3.0e9],
    "source_grid_arcmin": [0.01, 30.0, 0.05],
}
DEMO_QUADRATURE = {
    "granddaughter_shell_stride": 8,
    "granddaughter_radial_stride": 8,
    "psf_azimuth_order": 192,
}


def profile_figure(table, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with open(table, newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    lo = np.array([float(row["r_lo_arcsec"]) for row in rows]) / 60
    hi = np.array([float(row["r_hi_arcsec"]) for row in rows]) / 60
    center = (2 / 3) * (hi**3 - lo**3) / (hi**2 - lo**2)
    data = np.array([float(row["data_ct_s_arcmin2"]) for row in rows])
    error = np.array([float(row["total_error_ct_s_arcmin2"]) for row in rows])
    model = np.array([float(row["prediction_ct_s_arcmin2"]) for row in rows])
    edges = np.r_[center[0] / 2, hi]
    fig, ax = plt.subplots(figsize=(5.4, 3.7), layout="constrained")
    ax.stairs(model, edges, color="#c51b7d", lw=1.5, label="posterior prediction (MAP)")
    ax.errorbar(
        center,
        data,
        xerr=(center - lo, hi - center),
        yerr=error,
        fmt="o",
        ms=3.5,
        lw=0.9,
        color="#202020",
        capsize=0,
        label="background-subtracted annuli",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Radius (arcmin)")
    ax.set_ylabel(r"Surface brightness (count s$^{-1}$ arcmin$^{-2}$)")
    ax.legend(frameon=False, fontsize=9)
    fig.savefig(path, dpi=250)
    plt.close(fig)


def main():
    output = stage_dir(8)
    run8.PROFILE = OUT / "05"
    run8.SPECTRUM = OUT / "03" / "result.json"
    run8.WIND = OUT / "06" / "result.json"
    run8.GEOMETRY = OUT / "04" / "result.json"
    run8.MORPHOLOGY_IMAGE = OUT / "04" / "image.fits"
    run8.OUTPUT = output
    run8.HERE = DATA / "chemistry"
    products.OUTPUT = output
    products.PROFILE = OUT / "05"
    figure.draw = lambda *args, **kwargs: None

    run8.CONFIG["mcmc"] = dict(DEMO_MCMC)
    run8.CONFIG.update(DEMO_MODEL)
    run8.CONFIG["quadrature"].update(DEMO_QUADRATURE)

    inputs = run8.load_inputs()
    with np.load(DATA / "chemistry" / "psf.npz", allow_pickle=False) as archive:
        if not np.allclose(archive["annulus_edges_arcsec"], inputs["edges"]):
            raise RuntimeError("archived PSF annuli differ from the demo annuli")
        source = np.arange(*run8.CONFIG["source_grid_arcmin"])
        operator = psf_operator.radial_operator(
            archive["radius_arcsec"],
            archive["density_arcsec2"],
            source,
            archive["annulus_edges_arcsec"],
            run8.CONFIG["quadrature"]["psf_radial_order"],
            run8.CONFIG["quadrature"]["psf_azimuth_order"],
            archive["exposure_radius_arcmin"],
            archive["exposure_mean_s"],
        )

    posterior, reactions = run8.build_model(inputs, operator)
    settings = run8.CONFIG["mcmc"]
    rng = np.random.default_rng(settings["seed"])
    start = inference.draw_walkers(posterior, settings["walkers"], rng)
    chain, logp, logl, acceptance = inference.run_sampler(
        posterior,
        start,
        settings["burn_steps"],
        settings["retained_steps"],
        settings["processes"],
        settings["seed"],
    )
    products.write(inputs, posterior, reactions, chain, logp, logl, acceptance)
    profile_figure(output / "profile.tsv", output / "profile_fit.png")

    posterior_summary = read_json(output / "result.json")["posterior"]
    line = ", ".join(
        f"{name} {posterior_summary[name]['median']:.2e}"
        for name in ("Q_H2O", "Q_CO2", "Q_CO", "Q_H2")
        if name in posterior_summary
    )
    print(f"stage 08: acceptance {acceptance:.3f}; posterior medians {line}")
    if "H2_fraction" in posterior_summary:
        h2 = posterior_summary["H2_fraction"]
        print(
            f"stage 08: H2 fraction median {h2['median']:.3f}"
            f" (q16 {h2['q16']:.3f}, q84 {h2['q84']:.3f})"
        )


if __name__ == "__main__":
    main()
