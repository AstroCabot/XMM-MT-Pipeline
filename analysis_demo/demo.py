#!/usr/bin/env python
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

DEMO = Path(__file__).resolve().parent
REPRO = DEMO.parent
OUT = DEMO / "outputs"
STAGES = (
    ("02", "PPS luminosity geometry"),
    ("03", "PN spectral fit (PyXspec + ACX2)"),
    ("04", "comet-frame morphology"),
    ("05", "radial profile"),
    ("06", "solar wind"),
    ("07", "production-rate scaling"),
    ("08", "chemistry posterior (reduced MCMC)"),
    ("09", "recoil and radiolytic capacity"),
)
FIGURES = (
    ("04/figure1.png", "figure1_comet_image.png"),
    ("03/figure2.png", "figure2_spectrum.png"),
    ("06/figure4.png", "figure4_solar_wind.png"),
    ("08/profile_fit.png", "posterior_profile_fit.png"),
    ("08/corner.png", "posterior_corner.png"),
)


def headas_command(script):
    settings = json.loads(
        (REPRO / "analysis_01_reduction" / "settings.json").read_text()
    )
    headas = os.environ.get("HEADAS", settings["headas"])
    atomdb = os.environ.get("ATOMDB", settings["atomdb"])
    lib = Path(sys.executable).resolve().parents[1] / "lib"
    shell = (
        f'export HEADAS="{headas}" && source "$HEADAS/headas-init.sh" && '
        f'export ATOMDB="{atomdb}" && '
        f'export LD_LIBRARY_PATH="{lib}:${{LD_LIBRARY_PATH:-}}" && '
        f'exec "{sys.executable}" "{script}"'
    )
    return ["bash", "-c", shell]


def load(number):
    return json.loads((OUT / number / "result.json").read_text())


def summarize():
    results = {number: load(number) for number, _ in STAGES}
    r3, r5, r7, r8 = results["03"], results["05"], results["07"], results["08"]
    summary = {
        "pps_luminosity": {
            "epoch_utc": results["02"]["epoch_utc"],
            "full_aperture_luminosity_erg_s": results["02"]["flux_luminosity"][
                "full_aperture"
            ]["luminosity_erg_s"],
        },
        "spectrum_pn_only": {
            "chi2": r3["fit"]["chi2"],
            "dof": r3["fit"]["dof"],
            "cx_temperature_keV": r3["fit"]["parameters"]["cx_temperature_keV"],
            "cx_collision_keV_u": r3["fit"]["parameters"]["cx_collision_keV_u"],
            "mean_photon_energy_eV": r3["mean_photon_energy_eV"],
            "inner_luminosity_erg_s": r3["regions"]["inner"]["luminosity_erg_s"],
            "inner_luminosity_1sigma_erg_s": r3["regions"]["inner"][
                "luminosity_1sigma_erg_s"
            ],
        },
        "morphology": {
            "lowest_contour_extent_arcmin": results["04"][
                "lowest_contour_extent_arcmin"
            ],
        },
        "radial_profile": {
            "slope": r5["epic"]["slope"],
            "sector_standard_error": r5["epic"]["sector_standard_error"],
            "broken": {
                key: r5["epic"]["broken"][key]
                for key in ("alpha_inner", "alpha_outer", "break_radius_arcmin")
            },
            "luminosity_370_arcsec_erg_s": r5["apertures"]["370_arcsec"][
                "luminosity_erg_s"
            ],
        },
        "solar_wind": {
            "proton_flux_medians_cm-2_s-1": {
                name: results["06"]["proton_flux"][name]["median"]
                for name in ("stereo_a", "solar_orbiter", "wind", "ace_scaled")
                if name in results["06"]["proton_flux"]
            },
            "f107p_sfu": results["06"]["f107"]["f107p_sfu"],
        },
        "scaling": {
            "xmm_parent_rate_sum_s-1": r7["parent_rates"]["xmm_sum_s-1"],
            "cabot_370_missing_fraction": r7["cabot_370_slow_speed_check"][
                "missing_fraction"
            ],
            "n2_alternative_Q_N2_over_Q_CO": r7["n2_alternative"]["Q_N2_over_Q_CO"],
        },
        "chemistry_posterior": {
            "mean_acceptance_fraction": r8["mean_acceptance_fraction"],
            "medians": {
                name: r8["posterior"][name]["median"] for name in r8["posterior"]
            },
            "H2_fraction_interval": (
                [
                    r8["posterior"]["H2_fraction"]["q16"],
                    r8["posterior"]["H2_fraction"]["q84"],
                ]
                if "H2_fraction" in r8["posterior"]
                else None
            ),
        },
        "recoil": {
            "radial_recoil_headline": results["09"]["radial_recoil"]["headline"],
            "oumuamua_H2_requirement_molecules_s": results["09"]["oumuamua"][
                "H2_requirement_molecules_s"
            ],
            "radiolytic_gas_erosion": results["09"]["radiolytic_capacity"][
                "gas_erosion"
            ],
        },
    }
    (OUT / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    figures = OUT / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    for source, name in FIGURES:
        shutil.copyfile(OUT / source, figures / name)


def main():
    start = time.time()
    for number, label in STAGES:
        script = DEMO / "stages" / f"stage{number}.py"
        command = (
            headas_command(script)
            if number == "03"
            else [sys.executable, str(script)]
        )
        print(f"[{time.time() - start:7.1f}s] stage {number}: {label}", flush=True)
        tick = time.time()
        completed = subprocess.run(command, cwd=DEMO / "stages")
        if completed.returncode:
            raise SystemExit(
                f"stage {number} failed with exit code {completed.returncode}"
            )
        print(
            f"[{time.time() - start:7.1f}s] stage {number} finished"
            f" ({time.time() - tick:.1f}s)",
            flush=True,
        )
    summarize()
    print(
        f"[{time.time() - start:7.1f}s] complete:"
        " outputs/results.json, outputs/figures/"
    )


if __name__ == "__main__":
    main()
