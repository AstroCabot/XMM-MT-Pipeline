import json
from pathlib import Path

import chemistry
from config import CONFIG, SHARED

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "outputs"
SPECTRUM = ROOT / "analysis_03_spectrum/outputs/result.json"
PROFILE = ROOT / "analysis_05_radial_profile/outputs/result.json"
WIND = ROOT / "analysis_06_solar_wind/outputs/result.json"
SPHEREX = r"SPHEREx\cite{Lisse2026}"
FRINCKE = r"Frincke \& Seligman\cite{Frincke2026}"
MOORE = r"Moore \& Mapaye\cite{MooreMapaye2025}"
MACHACEK = r"Machacek et al.\cite{Machacek2015}"
COMBI = r"Combi et al.\cite{Combi2004}"
LABELS = {
    "H2O": r"\HtwoO",
    "CO2": r"\COtwo",
    "H2": r"\Htwo",
    "O2": r"\Otwo",
    "OH": "OH",
    "CO": "CO",
    "O": "O",
    "C": "C",
    "H": "H",
}
LIFETIMES = ("H2O", "OH", "CO2", "CO", "H2", "O", "C", "H", "O2")
YIELDS = (
    ("H2O", ("OH", "H", "O", "H2")),
    ("OH", ("O", "H")),
    ("CO2", ("CO", "O", "O2")),
    ("CO", ("C", "O")),
    ("H2", ("H",)),
    ("O2", ("O",)),
)


def science(value, digits):
    mantissa, exponent = f"{value:.{digits}e}".split("e")
    return rf"{mantissa}\times10^{{{int(exponent)}}}"


def share(rows, parent, selected):
    rates = [row for row in rows if row["parent"] == parent]
    return sum(row["rate_s-1"] for row in rates if selected(row)) / sum(
        row["rate_s-1"] for row in rates
    )


def photochemistry(f107p):
    rows = chemistry.reactions(HERE / "reactions.tsv", f107p)
    network = chemistry.network(rows, CONFIG["cross_section_cm2"])
    oxygen = [row for row in rows if row["parent"] == "H2O" and row["O"] > 0]
    molecular = sum(row["rate_s-1"] for row in oxygen if row["H2"] > 0)
    return {
        "lifetime_s": {name: 1 / network[name]["rate_1au_s-1"] for name in LIFETIMES},
        "branching": {
            "CO_to_C_O": share(rows, "CO", lambda r: r["C"] > 0 and r["O"] > 0),
            "CO2_to_CO_O": share(rows, "CO2", lambda r: r["CO"] > 0 and r["O"] > 0),
            "H2O_oxygen": share(rows, "H2O", lambda r: r["O"] > 0),
            "H2O_H2_split": molecular / sum(row["rate_s-1"] for row in oxygen),
            "H2O_ionization": share(
                rows, "H2O", lambda r: r["products"].rstrip().endswith("+ e")
            ),
            "OH_to_O_H": share(rows, "OH", lambda r: r["O"] > 0 and r["H"] > 0),
            "H2_to_H": dict(network["H2"]["children"])["H"],
        },
        "yields": {parent: dict(network[parent]["children"]) for parent, _ in YIELDS},
    }


def section(title):
    return [
        "\\addlinespace[4pt]",
        rf"\multicolumn{{3}}{{@{{}}l}}{{\itshape {title}}}\\",
        "\\addlinespace[2pt]",
    ]


def row(quantity, value, source):
    return [rf"{quantity} & {value} & {source}\\"]


def yield_text(value):
    return f"{value:.3g}" if value >= 1 else f"{value:.2g}"


def model_rows(wind, spectrum, background):
    priors = SHARED["spherex"]["parent_priors"]
    p = CONFIG["priors"]
    sigma = CONFIG["cross_section_cm2"]
    speed = CONFIG["fixed_velocity_km_s"]
    flux = wind["heavy_ion_prior"]
    steep = priors["H2O"]["beta_steep"]
    lines = section(
        f"Neutral production at $r_h={SHARED['spherex']['reference_distance_au']}$~au"
    )[1:]
    for name, tag in (("H2O", "H_2O"), ("CO2", "CO_2"), ("CO", "CO")):
        rate = priors[name]["reference_rate_s-1"]
        lines += row(
            rf"\(Q_{{\rm ref,{tag}}}\) (s$^{{-1}}$)",
            rf"\(\mathcal{{LN}}({science(rate['median'], 2)},"
            rf"{rate['sigma_ln']:.2f})\)",
            SPHEREX,
        )
    low, high = p["log10_H2_ratio"]
    lines += row(
        r"\(\log_{10}r_{\rm H_2}\)",
        rf"\(\mathcal U({low:g},{high:g})\)",
        "ratio to observed parent species",
    )
    lines += section("Post-perihelion production histories")
    lines += row(
        r"\(\beta_{\rm H_2O}\)",
        rf"\(\mathcal{{N}}({p['beta_H2O'][0]},{p['beta_H2O'][1]})\)",
        r"Tan et al.\cite{Tan2026}",
    )
    lines += row(
        r"\(\beta_{\rm H_2O}\), comparison fit",
        rf"\(\mathcal{{N}}({steep['mean']:.3f},{steep['sigma']:.3f})\)",
        FRINCKE,
    )
    for name, tag in (("CO2", "CO_2"), ("CO", "CO")):
        mean, width = p[f"beta_{name}"]
        lines += row(
            rf"\(\beta_{{\rm {tag}}}\)",
            rf"\(\mathcal{{N}}({mean:.3f},{width:.3f})\)",
            FRINCKE,
        )
    low, high = p["beta_H2"]
    lines += row(
        r"\(\beta_{\rm H_2}\)",
        rf"\(\mathcal{{U}}({low:g},{high:g})\)",
        "shallow evolution profile",
    )
    lines += section("Solar wind and X-ray emission")
    lines += row(
        r"\(\Phi_{\rm heavy}\) (cm\(^{-2}\)\,s\(^{-1}\))",
        rf"\(\mathcal{{U}}({science(flux['lower'], 2)},"
        rf"{science(flux['upper'], 2)})\)",
        r"mapped proton-flux envelope \(\times f\)",
    )
    lines += row(
        r"\(\Delta E\) (eV)",
        f"{spectrum['mean_photon_energy_eV']:.1f}, fixed",
        "count-weighted band mean",
    )
    lines += row(
        r"ECF\(_{\rm MOS2}\) (erg\,cm\(^{-2}\)\,cts\(^{-1}\))",
        rf"\({science(spectrum['mos2_ecf_erg_cm2_per_count'], 4)}\), fixed",
        "MOS2 energy conversion factor",
    )
    lines += section("Radial-profile nuisance parameters")
    lines += row(
        r"\(a\)",
        rf"\(\mathcal{{LN}}(1,{spectrum['mos2_ecf_sigma_ln']:.2f})\)",
        r"MOS2 ECF uncertainty\cite{XMMCal2025}",
    )
    lines += row(
        r"\(b\) (count\,s\(^{-1}\)\,arcmin\(^{-2}\))",
        rf"\(\mathcal{{N}}(0,{science(background, 1)})\)",
        "residual flat-sky term",
    )
    lines += row(
        r"\(j\)",
        rf"\(\mathcal{{HN}}({p['jitter_sigma']:.2f})\)",
        "fractional profile scatter",
    )
    lines += section(r"Charge-exchange cross sections (cm\(^2\))")
    for name, tag, source in (
        ("H2O", "H_2O", MACHACEK),
        ("CO", "CO", MACHACEK),
        ("CO2", "CO_2", MACHACEK),
        ("H2", "H_2", r"Machacek et al.\cite{Machacek2014}"),
    ):
        lines += row(
            rf"\(\sigma_{{\rm {tag}}}\)",
            rf"\({science(sigma[name], 2)}\), fixed",
            source,
        )
    lines += row(
        r"\(\sigma_{\rm OH,O,C,O_2,H}\)",
        rf"\({science(sigma['daughter_proxy'], 2)}\), fixed",
        r"molecule-average as proxy\cite{Machacek2015}",
    )
    lines += section(r"Outflow and recoil speeds (km\,s$^{-1}$)")
    mean, width = p["parent_velocity"]
    lines += row(
        r"\HtwoO; \COtwo; CO; \Htwo, nucleus",
        rf"\(\mathcal{{N}}({mean:.2f},{width:.2f})\)",
        r"IRAM\cite{Biver2026}",
    )
    low, high = p["OH_velocity"]
    lines += row(
        "OH",
        rf"\(\mathcal{{U}}({low:.2f},{high:.2f})\)",
        rf"IRAM\cite{{Biver2026}}, {COMBI}",
    )
    low, high = p["O_velocity"]
    lines += row(
        r"O from \HtwoO or OH", rf"\(\mathcal{{U}}({low:.2f},{high:.2f})\)", COMBI
    )
    low, high = p["H_velocity"]
    lines += row("H", rf"\(\mathcal{{U}}({low:.1f},{high:.1f})\)", COMBI)
    lines += row(
        r"CO; O; \Otwo; C, daughters",
        f"{speed['heavy_daughter']:.1f}, fixed",
        r"midpoint of range\cite{Combi2004}",
    )
    lines += row(
        r"Photolytic \Htwo",
        f"{speed['photolytic_H2']:.1f}, fixed",
        r"midpoint of range\cite{Combi2004}",
    )
    return lines


def rate_rows(f107p):
    values = photochemistry(f107p)
    branching = values["branching"]
    lines = section(
        rf"Photodestruction lifetimes at 1 au and \(F_{{\rm 10.7P}}={f107p:.1f}\)~sfu"
        r" (s), scaled as \(r_h^2\)"
    )[1:]
    for name in LIFETIMES:
        lines += row(
            LABELS[name], rf"\({science(values['lifetime_s'][name], 3)}\)", MOORE
        )
    lines += section("Branching yields")
    for label, name in (
        (r"CO \(\rightarrow\) C + O", "CO_to_C_O"),
        (r"\COtwo\ \(\rightarrow\) CO + O", "CO2_to_CO_O"),
        (r"\HtwoO, O-producing", "H2O_oxygen"),
        (r"\HtwoO, \(\mathrm{H_2}\!:\!\mathrm{2H}\) split", "H2O_H2_split"),
        (r"\HtwoO, ionization", "H2O_ionization"),
        (r"OH \(\rightarrow\) O + H", "OH_to_O_H"),
    ):
        lines += row(label, f"{branching[name]:.4f}", MOORE)
    lines += row(
        r"\Htwo\ \(\rightarrow\) H multiplicity",
        f"{branching['H2_to_H']:.3f}",
        MOORE,
    )
    lines += section(
        rf"Mean neutral yields per destruction (\(F_{{\rm 10.7P}}={f107p:.1f}\)~sfu)"
    )
    for parent, children in YIELDS:
        amounts = values["yields"][parent]
        if parent == "OH":
            text = f"{yield_text(amounts['O'])} (O + H)"
        else:
            text = " + ".join(
                f"{yield_text(amounts[child])} {LABELS[child]}" for child in children
            )
        lines += row(LABELS[parent], rf"\(\rightarrow\) {text}", MOORE)
    return lines


def main():
    wind = json.loads(WIND.read_text())
    spectrum = json.loads(SPECTRUM.read_text())
    profile = json.loads(PROFILE.read_text())
    background = profile["epic"]["background_prior"]["sigma_ct_s_arcmin2"]
    f107p = float(wind["f107"]["f107p_sfu"])
    lines = [
        "% Extended Data Table 3",
        *model_rows(wind, spectrum, background),
        "",
        f"% Extended Data Table 4, F10.7P = {f107p:.1f} sfu",
        *rate_rows(f107p),
    ]
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "prior_table.tex").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
