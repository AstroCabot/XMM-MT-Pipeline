import csv
import json
from pathlib import Path

import prior_table

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "outputs"
RESULT_ROWS = (
    ("phi_heavy", r"\(\Phi_{\rm heavy}\) (cm\(^{-2}\)\,s\(^{-1}\))"),
    ("log10_H2_ratio", r"\(\log_{10}r_{\rm H_2}\)"),
    ("Q_H2O", r"\(Q_{\rm H_2O}\) (s\(^{-1}\))"),
    ("Q_CO2", r"\(Q_{\rm CO_2}\) (s\(^{-1}\))"),
    ("Q_CO", r"\(Q_{\rm CO}\) (s\(^{-1}\))"),
    ("Q_H2", r"\(Q_{\rm H_2}\) (s\(^{-1}\))"),
    ("H2_fraction", r"\(f_{\rm H_2}\)"),
    ("Q_total", r"\(Q_{\rm total}\) (s\(^{-1}\))"),
    ("beta_H2O", r"\(\beta_{\rm H_2O}\)"),
    ("beta_CO2", r"\(\beta_{\rm CO_2}\)"),
    ("beta_CO", r"\(\beta_{\rm CO}\)"),
    ("beta_H2", r"\(\beta_{\rm H_2}\)"),
    ("scale", r"profile scale $a$"),
    ("offset", r"profile offset $b$ (count\,s\(^{-1}\)\,arcmin\(^{-2}\))"),
    ("jitter", r"jitter $j$"),
)


def table(rows, caption, label, columns, header, wide=False):
    environment = "table*" if wide else "table"
    return "\n".join((
        rf"\begin{{{environment}}}[!ht]",
        r"\centering", r"\singlespacing", r"\small",
        r"\renewcommand{\arraystretch}{0.9}" if wide else "",
        rf"\setlength{{\tabcolsep}}{{{5 if wide else 8}pt}}",
        rf"\caption{{{caption}}}",
        rf"\begin{{tabular}}{{{columns}}}",
        r"\toprule", header + r"\\", r"\midrule",
        *rows,
        r"\bottomrule", r"\end{tabular}", rf"\label{{{label}}}",
        rf"\end{{{environment}}}", "",
    ))


def posterior_rows(path):
    with path.open() as stream:
        values = {row["parameter"]: row for row in csv.DictReader(stream, delimiter="\t")}
    rows = []
    for name, label in RESULT_ROWS:
        if name == "scale":
            rows.append(r"\addlinespace[3pt]")
        cells = []
        for key in ("map", "q16", "median", "q84"):
            value = float(values[name][key])
            number = prior_table.science(value, 2) if name.startswith("Q_") or name in ("phi_heavy", "offset") else f"{value:.2f}"
            cells.append(rf"\({number}\)")
        rows.append(" & ".join((label, *cells)) + r"\\")
    return rows


def main():
    wind = json.loads(prior_table.WIND.read_text())
    spectrum = json.loads(prior_table.SPECTRUM.read_text())
    profile = json.loads(prior_table.PROFILE.read_text())
    background = profile["epic"]["background_prior"]["sigma_ct_s_arcmin2"]
    f107p = float(wind["f107"]["f107p_sfu"])
    captions = (
        r"Parameter specification for the chemical transport model (Part I). Notation: \(\mathcal U\), uniform; \(\mathcal N(\mu,\sigma)\), normal; \(\mathcal HN(\sigma)\), half-normal; and \(\mathcal{LN}(m,\sigma_{\ln})\), lognormal with median \(m\) and logarithmic width \(\sigma_{\ln}\). Normal distributions on strictly positive physical quantities are truncated at zero.",
        r"Parameter specification for the chemical transport model (Part II). This table contains values which were frozen in the model specification, mainly pertaining to photodestruction and photoionization. Rates are calculated by inputting "
        + rf"$F_{{\rm 10.7P}}={f107p:.1f}$ sfu "
        + r"into the power-law fits of Moore and Mapaye\cite{MooreMapaye2025}.",
        r"Sampled parameter distribution for the full-profile chemical transport model. The columns give the MAP and central percentiles of the full retained chain. Production rates and $f_{\rm H_2}$ refer to the \xmm\ epoch.",
    )
    tables = (
        table(prior_table.model_rows(wind, spectrum, background), captions[0], "tab:model",
              r"@{}p{0.30\textwidth}p{0.30\textwidth}p{0.34\textwidth}@{}",
              "Quantity & Prior & Source", wide=True),
        table(prior_table.rate_rows(f107p), captions[1], "tab:photorates", "@{}llr@{}",
              "Quantity & Value & Source"),
        table(posterior_rows(OUTPUT / "posterior.tsv"), captions[2], "tab:results",
              "@{}lcccc@{}", "Quantity & MAP & 16th & 50th & 84th"),
    )
    for number, content in enumerate(tables, 3):
        (OUTPUT / f"table{number}.tex").write_text(content)


if __name__ == "__main__":
    main()
