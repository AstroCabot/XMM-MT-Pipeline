# Analysis 08: chemical profile

This task propagates the four released parents and their neutral photoproducts,
applies the exposure-weighted EPIC PSF, and fits the ten Task 05 annuli in count
rate surface brightness from their raw count and exposure sufficient statistics.

Inputs are the Analysis 01 comet-frame products, Analysis 03 response
calibration, Analysis 04 exposure plane, Analysis 05 annuli and background, and
Analysis 06 solar inputs. Run under WSL:

```bash
export HEADAS=/home/scabot/Programs/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39
source "$HEADAS/headas-init.sh"
export SAS_CCFPATH=/home/scabot/Programs/sas/ccf
export SAS_CCF=../analysis_01_reduction/outputs/01_init/ccf.cif
export SAS_ODF=../analysis_01_reduction/outputs/01_init/odf/4759_0963720201_SCX00000SUM.SAS
source /home/scabot/Programs/sas/xmmsas_22.1.0-a8f2c2afa-20250304/setsas.sh
export LD_LIBRARY_PATH=/home/scabot/miniconda3/envs/xmm/lib:${LD_LIBRARY_PATH:-}
/home/scabot/miniconda3/envs/xmm/bin/python run.py
```

The run builds the PSF/operator and samples 96 walkers for 5,000 burn-in and
20,000 retained steps with emcee differential-evolution moves and four workers.
Set BLAS/OpenMP thread counts to one when using multiple workers.
`inference.Gaussianized` maps each parameter through its prior CDF to a unit
normal for sampling; chains and log probabilities are reported in the physical
parameters. Outputs are `result.json`, `chemistry.tsv`, `posterior.tsv`,
`profile.tsv`, `intrinsic.tsv`, `chain.npz`, `psf.npz`,
`figure3.png`/`figure3.pdf`, and `corner.png`.

The baseline water history uses `beta_H2O ~ N(3.3, 0.3)` (Tan et al. 2026). The
steeper Frincke & Seligman (2026, arXiv:2606.01377v3) post-perihelion slope is a
second fit on the same PSF:

```bash
/home/scabot/miniconda3/envs/xmm/bin/python run.py fit --water-prior frincke
```

It samples `beta_H2O ~ N(8.253, 0.560)` and writes its products, except
`psf.npz`, to `outputs/water_frincke/`. Task 09 reads them through its
`--task08`, `--chain`, and `--output` options. Production stays at the
perihelion rate before perihelion in both fits. `result.json` records the
`beta_H2O` prior used.

The conventional model fixes directly released H2 to zero, retaining all
photochemical daughters, including H2 produced from water. It samples 17
parameters, omitting the direct-H2 abundance, history and release speed:

```bash
/home/scabot/miniconda3/envs/xmm/bin/python run.py fit --no-direct-h2
/home/scabot/miniconda3/envs/xmm/bin/python run.py fit --no-direct-h2 --water-prior frincke
```

These fits share `outputs/psf.npz` and write to `outputs/no_direct_h2/` and
`outputs/no_direct_h2/water_frincke/`. Use the same flags with `corner`.
Their profile uses the overall MAP; H2 mode comparisons are omitted.

`fit` saves `outputs/chain.pkl`, in those unit-normal coordinates, every 500
steps. `run.py fit --resume` continues it exactly, and `--retained-steps N`
extends a finished run to N retained steps.
`fit` refuses to start while `chain.pkl` exists without `--resume`; remove the
file when any input changes, because resuming does not check them.

`python prior_table.py` writes the Extended Data Table 3 and 4 rows to
`outputs/prior_table.tex` from the inputs the fit reads. It needs no fit; rerun
it after Tasks 03, 05, or 06 change.

`python paper_tables.py` exports complete Extended Data Tables 3--5 to
`outputs/table3.tex`, `table4.tex`, and `table5.tex`, using the same current
inputs and `outputs/posterior.tsv`. The tables use the paper's labels and
macros and need `booktabs` and `setspace`.

`log10_H2_ratio >= -0.5` defines the high mode and `<= -1.5` the low mode.
`parameters.json` defines the priors; `transport.py` and `inference.py` implement
the physical model and likelihood.

`transport.py` evaluates daughter and granddaughter fluxes on nodes at every
second shell edge (`transport_node_stride`) and interpolates them to the 320
shell centers. The birth integrals use photochemical exposure, the integral of
`r_h^-2` along a flight, as the variable, so constant-rate daughter profiles are
exact. Against direct quadrature at 75 parameter points the annular predictions
agree within 6e-4 (median 3e-5).
