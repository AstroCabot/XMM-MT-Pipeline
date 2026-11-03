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
20,000 retained steps. It uses emcee differential-evolution moves in coordinates where
every prior is a unit normal: `inference.Gaussianized` maps each parameter
through its prior CDF, so the likelihood, priors, and posterior are unchanged,
and chains and log probabilities are reported in the physical parameters.
Outputs are `result.json`, `chemistry.tsv`,
`posterior.tsv`, `profile.tsv`, `intrinsic.tsv`, `chain.npz`, `psf.npz`,
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

`fit` saves `outputs/chain.pkl`, in those unit-normal coordinates, every 500
steps. `run.py fit --resume` continues it exactly, and `--retained-steps N`
extends a finished run to N retained steps.
`fit` refuses to start while `chain.pkl` exists without `--resume`; remove the
file when any input changes, because resuming does not check them.

`python prior_table.py` writes the Extended Data Table 3 and 4 rows to
`outputs/prior_table.tex` from the inputs the fit reads. It needs no fit; rerun
it after Tasks 03, 05, or 06 change.

`log10_H2_ratio >= -0.5` defines the high mode and `<= -1.5` the low mode.
`parameters.json` defines the priors; `transport.py` and `inference.py` implement
the physical model and likelihood.
