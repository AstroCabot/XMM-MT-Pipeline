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

The run builds the PSF/operator and samples 96 walkers for 500 burn-in and 1,500
retained steps. Outputs are `result.json`, `chemistry.tsv`,
`posterior.tsv`, `profile.tsv`, `intrinsic.tsv`, `chain.npz`, `psf.npz`,
`figure3.png`/`figure3.pdf`, and `corner.png`.

`log10_H2_ratio >= -0.5` defines the high mode and `<= -1.5` the low mode.
`parameters.json` defines the priors; `transport.py` and `inference.py` implement
the physical model and likelihood.
