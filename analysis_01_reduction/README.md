# Analysis 01: XMM reduction

This task reduces ObsID 0963720201 into the masked EPIC event, background,
exposure, and comet-frame products used by the paper analyses.

Run under WSL with the recorded Python environment:

```bash
export HEADAS=/home/scabot/Programs/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39
source "$HEADAS/headas-init.sh"
source /home/scabot/Programs/sas/xmmsas_22.1.0-a8f2c2afa-20250304/setsas.sh
export LD_LIBRARY_PATH=/home/scabot/miniconda3/envs/xmm/lib:${LD_LIBRARY_PATH:-}
/home/scabot/miniconda3/envs/xmm/bin/python run.py
```

The original ODF and PPS trees are read-only inputs.

PN background products use an OoT image rebuilt from the OoT event list. Hard
soft-proton maps use the signed all-filter FLARE B4--B6 calibration and fitted
`proton` total.
