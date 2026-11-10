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
Set their paths and the software paths in `settings.json` before running.
With `summary: null`, initialization builds the calibration index and uses
SAS `odfingest` to create the ODF summary in the copied ODF tree; an existing
summary is reused. An externally prepared summary can instead be supplied
as an absolute path in `summary`. This follows the
[SAS setup sequence](https://xmm-tools.cosmos.esa.int/external/xmm_user_support/documentation/sas_usg/USG/odfingest.html).

PN background products use an OoT image rebuilt from the OoT event list. Hard
soft-proton maps use the signed all-filter FLARE B4--B6 calibration and fitted
`proton` total.
