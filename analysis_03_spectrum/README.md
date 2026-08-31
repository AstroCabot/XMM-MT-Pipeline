# Analysis 03: EPIC spectrum

This task extracts the masked 0–370 and 400–600 arcsec spectra for PN, MOS1,
and MOS2, fits the joint sky plus VACX2 model, and writes Figure 2.

Inputs are the Analysis 01 masked events, ephemeris, GTIs, background products,
PN OoT events, and frame responses.

Extract in a sourced SAS/HEASoft environment:

```bash
export LD_LIBRARY_PATH=/home/scabot/miniconda3/envs/xmm/lib:${LD_LIBRARY_PATH:-}
/home/scabot/miniconda3/envs/xmm/bin/python run.py extract
```

Fit in a fresh HEASoft-only shell with `ATOMDB` set:

```bash
/home/scabot/miniconda3/envs/xmm/bin/python run.py fit
```

The SP spectrum uses the Analysis 01 fitted normalization and supplied zero
`STAT_ERR`.
`live_solid_angle_arcmin2` is the source-masked geometric sky area before
detector coverage.
MOS spectra and responses use Task 01's soft-band anomalous-CCD exclusions.
`detector_K_ct_s_per_erg_cm2_s` is the inner on-axis conversion, including the
fitted detector constant. The regional ARF conversion excludes throughput
already represented in the image exposure planes.
The fit uses calibrated RMF energy bounds; each image conversion folds the
fitted CX model through the grouped PHA bins containing PI=350--1100 eV.
This is raw PHA channels 70--220 inclusive; the final 5 eV bin is kept whole.
`image_response_ratio_to_inner` records the outer detector-mixture correction in
inner-MOS2 image units. The ECF log uncertainty combines the MOS2 profile
interval with the XMM-SOC-CAL-TN-0018 absolute-area term `ln(1.10)`.
Outputs are `extraction.tsv`, `spectra.tsv`, `stack/`, `folded.tsv`,
`result.json`, and `figure2.png`/`figure2.pdf`.
