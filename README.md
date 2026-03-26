# XMM-Newton comet moving-target pipeline (lean package)

This is a further-compacted repack of the refactored 3I/ATLAS EPIC pipeline.
It keeps the same main entry points and stage outputs, but trims the QC wrapper layer and folds the trail-mask utilities into one module.

## Entry points

- `xmm_comet_pipeline.sh`
- `xmm_comet_run_qc.sh`

Use the included starting config:

```bash
./xmm_comet_pipeline.sh my_comet.lean.env all
./xmm_comet_run_qc.sh my_comet.lean.env
```

## Retained behavior

- staged reduction and QC layout
- moved-target event lists and moved-ATTHK handling
- per-instrument and EPIC-combined images, light curves, spectra, and merged event lists
- comet-frame trail masking and image post-processing
- still-sky detect mosaics and comet-frame spot-check products

## Lean-package changes

- consolidated the QC shell helpers into `scripts/xmm_qc_runner.py`
- reduced `xmm_comet_run_qc.sh` to a thin wrapper
- merged `trail_mask_utils.py` into `trail_mask_tools.py`
- kept the earlier science-facing fixes and stage interfaces intact
