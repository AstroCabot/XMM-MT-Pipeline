# Refactor notes

This package is the cleaned refactor of the comet pipeline.

## Main implementation changes

### 1. The detect stage now keeps two source lists

- `detect/field_sources_all.fits`
- `detect/field_sources_curated.fits`

The full list is the default input for contamination accounting and trail masking. The track-excluded list is retained as a comparison product for QC and manual inspection.

### 2. Contamination handling is policy-based

`build_contamination_products.py` writes:

- `gti_full.fits`
- `gti_src.fits`
- `gti_bkg.fits`
- `gti_strict.fits`
- `science_gti.fits`
- `contamination_timeline.csv`
- `contamination_report.csv`
- `contamination_summary.json`

The selected policy is controlled by `SCIENCE_GTI_POLICY`.

### 3. Spatial masking is applied to event lists

`scripts/trail_mask_tools.py` provides two subcommands:

- `trail_mask_tools.py make`  — build the 2-D stationary-source trail mask
- `trail_mask_tools.py apply` — filter events that fall inside the mask

This lets later science stages use masked event lists instead of relying only on time-domain vetoes.

### 4. Imaging and QC defaults are soft-band aware

The template defaults to:

- `IMAGE_BANDS=soft:200:1000;broad:300:2000;hard:1000:12000`
- `DETECT_QC_BANDS=soft:200:1000;hard:1000:12000`

### 5. Light-curve defaults are conservative for a diffuse source

The package defaults to relative-only correction unless the configuration explicitly requests absolute correction.

### 6. Empty source lists are handled explicitly

The helper scripts now treat empty detect/source tables as valid inputs and write empty-but-well-formed downstream products instead of failing late.

### 7. Spot-check masking is segment-limited

The spot-check QC step clips the trail mask to the segment time range so each panel reflects the selected frame interval.

### 8. QC orchestration consolidated into Python

The old `xmm_comet_spotcheck_frames.sh` and `xmm_qc_build_mosaics.sh` shell helpers were replaced by a single `scripts/xmm_qc_runner.py` driver, invoked via the thin `xmm_comet_run_qc.sh` wrapper.

### 9. Comet-frame attitude built as a standalone tool

The inline Python that was embedded in the pipeline shell script for building the moved ATTHK is now `scripts/build_moved_atthk.py`. This makes the comet-tracking attitude correction testable and inspectable outside the full pipeline.
