# XMM-Newton comet moving-target pipeline

This package is a cleaned refactor of an XMM-Newton EPIC reduction workflow for a moving, extended comet observed across multiple pointings. It is aimed at comet-centric, publishable imaging products together with conservative timing and spectral products, and it writes a large QC bundle so each stage can be inspected while the reduction is running.

## What the refactor does

The pipeline is organized around a few choices that are better suited to a diffuse moving target.

- Still-sky source detection stays SAS-based and mosaic-aware.
- Detection products are split into a full field-source list and a track-excluded comparison list.
- Contamination handling is policy-based instead of being hard-wired to one aggressive veto.
- A stationary-source trail mask is built explicitly and can be applied directly to comet-frame event lists.
- Imaging and QC are soft-band aware by default.
- QC products are produced at every major stage, including per-detector image panels and per-segment spot checks.

## Main products

The reducer writes products under `WORKDIR/`.

### Science products

- `clean/` — cleaned per-exposure files and merged per-instrument event lists
- `track/` — comet ephemeris FITS and motion-segmentation CSV
- `detect/` — detect products, pseudo-exposure manifest, still-sky mosaics, source lists
- `comet/` — comet-frame per-instrument event lists
- `contam/` — policy GTIs, contamination reports, trail mask, masked science-base event lists
- `images/<band>/` — counts, exposure, rate, and masked-rate images
- `lcurve/` — per-instrument and combined light curves
- `spectra/` — per-instrument spectra/responses and combined EPIC products
- `final/` — merged EPIC event lists for the full comet frame and the selected science GTI

### QC products

- `qc/00_manifest.*` — stage manifest and package health overview
- `qc/00_clean.png` — cleaning summary by instrument/exposure
- `qc/01_track.png` — track and motion-segmentation check
- `qc/02_detect.png` — soft/hard still-sky detect diagnostics with all-vs-curated sources
- `qc/03_image.png` plus `qc/03_image_<band>.png` — comet-frame image diagnostics with per-detector panels
- `qc/04_contam.png` — contamination-policy timeline and most problematic sources
- `qc/05_lcurve.png` — per-instrument and combined light curves
- `qc/06_spectrum.png` — spectrum quick look
- `qc/07_spotchecks_<band>.png` — per-segment spot checks in the comet frame
- `qc/region_support_*.json` — aperture-support diagnostics
- `qc/STAGE_QC.md` — stage-by-stage QC navigation page

## Recommended defaults for a diffuse soft comet

The included configuration template is tuned for a diffuse, soft source.

- `IMAGE_BANDS=soft:200:1000;broad:300:2000;hard:1000:12000`
- `DETECT_QC_BANDS=soft:200:1000;hard:1000:12000`
- `SCIENCE_GTI_POLICY=full`
- `SCIENCE_APPLY_TRAIL_MASK=yes`
- `LC_APPLY_ABSOLUTE_CORRECTIONS=no`

That default combination keeps exposure unless there is a strong reason to cut it, masks stationary-source trails spatially, and emphasizes soft-band imaging diagnostics.

## Quick start

Before running the helper scripts, install the Python requirements in `requirements/python.txt`.

1. Copy the template configuration.

```bash
cp config/xmm_comet_config.template.env my_comet.env
```

2. Edit at least:

- `WORKDIR`
- `ODFDIR`
- either `TARGET_ID` or `TRACK_INPUT`

3. Run the full reducer.

```bash
./xmm_comet_pipeline.sh my_comet.env all
```

4. Or rerun just the science stages once `init/` and `repro/` already exist.

```bash
FORCE=1 ./xmm_comet_pipeline.sh my_comet.env science
```

5. Build the QC bundle.

```bash
./xmm_comet_run_qc.sh my_comet.env
```

6. Start inspection at:

```text
$WORKDIR/qc/STAGE_QC.md
```

## Configuration highlights

### Detection and source lists

The detect stage produces both:

- `detect/field_sources_all.fits`
- `detect/field_sources_curated.fits`

The full list is the default input for contamination accounting and trail masking. The track-excluded list is retained as a comparison product for QC and manual review.

### Contamination policies

The contamination builder writes four GTI policies:

- `full` — keep all sampled times
- `src` — reject only source-aperture overlaps
- `bkg` — reject only background-aperture overlaps
- `strict` — reject either overlap type

The selected policy is copied to `contam/science_gti.fits`.

### Spatial trail mask

When `SCIENCE_APPLY_TRAIL_MASK=yes`, the pipeline writes:

- `contam/EPIC_trail_mask.fits`
- `contam/<INST>_science_base.fits`

Those masked event lists become the default inputs for light-curve, spectrum, and science-merge stages.

## Package layout

- `xmm_comet_pipeline.sh` — main pipeline entry point
- `xmm_comet_run_qc.sh` — QC wrapper
- `config/` — configuration template
- `docs/` — concise operational notes
- `requirements/` — Python helper requirements
- `scripts/` — helper scripts
- `tests/` — synthetic tests for the refactor helpers
- `PACKAGE_MANIFEST.json` — package manifest

## Running the synthetic helper tests

The included tests cover the refactor logic that is easiest to validate without SAS:

- empty-source-list handling
- policy-based contamination products
- trail-mask generation and event filtering

Run them from the package root with:

```bash
python -m unittest discover -s tests -v
```

## Notes on QC interpretation

A short guide to the main QC plots is in `docs/QC_GUIDE.md`.
