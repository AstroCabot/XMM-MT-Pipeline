# QC guide

This is a short interpretation guide for the main QC products.

## `02_detect.png`

This plot shows still-sky diagnostic mosaics, typically in both a soft band and a hard band.

- The soft panel is the more comet-relevant diagnostic if the coma is strongest around 0.2–1.0 keV.
- The hard panel is mainly useful for stationary field-source identification and background contamination checks.
- White circles mark the full detected-source list.
- Orange squares mark the track-excluded comparison list.
- The plotted line is the visible part of the comet track that falls inside the rendered mosaic.

A mismatch between the soft and hard source patterns is normal. A large mismatch between the track and the visible exposure footprint is not.

## `03_image.png` and `03_image_<band>.png`

These panels are comet-frame image diagnostics.

Top row:

- combined counts
- combined exposure
- combined rate
- combined masked-rate image

Bottom row:

- PN rate
- MOS1 rate
- MOS2 rate
- trail mask or a compact numerical summary

The masked-rate panel and the trail-mask panel often show streaks. In the comet frame those streaks usually represent stationary field sources trailing across the moving-target image plane. That is expected if the moving-target transform is behaving sensibly. The key question is whether those trails are being masked in the science products rather than forcing excessive GTI cuts.

## `04_contam.png`

This plot is a policy comparison, not just a pass/fail chart.

- `full` shows the exposure available if no contamination veto is applied.
- `src` rejects only source-aperture overlaps.
- `bkg` rejects only background-aperture overlaps.
- `strict` rejects either type.

For a diffuse coma, `full` plus a spatial trail mask is often the most informative starting point.

## `00_clean.png`

Use this to see which exposures lost a large fraction of time during flare cleaning. If one late exposure is heavily cut while the others are healthy, that is a sign that soft-proton filtering is doing work rather than a reason by itself to discard the entire observation.

## `07_spotchecks_<band>.png`

These are the most direct moving-target checks. The comet should stay approximately fixed between segments while stationary field sources trail.
