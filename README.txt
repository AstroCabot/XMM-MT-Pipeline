XMM-Newton diffuse moving-target pipeline, reduction_v6
======================================================

Scope
-----

The active source files are:

  pipeline.sh
  tools.py
  config.json
  README.txt

Generated products live under output/. The PPS directory is reference input and
must not be edited by this pipeline.


Running
-------

Run a production stage:

  ./pipeline.sh init
  ./pipeline.sh repro
  ./pipeline.sh clean
  ./pipeline.sh exposure

Run QC from the same entry point:

  ./pipeline.sh qc-init
  ./pipeline.sh qc-repro
  ./pipeline.sh qc-clean
  ./pipeline.sh qc-exposure
  ./pipeline.sh qc

Run matching QC after a production stage:

  ./pipeline.sh exposure --force --qc

Use --force only when intentionally rebuilding a stage.


Current Data Layout
-------------------

The current cleaned PN products on disk use this layout:

  output/clean/events/PN/*_clean.fits
  output/clean/manifest/PN_clean_files.txt
  output/clean/flare_summary.tsv
  output/clean/gti/PN/*_flare_gti.fits
  output/clean/lightcurves/PN/*_flare_lc.fits

The manifest is authoritative. Extra clean event files may exist on disk from
earlier trials, but downstream stages and QC use only the files listed in
output/clean/manifest/PN_clean_files.txt. Running clean without --force should
skip when the selected detector manifests exist and all listed event files are
present.

The current output/exposure directory contains earlier two-band trial products
under output/exposure/PN/. The repaired exposure stage does not treat those as
complete PPS-band products. Rebuilding exposure writes the PPS-band layout
described below.


Exposure Strategy
-----------------

The public images stage is intentionally absent. The exposure stage creates the
count images it needs internally.

The configured PN exposure bands follow the PPS band boundaries:

  1000  eexpmap 200-500 eV      counts PI 201-500, PATTERN==0
  2000  eexpmap 500-1000 eV     counts PI 501-1000, PATTERN<=4
  3000  eexpmap 1000-2000 eV    counts PI 1001-2000, PATTERN<=4
  4000  eexpmap 2000-4500 eV    counts PI 2001-4500, PATTERN<=4
  5000  eexpmap 4500-12000 eV   counts PI 4501-7800 or 8201-12000, PATTERN<=4

The count filters also require RAWY>=13 and PPS-style FLAG masks.

For each band and slice, exposure builds:

  *_counts.fits
  *_exposure_ref.fits
  *_exposure.fits
  *_novig_exposure.fits

The *_exposure_ref.fits files are internal eexpmap reference images. They keep
the slice TIME and FLAG selections in DSS even when the science counts image is
empty in a narrow band.

The vignetted exposure map is for sky/comet photons. The non-vignetted exposure
map is retained for later background modeling support.


Stage Summary
-------------

init
  Builds output/init/ccf.cif, output/init/*SUM.SAS, ODF symlinks, and
  output/sas_setup.env.

repro
  Runs epproc/emproc for selected detectors and writes output/repro manifests
  plus output/repro/atthk.dat and output/attitude.env.

clean
  Applies event-quality and flare-GTI filtering into output/clean/events and
  writes output/clean/manifest. Existing clean products are detected from that
  manifest.

exposure
  Builds PPS-band counts, vignetted exposure, non-vignetted exposure, and rate
  maps from the cleaned event lists.

background, sources, stack
  Later stages are still present but should be treated as downstream work. They
  depend on the exposure products.


Important Config Keys
---------------------

detectors
  Currently PN.

image_bands
  PPS exposure bands, encoded as label:pimin:pimax entries.

image_bin_phys
  SKY image bin size. Current value is 80.

eexpmap_attrebin
  Passed to eexpmap.

source_detect_band
  Current downstream source-detection band. With PPS bands this is 4000.

clean_gti_*
  High-energy flare-GTI settings for clean when rebuilding. Current values
  match the existing PN flare-lightcurve provenance: PI 7000-15000 eV, 10 s
  bins, and 4.80001211 ct/s for the accepted S003 exposure.


Validation Commands
-------------------

These checks do not run SAS production tasks:

  bash -n pipeline.sh
  python -m py_compile tools.py
  python -m json.tool config.json
