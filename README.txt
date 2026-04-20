XMM-Newton diffuse moving-target pipeline, reduction_v6
======================================================

Scope
-----

The active source files are:

  pipeline.sh
  tools.py
  config.json
  README.txt

For now this version implements init, repro, and clean. Generated products live
under output/.


Running
-------

Production:

  ./pipeline.sh init
  ./pipeline.sh repro
  ./pipeline.sh clean
  ./pipeline.sh all

QC:

  ./pipeline.sh qc-init
  ./pipeline.sh qc-repro
  ./pipeline.sh qc-clean
  ./pipeline.sh qc

Use --force only when intentionally rebuilding products.


Config
------

qc_soft_hard_split_ev
  PI/eV threshold for the repro QC soft and hard mosaics. The current default is
  1000, so soft uses PI < 1000 and hard uses PI > 1000.

clean_gti_*
  Fixed high-energy flare-GTI settings. The current defaults use 7-15 keV,
  PATTERN <= 0, 10 s light-curve bins, and RATE <= 4.80001211 count/s.

clean_bands
  PPS-style spectral event-filter definitions used by clean. Each band defines
  PI ranges plus pn and MOS pattern, flag, and extra terms. The default pn band
  filters follow the PPS image-band convention:

    1000  0.2-0.5 keV   PI 201-500       PATTERN == 0   FLAG mask 0x2fb002c
    2000  0.5-1.0 keV   PI 501-1000      PATTERN <= 4   FLAG mask 0x2fb002c
    3000  1.0-2.0 keV   PI 1001-2000     PATTERN <= 4   FLAG mask 0x2fb0024
    4000  2.0-4.5 keV   PI 2001-4500     PATTERN <= 4   FLAG mask 0x2fb0024
    5000  4.5-12.0 keV  PI 4501-7800 or 8201-12000,
                        PATTERN <= 4, FLAG mask 0x2fb0024

  The default pn filters also include RAWY >= 13. The default MOS filters use
  #XMMEA_EM and PATTERN <= 12 for each configured PI band.


Init
----

init creates or adopts:

  output/init/ccf.cif
  output/init/*SUM.SAS
  output/init/<ODF constituent symlinks>
  output/sas_setup.env

The ODF constituent symlinks are intentional. The SUM.SAS PATH line is rewritten
to output/init/ so later SAS tasks read the mirrored ODF directory rather than
the original ODF directory. This preserves the old v6 convention and keeps the
working products relocatable within this reduction tree.

If ccf.cif and *SUM.SAS already exist, init writes the missing env file and
skips cifbuild/odfingest unless --force is set.


Repro
-----

repro creates or adopts:

  output/repro/*_ImagingEvts.ds
  output/repro/manifest/PN_raw.txt
  output/repro/manifest/M1_raw.txt
  output/repro/manifest/M2_raw.txt
  output/repro/atthk.dat
  output/attitude.env

The selected detectors are controlled by config.json. Accepted detector names
are PN, M1, M2, MOS1, MOS2, EMOS1, EMOS2, EPIC, and ALL.

Before running epproc or emproc, repro rebuilds manifests from existing event
lists. If the selected manifests, atthk.dat, and attitude.env are already valid,
the stage skips without requiring SAS in PATH.


Clean
-----

clean creates or adopts:

  output/clean/events/<detector>/<band>/*_<band>_clean.fits
  output/clean/manifest/<detector>_<band>_clean_files.txt
  output/clean/manifest/<detector>_clean_files.txt
  output/clean/clean_band_filters.tsv
  output/clean/gti/<detector>/*_flare_gti.fits
  output/clean/lightcurves/<detector>/*_flare_lc.fits
  output/clean/flare_gti_summary.tsv

One flare GTI is built per raw event list and applied to each configured band.
The aggregate manifest is kept for compatibility with older products; the
per-band manifests carry the configured PPS-style filters forward explicitly.


Shortlink Convention
--------------------

When the configured output path is longer than shortlink_max_path, the pipeline
runs SAS tasks through /tmp/_xmm_<shortlink_name>. Env files and manifests still
store stable real paths under output/.


Validation
----------

These checks do not run SAS production tasks:

  bash -n pipeline.sh
  python3 -m py_compile tools.py
  python3 -m json.tool config.json
  python3 tools.py clean-band-table config.json
