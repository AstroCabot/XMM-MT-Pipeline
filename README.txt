XMM-Newton diffuse moving-target pipeline, reduction_v6
======================================================

Overview
--------

This pipeline is being built incrementally. Source files stay limited to:

  pipeline.sh
  qc.sh
  tools.py
  config.json
  README.txt

Generated products live under output/.

Python QC plots require astropy, numpy, and matplotlib.


General execution model
-----------------------

Run a stage:

  ./pipeline.sh <stage>

Run QC for a stage:

  ./qc.sh <stage>

Rerun a stage from scratch:

  ./pipeline.sh <stage> --force

Use a different config:

  ./pipeline.sh <stage> --config other_config.json

Choose detectors in config.json:

  "detectors": "PN"

Allowed values are PN, M1, M2, any comma/space separated combination, or all.
The default is PN.

Print the saved SAS environment after initialization:

  ./pipeline.sh env


Stage: init
-----------

Purpose:

  Prepare the SAS ODF state used by all later stages.

What it does:

  1. Sets SAS_ODF to the raw ODF directory.
  2. Uses the configured or inherited SAS_CCFPATH.
  3. Runs cifbuild in output/init.
  4. Sets SAS_CCF to output/init/ccf.cif.
  5. Runs odfingest with odfdir=<raw ODF> and outdir=output/init.
  6. Symlinks decompressed ODF files into output/init.
  7. Rewrites the generated *SUM.SAS PATH record to output/init.
  8. Writes output/sas_setup.env.

Outputs:

  output/init/ccf.cif
  output/init/*SUM.SAS
  output/init/<raw ODF symlinks>
  output/logs/init_cifbuild.log
  output/logs/init_odfingest.log
  output/sas_setup.env

Run it:

  cd /mnt/c/Users/cabot/OneDrive/Documents/Research/Cambridge/X3I/reduction_v6
  ./pipeline.sh init

QC:

  ./qc.sh init

QC outputs:

  output/qc/init/output_files.txt
  output/qc/init/file_type_counts.txt


Stage: repro
------------

Purpose:

  Reprocess EPIC ODF data into calibrated event lists and attitude products.

What it does:

  1. Sources output/sas_setup.env.
  2. Runs epproc and/or emproc in output/repro for the selected detectors.
  3. Writes raw event-list manifests for the selected detectors.
  4. Runs atthkgen to create output/repro/atthk.dat.
  5. Writes output/attitude.env for later moving-target stages.

Outputs:

  output/repro/*ImagingEvts.ds
  output/repro/manifest/PN_raw.txt
  output/repro/manifest/M1_raw.txt
  output/repro/manifest/M2_raw.txt
  output/repro/atthk.dat
  output/attitude.env

Run it:

  ./pipeline.sh repro

QC:

  ./qc.sh repro

QC outputs:

  output/qc/repro/output_files.txt
  output/qc/repro/manifest_counts.txt
  output/qc/repro/status.txt
  output/qc/repro/soft_lt1kev_mosaic.png
  output/qc/repro/hard_gt1kev_mosaic.png
  output/qc/repro/mosaic_summary.txt

The repro mosaics use the selected detectors. They split the event lists into
PI < 1000 eV and PI > 1000 eV bands.


Stage: clean
------------

Purpose:

  Apply standard EPIC event-quality and energy cuts. No espfilt or flare-GTI
  filtering is used.

What it does:

  1. Sources output/sas_setup.env and output/attitude.env.
  2. Filters PN with FLAG==0, PATTERN<=4, and the configured PI range.
  3. Filters MOS with #XMMEA_EM, PATTERN<=12, and the configured PI range.
  4. Runs attcalc on each cleaned event list.
  5. Writes per-instrument clean manifests and merged clean event files for the
     selected detectors.

Outputs:

  output/clean/PN/*.clean.fits
  output/clean/M1/*.clean.fits
  output/clean/M2/*.clean.fits
  output/clean/PN_clean_files.txt
  output/clean/M1_clean_files.txt
  output/clean/M2_clean_files.txt
  output/clean/PN_clean_merged.fits
  output/clean/M1_clean_merged.fits
  output/clean/M2_clean_merged.fits

Run it:

  ./pipeline.sh clean

QC:

  ./qc.sh clean

QC outputs:

  output/qc/clean/output_files.txt
  output/qc/clean/manifest_counts.txt
  output/qc/clean/status.txt
  output/qc/clean/soft_lt1kev_mosaic.png
  output/qc/clean/hard_gt1kev_mosaic.png
  output/qc/clean/mosaic_summary.txt

The clean mosaics are quick-look images made from the selected cleaned event
lists.


Stage: exposure
---------------

Purpose:

  Build matched EPIC counts, vignetted exposure, and count-rate maps from the
  cleaned event lists.

What it does:

  1. Sources output/sas_setup.env and output/attitude.env.
  2. Builds one shared SKY grid from the selected merged clean event lists.
  3. For each configured band, runs evselect to make selected-detector counts
     images.
  4. Runs eexpmap on those images using output/repro/atthk.dat.
  5. Combines available selected-detector counts and exposure maps into EPIC
     maps.

Outputs:

  output/exposure/grid.env
  output/exposure/grid.json
  output/exposure/soft/PN_counts.fits
  output/exposure/soft/PN_exposure.fits
  output/exposure/soft/EPIC_counts.fits
  output/exposure/soft/EPIC_exposure.fits
  output/exposure/soft/EPIC_rate.fits
  output/exposure/hard/<same pattern>

Run it:

  ./pipeline.sh exposure

QC:

  ./qc.sh exposure

QC outputs:

  output/qc/exposure/output_files.txt
  output/qc/exposure/status.txt
  output/qc/exposure/soft_exposure.png
  output/qc/exposure/soft_rate.png
  output/qc/exposure/hard_exposure.png
  output/qc/exposure/hard_rate.png

The default bands are soft=200-999 eV and hard=1001-12000 eV, matching the
PI < 1 keV and PI > 1 keV quick-look split without double-counting PI=1000.


Configuration
-------------

config.json currently defines:

  workdir
  odfdir
  shortlink_name / keep_shortlink / shortlink_max_path
  sas_setup_script
  sas_ccfpath
  sas_verbosity
  detectors
  ahf_input
  clean_pi_min / clean_pi_max
  premerge_image_size_deg
  image_bands
  image_bin_phys / image_pad_frac
  eexpmap_attrebin
  link_odf_constituents / skip_odf_link_patterns

sas_setup_script should point to the SAS setsas.sh file. If it is blank, the
pipeline will also try $SAS_DIR/setsas.sh when SAS_DIR is set.

No SAS warning suppression is set. SAS_RAND_SEED is not set; the current
stages do not need it.


Downstream contract
-------------------

Later stages should source:

  source output/sas_setup.env

That provides:

  SAS_CCF=output/init/ccf.cif
  SAS_ODF=output/init/*SUM.SAS
  SAS_CCFPATH=<configured or inherited value, if available>


References
----------

  ESA SAS User Guide, sections 2.3.1-2.3.4
  https://xmm-tools.cosmos.esa.int/external/xmm_user_support/documentation/sas_usg/USG/

  odfingest
  https://xmm-tools.cosmos.esa.int/external/sas/current/doc/odfingest/node3.html

  cifbuild
  https://xmm-tools.cosmos.esa.int/external/sas/current/doc/cifbuild.pdf

  eexpmap
  https://xmm-tools.cosmos.esa.int/external/sas/current/doc/eexpmap/
