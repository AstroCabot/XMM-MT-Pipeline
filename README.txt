XMM-Newton diffuse moving-target pipeline, reduction_v6
======================================================

Scope
-----

The active source files are:

  pipeline.sh
  tools.py
  config.json
  README.txt

For now this version implements init, repro, clean, maps, and cheese.
Generated products live under output/.


Running
-------

Production:

  ./pipeline.sh init
  ./pipeline.sh repro
  ./pipeline.sh clean
  ./pipeline.sh maps
  ./pipeline.sh cheese
  ./pipeline.sh all

QC:

  ./pipeline.sh qc-init
  ./pipeline.sh qc-repro
  ./pipeline.sh qc-clean
  ./pipeline.sh qc-maps
  ./pipeline.sh qc-cheese
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
  Optional per-band exposure_pi_range can override the pimin/pimax passed to
  exposure and detection tasks when a broad or split band needs a different
  representative exposure energy.

map_bin_phys, map_pad_frac
  Shared sky-image grid controls for the maps stage. The default bin size is 80
  physical X/Y units, corresponding to 4 arcsec EPIC image pixels.

map_soft_hard_split_ev
  PI/eV boundary for maps production bands. maps only builds soft and hard
  products. The default is 1000, giving soft PI 201-1000 and hard PI
  1001-12000. map_pi_min and map_pi_max set the outer bounds.

map_eexpmap_attrebin
  eexpmap attitude rebinning parameter. The default is 0.020626481 arcsec, the
  SAS-documented value for sky exposure maps that match input images closely.

map_source_detection_bands
  Bands used to build the shared point-source list. The default is hard, i.e.
  the >1 keV maps band where comet charge-exchange emission should be weak.
  That hard-band source list is reused by cheese for every band, including the
  soft CX band.

map_ebox_*
  eboxdetect controls for the shared hard-band source list. The same threshold
  is used for the local and map-mode passes.

map_bkg_*
  Hard-band source-refinement background controls used inside maps. maps now
  runs local eboxdetect, then esplinemap on the configured source-detection
  band(s), then map-mode eboxdetect, and cheese uses that refined source list.


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


Maps
----

maps creates or adopts:

  output/maps/grid.env
  output/maps/grid_summary.tsv
  output/maps/maps_band_table.tsv
  output/maps/maps_manifest.tsv
  output/maps/source_manifest.tsv
  output/maps/events/<band>/<detector>/*_<band>_map_events.fits
  output/maps/<band>/<detector>/<event>/*_<band>_counts.fits
  output/maps/<band>/<detector>/<event>/*_<band>_exp_vig.fits
  output/maps/<band>/<detector>/<event>/*_<band>_exp_unvig.fits
  output/maps/sources/<detector>/<event>/*_source_local.fits
  output/maps/sources/<detector>/<event>/*_<band>_source_background.fits
  output/maps/sources/<detector>/<event>/*_source_list.fits

maps consumes the clean-stage event products under output/clean/events/ and
output/clean/manifest/. It does not rebuild map inputs from raw repro event
lists.

The stage now stops at the light-weight SAS source-detection preparation step:

  derive soft/hard map event lists from clean event products
  evselect counts image
  eexpmap withvignetting=yes
  eexpmap withvignetting=no
  eboxdetect local on the configured hard/source band
  esplinemap on the configured hard/source band
  eboxdetect map mode on the configured hard/source band

maps does not exclude the central comet zone. The point is to obtain the
shared hard-band source list and the per-band count/exposure products on a
common grid before any soft-band background modelling. The hard-band
background map is only used to improve point-source finding for cheese.

Maps can be rerun from a chosen dependency point without rebuilding upstream
products:

  ./pipeline.sh maps --maps-from grid
  ./pipeline.sh maps --maps-from counts
  ./pipeline.sh maps --maps-from exposure
  ./pipeline.sh maps --maps-from sources

Use grid after changing map_bin_phys, map_pad_frac, or clean products. Use
counts after changing map_soft_hard_split_ev, map_pi_min, map_pi_max, or the
clean event inputs. Use exposure after changing eexpmap energy or
attitude settings. Use sources after changing the hard source-detection bands
or local eboxdetect thresholds. The old mask keyword is still accepted as an
alias for sources.

qc-maps also writes PNG mosaics from maps_manifest.tsv:

  output/qc/maps/soft_counts_mosaic.png
  output/qc/maps/soft_exposure_vig_mosaic.png
  output/qc/maps/soft_exposure_unvig_mosaic.png
  output/qc/maps/hard_counts_mosaic.png
  output/qc/maps/hard_exposure_vig_mosaic.png
  output/qc/maps/hard_exposure_unvig_mosaic.png


Cheese
------

cheese creates or adopts:

  output/cheese/cheese_manifest.tsv
  output/cheese/regions/<detector>/<event>/*_hard_global_region.fits
  output/cheese/regions/<detector>/<event>/*_hard_base_mask.fits
  output/cheese/<band>/<detector>/<event>/*_<band>_cheesemask.fits
  output/cheese/<band>/<detector>/<event>/*_<band>_cheesed_counts.fits
  output/cheese/<band>/<detector>/<event>/*_<band>_cheesed_exp_vig.fits

cheese consumes output/maps/source_manifest.tsv and output/maps/maps_manifest.tsv.
It uses the shared hard-band source list to define source holes, and it does
not apply the central comet exclusion. The source holes are turned into a
per-band image mask by:

  region operationstyle=global
  regionmask
  copy the resulting keep-mask to each band
  multiply counts and vignetted exposure by the resulting cheese mask

This stage is intentionally lighter than a full SAS diffuse-background fit. It
is meant to provide masked images for later background estimation or smoothing
while reusing the hard-band-only source-refinement products already created by
maps.

cheese can be rerun from a chosen dependency point:

  ./pipeline.sh cheese --cheese-from regions
  ./pipeline.sh cheese --cheese-from mask
  ./pipeline.sh cheese --cheese-from images

Use regions after changing how the hard-band source list is turned into source
holes. Use mask after changing the per-band cheese-mask construction. Use images after
changing only the masked-image writing step.

qc-cheese also writes PNG mosaics from cheese_manifest.tsv:

  output/qc/cheese/soft_cheese_mask_mosaic.png
  output/qc/cheese/soft_cheesed_counts_mosaic.png
  output/qc/cheese/soft_cheesed_exposure_vig_mosaic.png
  output/qc/cheese/hard_cheese_mask_mosaic.png
  output/qc/cheese/hard_cheesed_counts_mosaic.png
  output/qc/cheese/hard_cheesed_exposure_vig_mosaic.png


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
  python3 tools.py maps-band-table config.json
  python3 tools.py maps-expr config.json PN soft
