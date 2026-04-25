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
  output/cheese/inputs/<detector>/<event>/*_cheese_events.fits       (PI 1000-4000)
  output/cheese/inputs/<detector>/<event>/*_cheese_counts.fits
  output/cheese/inputs/<detector>/<event>/*_cheese_exp_vig.fits
  output/cheese/inputs/<detector>/<event>/*_cheese_bkg.fits          (esplinemap)
  output/cheese/inputs/<detector>/<event>/*_cheese_eboxlist_{l,m}.fits
  output/cheese/regions/<detector>/<event>/*_hard_detmask.fits       (emask)
  output/cheese/regions/<detector>/<event>/*_hard_emllist.fits       (emldetect)
  output/cheese/regions/<detector>/<event>/*_hard_global_region.fits (region sky)
  output/cheese/regions/<detector>/<event>/*_hard_detxy_region.fits  (region detxy)
  output/cheese/regions/<detector>/<event>/*_hard_regmask.fits       (regionmask)
  output/cheese/regions/<detector>/<event>/*_hard_base_mask.fits     (combined)
  output/cheese/<band>/<detector>/<event>/*_<band>_cheesemask.fits
  output/cheese/<band>/<detector>/<event>/*_<band>_cheesed_counts.fits
  output/cheese/<band>/<detector>/<event>/*_<band>_cheesed_exp_vig.fits

cheese mirrors the reduction_v8 cheese stage: it builds its own detection
images at PI [cheese_emin, cheese_emax] (default 1000-4000 eV) by filtering
the clean-stage event files in that range, binning onto the maps grid, and
running eexpmap + emask. It then performs v8's full source-detection
sequence:

  eboxdetect usemap=no   local box pass with likemin=cheese_mlmin, boxsize=5,
                         nruns=3, over the cheese band.
  esplinemap             spline background with nsplinenodes=20,
                         excesssigma=4, nfitrun=3, snrmin=30, smoothsigma=15,
                         scut=0.01, mlmin=1.
  eboxdetect usemap=yes  map-mode box pass seeded by the spline background.
  emldetect              ellbeta PSF with beta-model extent
                         (fitextent=yes, extentmodel=beta, dmlextmin=6,
                         minextent=1.5, maxextent=20, withtwostage=yes).
  region x2              outunit=detxy and outunit=xy on the emllist, with
                         expression (ID_INST==N)&&(ID_BAND==1)&&
                         (DET_ML>=MLMIN.0)&&(FLUX>=FLUX*1e-14),
                         shrinkconfused=yes, radiusstyle=contour,
                         fixedradius=12, energyfraction=0.9,
                         bkgfraction=0.6, withboresightfudge=yes.
  regionmask             rasterize the sky region into a regmask on the
                         cheese-band counts grid.
  combine-masks          base_mask = (detmask > 0) & (regmask > 0)
                         (python helper in tools.py, v8-identical).
  per-band               copy base_mask to each maps band's cheesemask and
                         multiply counts and vignetted exposure by it.

If emldetect (or a subsequent region call) fails, cheese degrades gracefully
by copying the detmask into base_mask — the frame still gets a FOV keep-mask
so downstream bands are not blocked. This matches the v8 fallback behaviour.

PN is instid=1. For MOS1/MOS2 the stage sets instid=2 and instid=3 and
switches to hrm1def/hrm2def and xidm1def/xidm2def automatically; the
detectors list in config.json controls which are processed.

Relevant config keys (all optional; listed with defaults):

  cheese_emin                 1000    cheese detection band lower PI
  cheese_emax                 4000    cheese detection band upper PI
  cheese_mlmin                10      emldetect mlmin and region DET_ML floor
  cheese_flux                 0.2     region FLUX floor in units of 1e-14
  cheese_ecf_pn               3.2     energy conversion factor (pn)
  cheese_ecf_mos              1.65    energy conversion factor (mos)
  cheese_emask_thresh1        0.3     emask low threshold
  cheese_emask_thresh2        0.5     emask high threshold
  cheese_region_energyfraction 0.9    region energyfraction
  cheese_region_fixedradius   12      region fixedradius (sky pixels)
  cheese_region_bkgfraction   0.6     region bkgfraction
  cheese_esplinemap_nsplinenodes 20
  cheese_esplinemap_excesssigma  4
  cheese_esplinemap_nfitrun      3
  cheese_esplinemap_snrmin       30
  cheese_esplinemap_smoothsigma  15
  cheese_esplinemap_scut         0.01
  cheese_esplinemap_mlmin        1
  cheese_eboxdetect_nruns        3
  cheese_eboxdetect_boxsize      5

Note: cheese runs its own eexpmap for the cheese band (once per inst/base).
This is additional work relative to the original v6 cheese, but preserves
the v8 detection-chain fidelity -- emldetect / esplinemap / eboxdetect all
require an exposure map and detmask matched to the detection band. The
vignetted maps-stage exposures at 1000-2000 / 2000-7200 eV (soft / hard) are
reused unchanged for the cheesed_exp_vig products.

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
