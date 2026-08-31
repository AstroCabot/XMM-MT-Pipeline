# Analysis 05: radial profile

This task fits the unsmoothed 0.35--1.10 keV EPIC counts in 48 radial rings.
The likelihood is Poisson in raw counts with the QPB, soft-proton, and PN
out-of-time components fixed. Twelve leave-one-sector-out fits set the slope
error. A separate moving-target PN fit profiles a hard-band particle template.

Aperture luminosities integrate the continuous broken-law comet term, excluding
the fitted flat sky. The unresolved center is regularized inside 2 arcsec (half
an image pixel).

Run:

    python -m analysis_05_radial_profile.run

`profile.tsv` contains the 48-ring EPIC profile; `annuli.tsv`, the ten chemistry
annuli; `pps_profile.tsv`, the PPS profile; `apertures.tsv`, the aperture
luminosities; and `result.json`, the fitted values and 370-arcsec
profile-to-spectrum luminosity ratio.
