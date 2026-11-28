# Analysis 05: radial profile

One continuous broken power law plus flat sky fits the reduced 0.35–1.10 keV
EPIC profile, including all 48 rings from 0 to 12 arcmin.
The Poisson likelihood uses raw counts, with the QPB,
soft-proton and pn out-of-time components fixed. The fitted parameters are the
normalization at the break, inner slope, outer slope, break radius and sky.
The inner power law extends to zero radius, with no flat core or central
exclusion. The pixel centered on the nucleus uses the finite integral over
its square area; aperture integrals also start at zero (the inner slope is
constrained below 2). This is an empirical fit, not a PSF deconvolution.

The same fitted comet component supplies aperture luminosities, with the
Task 03 conversion. Twelve leave-one-sector-out refits of the same broken model
give the jackknife errors on both slopes. The saved `slope` and
`sector_standard_error` fields refer to the outer slope. The saved
`normalization_ct_s_arcmin2` is the outer branch's amplitude at 1 arcmin.

The PPS comparison reads Task 02's saved 0.35–1.1 keV broken-profile fit, which
jointly profiles its hard-band particle template.

Inputs are the Task 04 `image.fits` and `result.json`, the Task 03
`result.json`, and the Task 02 outputs.

    python -m analysis_05_radial_profile.run
    python -m analysis_05_radial_profile.test

`profile.tsv` contains the EPIC profile; `annuli.tsv`, the chemistry annuli;
`pps_profile.tsv`, the PPS comparison; `apertures.tsv`, the EPIC aperture
luminosities; and `result.json`, the parameters and 370-arcsec
profile-to-spectrum luminosity ratio. `profiles.png` and `profiles.pdf` compare
the fits, in each data set's own count-rate units.

Synthetic tests recover both slopes, break and amplitude, check the joint
soft/hard Poisson nuisance solution, and compare the aperture integral with
numerical quadrature. The run also checks that the morphology exposure uses
the current spectral conversion.
