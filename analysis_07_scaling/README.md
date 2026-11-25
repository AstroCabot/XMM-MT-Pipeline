# Analysis 07: physical scaling

This task evaluates Wegmann et al. (2004) Equation 30 within 110,000 km, the
Cabot et al. (2023) Appendix B charge-exchange model, and the Wegmann outer-coma
thin limit.

Run:

    python -m analysis_07_scaling.run

`result.json` contains the Wegmann estimates, the Cabot 370-arcsec inversion,
propagated SPHEREx parent rates, pure-ice comparisons, and the 6-arcmin H2 thin
result. `scaling.tsv` contains the Wegmann, aperture and end-member table. Shared SPHEREx and H2 inputs are in
`../parameters.json`.

The SPHEREx rates cover 8--15 December. Their 2.0699 au reference distance is
the JPL value at the adopted midpoint, 2025-12-11 12:00 UTC.

The 110,000 km luminosity integrates the Task 05 broken-law comet term and
converts it with the Task 03 MOS2 ECF, whose log uncertainty it carries. The
370-arcsec luminosity is Task 03's full-circle spectral value. The Cabot
inversion uses 1 km/s; pure-ice comparisons use 0.37 km/s. B12 end-member
luminosities use the empirical heavy-ion energy flux and are total-model values
compared with the 110,000 km measurement. The luminosity-ratio uncertainty adds
the luminosity, C and flux log uncertainties in quadrature, with the nucleus
radius, speed and ice properties held fixed.

Primary equations: Cabot et al. (2023), B2--B9 and B11--B12; Wegmann et al. (2004),
10--11 and 30. Pure-ice rows use Cabot's generic cross section; the H2 thin
calculation uses the Machacek et al. (2014) O6+ single-capture measurement. B12
rows use Cabot's nominal `lambda=5e10 cm`; speed is specified separately.
The thin calculation uses the fitted outer power law. End-member comparisons use
a 1.3 km radius.
Production deficits are floored at zero.

`wegmann2004_fig19.tsv` holds the 48 points of Wegmann et al. (2004) Fig. 19:
H G^2 (1e54 erg s^-3 cm^-2) and the best-fit model luminosity (1e16 erg s^-1).
The log10 uncertainty of C = 1e-38 s^2 cm^2 is their root-mean-square offset from
that relation. Equation 30 is evaluated with three heavy-ion energy fluxes: the
Wegmann fast-wind value, 1.10e-4 erg cm^-2 s^-1; the same scaled by r_h^-2; and
f n u dE from our measurements. The measured flux uses the log-mean of the
STEREO-A, Solar Orbiter and Wind median proton fluxes (Task 06), the log-mean of
the ACE and Solar Orbiter median ion fractions with half their 16th-84th
percentile log range as uncertainty, and the Task 03 band-mean photon energy.
Rate ranges add the luminosity, C and flux log uncertainties in quadrature and
halve the sum.
