# Analysis 07: physical scaling

This task evaluates the Cabot et al. (2023) Appendix B charge-exchange model,
an approximate Hyakutake-scale application of Wegmann et al. (2004) Equation 30,
and the Wegmann outer-coma thin limit.

Run:

    python -m analysis_07_scaling.run

`result.json` contains the three production-rate inferences, propagated SPHEREx
parent rates, pure-ice comparisons, and the 6-arcmin H2 thin result. `scaling.tsv`
contains the aperture/end-member table. Shared SPHEREx and H2 inputs are in
`../parameters.json`.

The SPHEREx rates cover 8--15 December. Their 2.0699 au reference distance is
the JPL value at the adopted midpoint, 2025-12-11 12:00 UTC.

The 120,000 and 135,000 km luminosities are Task 05 integrals of the comet term;
the 370-arcsec luminosity is Task 03's full-circle spectral value. The Cabot
inversion uses 1 km/s; pure-ice comparisons use 0.37 km/s. B12 end-member
luminosities are total-model values compared with the 370-arcsec measurement.

Primary equations: Cabot et al. (2023), B2--B9 and B11--B12; Wegmann et al. (2004),
10--11 and 30. Pure-ice rows use Cabot's generic cross section; the H2 thin
calculation uses the Machacek et al. (2014) O6+ single-capture measurement. B12
rows use Cabot's nominal `lambda=5e10 cm`; speed is specified separately.
The thin calculation uses the fitted outer power law. End-member comparisons use
a 1.3 km radius.
Production deficits are floored at zero.
The thin-limit calculation pairs the Cabot ion number flux with Task 03's fitted
mean photon energy. Its optical depth is speed-independent; the reported
production rate scales linearly between 1.0 and 1.2 km/s. The 370-arcsec result
also gives the heavy-ion-flux multiple required by B12 using the measured
parents.
