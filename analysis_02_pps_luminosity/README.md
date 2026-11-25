# Analysis 02: preliminary PPS luminosity

Run before the spectral analysis:

    python -m analysis_02_pps_luminosity.run

Inputs are the moving-target PPS pn event list, per-CCD GTIs from IMAGE_2000,
and EXPMAP1000/2000/3000. Task 01's settings file only
supplies the PPS directory. No custom-reduced events, source masks, background
maps or spectral-fit parameters are read. The normal run reads the saved
standard-response conversion factors in `ecfs.json`; SAS/XSPEC are not required.

The mean counts extend the submitted manuscript's Equation 4:

    soft_i = A sum_pixels(exposure * broken_shape * pixel_area)
             + B sum_pixels(exposure * pixel_area) + kappa h_i
    hard_i = h_i

The continuous broken power law has independent inner and outer slopes and a
fitted break radius. Its amplitude is defined at the break. Both count bands
enter a joint Poisson likelihood; each latent hard expectation h_i is profiled
analytically. All 48 annuli from 0 to 12 arcmin are fitted, including the center.
The observed inner profile is empirical, not an intrinsic PSF-deconvolved coma.
The inner power law extends to zero, with no flat core or central exclusion;
the aperture integral starts at zero, as in Task 05.

Independent fits use 0.2–0.5, 0.5–1.0, 1.0–2.0, 0.2–0.35, 0.2–0.30 and
0.35–1.1 keV. Selection is FLAG=0,
RAWY>12, native CCD GTIs, and PI in (low,high]. Bands wholly below 0.5 keV use
PATTERN=0; the others use PATTERN<=4, matching their conversion factors. The
hard template uses 2.5–8.5 keV and PATTERN<=4. The stricter FLAG cut differs
from the pipeline image bitmask. Field sources and out-of-time events are not
separately removed; these remain systematic limitations of this check.

For each band, integrate only the fitted comet term over a complete circle:

    rate(<R) = 2 pi A integral_0^R broken_shape(r) r dr
    flux = rate / K
    luminosity = 4 pi D^2 flux

Detector coverage enters the exposure-weighted likelihood. This replaces the
former annulus subtraction and its assumed 1/r self-subtraction correction.
K is the predicted count rate divided by model energy flux in the same band,
in counts cm²/erg, calculated directly from standard responses as below.

`parameters.json` specifies distance, apertures, bands, and sums. Default
radii are 110,000, 120,000, 135,000 and 190,000 km, 370 arcsec and 480 arcsec.
At the adopted distance, 110,000 km is 80.61 arcsec. The 0.2–1.0
and 0.2–2.0 keV luminosities sum the native-band fluxes; overlapping estimates
are not independent. The 1–2 keV component is weak and has a much flatter fitted
profile, so it warrants particular background scrutiny. No formal luminosity
error is claimed from these preliminary fits.

The 0.35–1.1 keV fit provides Task 05's profile comparison. Its exposure map
is an observed-rate-weighted combination of the native maps. The two narrow
bands reuse EXPMAP1000 (0.2–0.5 keV). These are approximate energy weightings.
Each band is fitted independently with the same spatial model; narrow-band
luminosities come from their own counts, not a fraction of a broader-band flux.

For every band, `ecfs.json` supplies one fixed K per band.
`calibrate.py` folds a fixed XSPEC `phabs*powerlaw` (Gamma=1.7, NH=3e20 cm^-2,
`angr` abundances, `bcmc` cross sections) through the off-the-shelf pn thin-filter
RMF/ARF in the [HEASARC AO26 bundle](https://heasarc.gsfc.nasa.gov/FTP/xmm/data/responses/).
The single-event bands use the canned epoch-e4/full-frame/Y9 singles RMF and
the CCF `EPN_QUANTUMEF_0018` single/double pattern-efficiency ratio applied to
the standard ARF in memory. All calibration inputs are public standard files;
no observation-specific response is generated or borrowed from Stage 03.
No comet spectrum is fitted. The model spans the response's full energy range,
including photons outside the selected PI band. K is predicted band count rate
divided by absorbed band energy flux; no absorption correction is applied to L.
The response uses 5-eV PI channels; calibrated EBOUNDS energies are not used to
select columns. Fractional edge bins approximate the integer-PI (low, high] cut.

The standard responses, including their point-source aperture losses, are used
directly; no catalogue normalization or additional aperture correction is applied.
Their use for the spatially extended comet is approximate.
The assumed power law describes photon energies, not the radial profile, and
is a reference conversion spectrum rather than a physical comet model.
Especially below 0.35 keV, these are spectrum-dependent preliminary estimates.
Broad-band luminosities sum the nonoverlapping native energy bands;
overlapping narrow-band estimates are never added to those sums.

To regenerate the fixed ECFs (requires PyXspec and downloads standard calibration
only; neither SAS nor the paper's spectral-fit outputs are used):

    python -m analysis_02_pps_luminosity.calibrate

The response-file hashes, model assumptions, predicted rates, fluxes and ECFs are saved in
`ecfs.json`. Normal profile reruns reuse this small file without response files.

Outputs: `result.json` (parameters, calibration, input hashes and luminosities),
`profiles.tsv` (observed and modeled ring counts), `luminosities.tsv` (aperture
table), and `samples.npz` (compact radial sufficient statistics for the demo).
Task 05 reads the saved 0.35–1.1 keV fit instead of fitting PPS a second time.
The shared likelihood and integral live in Task 05's `fit.py`; importing that
module does not read any Task 03–05 products.

    python -m analysis_02_pps_luminosity.test

Synthetic recovery tests for the joint Poisson fit are in Task 05's test module.
