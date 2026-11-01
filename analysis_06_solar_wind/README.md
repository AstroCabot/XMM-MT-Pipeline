# Analysis 06: solar wind and PN light curve

## Run

```bash
export LD_LIBRARY_PATH=/home/scabot/miniconda3/envs/xmm/lib:${LD_LIBRARY_PATH:-}
/home/scabot/miniconda3/envs/xmm/bin/python run.py
```

Inputs are the cached solar-wind series and the eight masked PN event and OoT
files, GTIs, ephemeris, attitude, calibration index, and ODF summary, plus the
ACE/SWICS and Solar Orbiter/HIS composition CDFs and the Task 03 VACX2 fit.

Outputs are `mapped_flux.tsv`, `lightcurve.tsv`, `result.json`, and Figure 4 as
PNG and PDF.

## Method

The hourly inputs are STEREO-A `STA_COHO1HR_MERGED_MAG_PLASMA`, Solar Orbiter
`SOLO_COHO1HR_MERGED_MAG_PLASMA`, ACE `AC_K1_SWE`, and Wind
`OMNI2_H0_MRG1HR`; OMNI spacecraft code 52 selects Wind. Horizons tracks use
3I/ATLAS, STEREO-A (`-234`), Solar Orbiter (`-144`), ACE (`-92`), and Earth
(`399`). OMNI times already represent arrival near Earth.

Each monitor sample is mapped by

`t_a = t_s + wrap180(lambda_c(t_a)-lambda_m(t_s))/14.1844 + (r_c(t_a)-r_m(t_s))*AU/v`

and its flux is scaled by `(r_m/r_c)^2`. Gaps over 1.5 hours remain empty;
segments are sampled every 10 minutes and averaged in a centered 120-minute
window. ACE density uses a fixed factor of 1.7 relative to Wind; the unscaled
series is omitted from the relative panel. The Parker panel uses the nearest
STEREO-A hourly sample within the 1.5-hour gap limit.

F10.7 is the Penticton 20:00 UT value adjusted to 1 au, interpolated to the
observation midpoint minus the solar-rotation lag from Earth to the comet,
`lag = wrap180(lambda_c(t)-lambda_E(t-lag))/14.1844` days, so that it describes
the hemisphere facing the comet. `F10.7A` is the centered 81-day mean and
`F10.7P=(F10.7+F10.7A)/2`; Task 08 reads F10.7P from `result.json`.

`ion_fraction.py` tabulates C5+,6+, N6+,7+, O7+,8+, Ne9+, and Mg11+ relative to
protons for ACE/SWICS and Solar Orbiter/HIS, using the same speed-dependent O/H,
0.501-dex FWHM lognormal O/H scatter, reconstructed Ne9+, and tabulated N6+,7+.
The 16th--84th percentiles are `(1.39--4.40)e-4` for ACE and `(1.13--3.55)e-4`
for HIS; the adopted ion-fraction range is `1e-4--5e-4`. Its limits times the
lowest and highest mapped proton flux across the monitors in the observation
window, with ACE scaled by 1.7, give the bounds of the uniform heavy-ion
number-flux prior, `4.3e3--1.6e5` cm^-2 s^-1, which Tasks 07 and 08 read from
`result.json`. These ions supply 97.3% of the intrinsic 0.35--1.10 keV energy
flux in the Task 03 VACX2 model. ACE CDFs come from SPDF and HIS V04 CDFs from
SOAR; `result.json` records both input digests. HAPI cannot replace them because
its CSV rounds the ratio floors and omits the HIS quality flag.

The light curve selects 200--1000 eV PN singles (`PATTERN==0`, `FLAG==0`). The
15 common 45-arcsec source holes apply to events, OoT events, and exposure. SAS
`eexpmap` is unvignetted (`VIGNET=F`) and includes filter transmission, quantum
efficiency, GTIs, gaps, and bad pixels. For every frame,

`rate = A370[(Nsrc-0.063 Osrc)/Esrc - (Nbg-0.063 Obg)/Ebg]`,

where `A370=pi(370/60)^2` arcmin2 and the background annulus is 480--720
arcsec.

Primary documentation: [CDAWeb HAPI](https://cdaweb.gsfc.nasa.gov/hapi),
[OMNI timing](https://omniweb.gsfc.nasa.gov/html/sc_merge_data1.html),
[Horizons](https://ssd.jpl.nasa.gov/horizons/manual.html),
[Penticton F10.7](https://www.spaceweather.gc.ca/forecast-prevision/solar-solaire/solarflux/sx-5-en.php),
[SAS eexpmap](https://xmm-tools.cosmos.esa.int/external/sas/current/doc/eexpmap/node3.html),
[Schwadron & Cravens (2000)](https://doi.org/10.1086/317176), and
[Koutroumpa et al. (2025)](https://doi.org/10.1029/2024GL114374), with the
[Solar Orbiter HIS archive](https://cdaweb.gsfc.nasa.gov/misc/NotesS.html).
