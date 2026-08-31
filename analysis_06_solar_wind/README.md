# Analysis 06: solar wind and PN light curve

## Run

```bash
export LD_LIBRARY_PATH=/home/scabot/miniconda3/envs/xmm/lib:${LD_LIBRARY_PATH:-}
/home/scabot/miniconda3/envs/xmm/bin/python run.py
```

Inputs are the cached solar-wind series and the eight masked PN event and OoT
files, GTIs, ephemeris, attitude, calibration index, and ODF summary.

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

F10.7 is the Penticton 20:00 UT value adjusted to 1 au. `F10.7A` is its centered
81-day mean and `F10.7P=(F10.7+F10.7A)/2`; Task 08 uses the mean F10.7P of the
two observation dates. The adopted ion-fraction range is `3e-4--1e-3`, from
Schwadron--Cravens, ACE/SWICS, and Solar Orbiter/HIS. The heavy-ion number-flux
prior is the zero-truncated `N(1.0e5, 3.8e4)` cm^-2 s^-1 distribution.

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
