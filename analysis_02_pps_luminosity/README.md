# Analysis 02: PPS luminosity

`run.py` measures the masked 0.35--1.10 keV PN luminosity inside 370 arcsec.
It applies the native moving-target CCD GTIs, masks the fifteen stationary-source
trails, transfers the 600--720 arcsec background with the three PPS exposure
maps, builds a flat extended-source ARF, and fluxes each sub-band with the canned
PN RMF. `efluxer` processes the full response-supported PI range before
integration over physical-energy bands. `outputs/result.json`
names three flux/luminosity products: raw masked,
annulus-self-subtraction-corrected masked, and source-mask-filled full aperture.
It also records the annulus, mask-fill, and combined `1/rho` correction factors.
The versioned canned RMF is in `calibration/`.

From `analysis_repro`, run in a sourced SAS WSL environment:

```bash
export LD_LIBRARY_PATH=/home/scabot/miniconda3/envs/xmm/lib:${LD_LIBRARY_PATH:-}
/home/scabot/miniconda3/envs/xmm/bin/python -m analysis_02_pps_luminosity.run
```

Inputs are the Analysis 01 ephemeris and source catalog.
