# 3I/ATLAS reproducible analysis

Each numbered directory produces one paper result in `outputs/`. Shared constants
are in `parameters.json`.

| Task | Result |
|---|---|
| `analysis_01_reduction` | Masked EPIC events, ESAS backgrounds, soft-proton models, and comet-frame stacks |
| `analysis_02_pps_luminosity` | Canned-response PPS luminosity |
| `analysis_03_spectrum` | Joint PN/MOS VACX2 fit, detector conversion, and Figure 2 |
| `analysis_04_morphology` | Comet image, Sun-axis slice, and Figure 1 |
| `analysis_05_radial_profile` | EPIC/PPS density profiles, aperture luminosities, and chemistry annuli |
| `analysis_06_solar_wind` | Ballistically mapped solar wind, masked PN light curve, and Figure 4 |
| `analysis_07_scaling` | Cabot and Wegmann production-rate inversions |
| `analysis_08_chemistry` | Neutral-chain, opacity, PSF, and H2 posterior calculation |
| `analysis_09_recoil` | 3I/ATLAS and 1I/ʻOumuamua recoil and radiolytic capacity |

## System requirements

Python >= 3.12 with the pinned packages in `requirements.txt` (tested with
Python 3.14.2 and 3.12.4 on Ubuntu 24.04 under WSL2; no non-standard
hardware). Tasks 01–06 additionally require HEASoft 6.36 with PyXspec,
SAS 22.1.0 with current CCFs, the ESAS calibration files, AtomDB 3.1.3, ACX2
v2.5.1 (github.com/AtomDB/ACX2), and roughly 20 GB of disk.

## Installation guide

In a fresh virtual environment:

    pip install -r requirements.txt

About one minute on a normal desktop computer. HEASoft, SAS, and AtomDB
(full pipeline only) install per their own documentation in one to two hours.

## Demo

    python analysis_09_recoil/run.py

Runs in under one second on a normal desktop computer from the committed demo
inputs (the Task 08 posterior draws in `analysis_08_chemistry/outputs/`
and `analysis_09_recoil/parameters.json`), and needs only NumPy and SciPy. Expected output:
`analysis_09_recoil/outputs/result.json` is rewritten identical to the
committed copy, so `git diff` stays clean.

## Instructions for use

Run the numbered tasks in order; later tasks read the earlier `outputs/`.
For task 01, download the ODF, PPS, and Moving-Target PPS trees for
ObsID 0963720201 from the XMM-Newton Science Archive and point the paths in
`analysis_01_reduction/settings.json` at those trees and the local SAS,
HEASoft, ESAS, and AtomDB installations. Per-task details, including the
required environment sourcing, are in each task's `README.md`.

## License

MIT (see `LICENSE`).
