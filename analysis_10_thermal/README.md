# Analysis 10: nucleus heating

`run.py` applies the 3I/ATLAS thermal model of Yaginuma et al. (2026) to the H2
budget and writes `outputs/result.json`. Heat conducts through 40 m of ice with
density 500 kg m^-3, heat capacity 550 J kg^-1 K^-1, and conductivity 1e-4,
1e-3, or 1e-2 W m^-1 K^-1, or the 0.1 W m^-1 K^-1 of Bergner & Seligman (2023),
under sunlight absorbed with Bond albedo 0.1 and rotation factor 1/pi, at a
fixed surface without sublimation. The ice starts at 20.9 K, 6500 days before
perihelion. With their density of 1000 kg m^-3, depth of 10 m, and emissivity
0.9, `test.py` reproduces their published peak temperatures and 80 and 150 K
depths; the paper states 0.95.

`water` integrates the water and H2 production of both Task 08 fits from
perihelion to XMM-Newton, and gives the depth of surface ice that holds that
water for an ice mass fraction of 0.5. It also gives the water mass produced
9--41 days after perihelion and, at 3.5 days, the H2 mass fraction of the H2O,
CO2, CO, and H2 production and the factor by which H2 raises the terminal speed
of that gas at fixed temperature.

`depth_reached_m` gives the deepest ice heated past 35 and 150 K by the
XMM-Newton epoch. `conduction` gives the thickest layer first heated past any
one temperature from 16 to 150 K between perihelion and XMM-Newton, its water,
and the share of the median H2 production it holds at 0.4 H2 per H2O. This
bounds the H2 that heating releases, whatever the release temperatures. It also
gives the share held by all ice heated past 35 K at 0.3 and 0.4 H2 per H2O.
`sunlight_H2O_molecules` is the most water that sunlight absorbed by the nucleus
could sublimate over the same interval.

Inputs are the Task 08 orbit and XMM-Newton epoch, the `chain.npz` draws of both
fits, and the Task 09 nucleus area. Run: `python run.py`.
