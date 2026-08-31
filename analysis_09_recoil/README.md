# Analysis 09: recoil

`run.py` reads Task 08 draws and writes the recoil and radiolytic-capacity
quantities to `outputs/result.json`.

The 3I calculation uses the radial Spada et al. `A1` coefficients, the Hui et al.
area product at `pV=0.04`, and thermal H2 recoil from Bergner & Seligman.
Conventional-parent speeds use the rowwise Task 08 posterior draws.
The zero-delay symmetric Spada solution is the headline; the 7.49-day solution
is a bracket.

The Oumuamua requirement uses its `pV=0.1` ellipsoid, density,
temperature, collimation, and radial acceleration. Force fractions use the central nucleus
mass; its propagated uncertainty is reported with each required force. The
area-scaled H2 comparison uses the Task 08 posterior interval.

The radiolytic rates divide each eroded-shell H2 inventory by the 249.791 days
between first reported activity and the XMM-Newton midpoint.

Run: `python run.py`.
