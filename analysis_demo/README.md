# 3I/ATLAS X-ray analysis demo

End-to-end demonstration of the post-reduction analysis pipeline
(stages 02–09 of the parent repository) from a ~19 MB self-contained data
bundle. This demo imports the analysis code from the parent repository.
Stage 02 reruns the PPS fits from bundled annular counts, exposure weights and
standard-response conversion factors; the full event list is not bundled.
This includes independent 0.2–0.35 and 0.2–0.30 keV profiles and the 110,000 km
aperture. No response files or spectral-fit outputs are needed for Stage 02.
Stage 06 uses the ion-fraction reference summaries in the shared parameters;
it does not repeat the ACE/SWICS or Solar Orbiter/HIS composition analysis.

## Run

```bash
python demo.py
```

Total runtime (several minutes) is dominated by the PN spectral 
fit (PyXspec + ACX2) and the reduced MCMCs.

More information on each stage is available on the main repo's 
GitHub page.
Deviations from the paper pipeline: PN only (no MOS); the SPHEREx comparison
is omitted; ESAS binadapt smoothing is
replaced by a scipy equivalent; the MCMC runs 48 walkers x 800 steps on a
coarsened grid.
