# 3I/ATLAS X-ray analysis demo

End-to-end demonstration of the post-reduction analysis pipeline
(stages 02–09 of the parent repository) from a ~11 MB self-contained data
bundle. This demo imports the analysis code from the parent repository.

## Run

```bash
python demo.py
```

Total runtime (several minutes) is dominated by the PN spectral 
fit (PyXspec + ACX2) and the reduced MCMCs.

More information on each stage is available on the main repo's 
GitHub page.
Deviations from the paper pipeline: PN only (no MOS); the PPS radial
cross-check and SPHEREx comparison are omitted; ESAS binadapt smoothing is
replaced by a scipy equivalent; the MCMC runs 48 walkers x 800 steps on a
coarsened grid.
