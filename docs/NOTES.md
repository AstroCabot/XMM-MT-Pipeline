# Notes

Measured against the previous compact bundle:

- files: 27 -> 25
- code lines (`*.py` + `*.sh`): 7313 -> 7096
- shipped text lines: 7589 -> 7340

Targeted QC-wrapper reduction:

- QC wrapper/helper files: 3 -> 2
- QC wrapper/helper lines: 493 -> 282

What changed:

- replaced the three QC shell/helper layers with one Python driver plus a thin shell wrapper
- folded trail-mask shared utilities into `scripts/trail_mask_tools.py`
- kept the reducer entry points and output layout unchanged
