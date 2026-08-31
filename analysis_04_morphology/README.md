# Analysis 04: morphology

Inputs are the Analysis 01 soft-band stacks and Analysis 03 detector conversion
factors. Run in a sourced SAS environment:

```bash
export LD_LIBRARY_PATH=/home/scabot/miniconda3/envs/xmm/lib:${LD_LIBRARY_PATH:-}
/home/scabot/miniconda3/envs/xmm/bin/python run.py
```

It combines the
soft-band detector stacks in MOS2-equivalent exposure units, writes the raw count
and background planes used by Task 05, and creates the display-only adaptively
smoothed image, Sun-axis slice, and Figure 1.

The display mask enters `binadapt`, and the reported contour extent is the
component connected to the comet center.
Outputs are `result.json`, `image.fits`, `slice.tsv`, and
`figure1.png`/`figure1.pdf`.
