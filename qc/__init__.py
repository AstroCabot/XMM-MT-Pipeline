"""XMM Moving-Target Pipeline — Quality-Control diagnostic package.

Modules
-------
checks      – Full diagnostic check suite (17 check_* functions).
manifest    – Stage-level manifest builder and QC index generator.
plot_utils  – Shared matplotlib helpers: figures, norms, colorbars, scalebars.
quick_diag  – Fast comet-frame time-bin diagnostics.
runner      – End-to-end QC orchestrator (calls checks/manifest/quickdiag).
utils       – Shared I/O helpers: env loading, FITS utilities, file readers.
"""
