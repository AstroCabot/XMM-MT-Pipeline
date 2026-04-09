#!/usr/bin/env python3
"""Compare atthkgen-produced ATTHK with attmove AHF to confirm consistency."""
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.coordinates import Angle
from astropy.time import Time
import astropy.units as u

base = Path("/mnt/c/Users/cabot/OneDrive/Documents/Research/Cambridge/X3I/reduction_v4/output")

# Read atthkgen ATTHK
with fits.open(str(base / "comet" / "test_atthkgen_atthk.dat")) as hdul:
    ext = hdul["ATTHK"]
    t_atthk = ext.data["TIME"].astype(np.float64)
    ra_atthk = ext.data["AHFRA"].astype(np.float64)
    dec_atthk = ext.data["AHFDEC"].astype(np.float64)

# Read AHF (sexagesimal)
with fits.open(str(base / "comet" / "moved_odf" / "odf" / "4759_0963720201_SCX00000ATS.FIT")) as hdul:
    ext = hdul["SCATS1"]
    valtime_str = ext.data['VALTIME']
    viewra_str = ext.data['VIEWRA']
    viewdecl_str = ext.data['VIEWDECL']

ra_ahf = np.array([Angle(s, unit=u.hourangle).deg for s in viewra_str])
dec_ahf = np.array([Angle(s, unit=u.deg).deg for s in viewdecl_str])

# Convert AHF ISO times to MJD for interpolation
t_ahf_time = Time(list(valtime_str), format='isot', scale='tt')
t_ahf_mjd = t_ahf_time.mjd

# Convert ATTHK mission times to MJD
xmm_epoch_mjd = 50814.0
t_atthk_mjd = xmm_epoch_mjd + t_atthk / 86400.0

# Filter to science-pointing AHF entries (within 0.5° of ATTHK median)
ra_mid = np.median(ra_atthk[ra_atthk < 300])  # exclude wraparound entries
dec_mid = np.median(dec_atthk[np.abs(dec_atthk) < 10])
science_mask = (np.abs(ra_ahf - ra_mid) < 0.5) & (np.abs(dec_ahf - dec_mid) < 0.5)

ra_ahf_sci = ra_ahf[science_mask]
dec_ahf_sci = dec_ahf[science_mask]
t_ahf_sci_mjd = t_ahf_mjd[science_mask]

# Interpolate atthkgen ATTHK onto AHF science times
ra_atthk_interp = np.interp(t_ahf_sci_mjd, t_atthk_mjd, ra_atthk)
dec_atthk_interp = np.interp(t_ahf_sci_mjd, t_atthk_mjd, dec_atthk)

dra = (ra_ahf_sci - ra_atthk_interp) * 3600 * np.cos(np.radians(dec_ahf_sci))
ddec = (dec_ahf_sci - dec_atthk_interp) * 3600
sep = np.sqrt(dra**2 + ddec**2)

print("=" * 60)
print("AHF vs atthkgen ATTHK (science-only, interpolated)")
print("=" * 60)
print(f"Samples: {len(ra_ahf_sci)}")
print(f"dRA*cos(dec): mean={dra.mean():.3f}\", median={np.median(dra):.3f}\", std={dra.std():.3f}\"")
print(f"dDec:         mean={ddec.mean():.3f}\", median={np.median(ddec):.3f}\", std={ddec.std():.3f}\"")
print(f"separation:   mean={sep.mean():.3f}\", median={np.median(sep):.3f}\", std={sep.std():.3f}\"")
print(f"separation:   min={sep.min():.3f}\", max={sep.max():.3f}\"")

# For comparison, also do AHF vs build_moved_atthk.py ATTHK
with fits.open(str(base / "comet" / "moved_atthk.dat")) as hdul:
    ext = hdul["ATTHK"]
    t_old = ext.data["TIME"].astype(np.float64)
    ra_old = ext.data["AHFRA"].astype(np.float64)
    dec_old = ext.data["AHFDEC"].astype(np.float64)

t_old_mjd = xmm_epoch_mjd + t_old / 86400.0
ra_old_interp = np.interp(t_ahf_sci_mjd, t_old_mjd, ra_old)
dec_old_interp = np.interp(t_ahf_sci_mjd, t_old_mjd, dec_old)

dra_old = (ra_ahf_sci - ra_old_interp) * 3600 * np.cos(np.radians(dec_ahf_sci))
ddec_old = (dec_ahf_sci - dec_old_interp) * 3600
sep_old = np.sqrt(dra_old**2 + ddec_old**2)

print()
print("=" * 60)
print("AHF vs build_moved_atthk.py ATTHK (for comparison)")
print("=" * 60)
print(f"Samples: {len(ra_ahf_sci)}")
print(f"dRA*cos(dec): mean={dra_old.mean():.3f}\", median={np.median(dra_old):.3f}\", std={dra_old.std():.3f}\"")
print(f"dDec:         mean={ddec_old.mean():.3f}\", median={np.median(ddec_old):.3f}\", std={ddec_old.std():.3f}\"")
print(f"separation:   mean={sep_old.mean():.3f}\", median={np.median(sep_old):.3f}\", std={sep_old.std():.3f}\"")
print(f"separation:   min={sep_old.min():.3f}\", max={sep_old.max():.3f}\"")

print()
print("=" * 60)
print("IMPROVEMENT RATIO")
print("=" * 60)
print(f"Median separation: atthkgen={np.median(sep):.3f}\" vs python={np.median(sep_old):.3f}\"  →  {np.median(sep_old)/np.median(sep):.1f}x better")
