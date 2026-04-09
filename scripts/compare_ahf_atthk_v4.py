#!/usr/bin/env python3
"""Compare attmove AHF (sexagesimal RA/DEC, ISO time) with moved ATTHK."""
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.coordinates import SkyCoord, Angle
from astropy.time import Time
import astropy.units as u

base = Path("/mnt/c/Users/cabot/OneDrive/Documents/Research/Cambridge/X3I/reduction_v4/output")

# ── 1. Read AHF and convert ──
ahf_path = base / "comet" / "moved_odf" / "odf" / "4759_0963720201_SCX00000ATS.FIT"
print("Reading AHF...")
with fits.open(str(ahf_path)) as hdul:
    ext = hdul["SCATS1"]
    valtime_str = ext.data['VALTIME']
    viewra_str = ext.data['VIEWRA']
    viewdecl_str = ext.data['VIEWDECL']

# Convert sexagesimal to degrees
# VIEWRA is in HH:MM:SS.ss (hours), VIEWDECL is in DD:MM:SS.s
print("Converting sexagesimal to degrees...")
ra_ahf_deg = np.array([Angle(s, unit=u.hourangle).deg for s in viewra_str])
dec_ahf_deg = np.array([Angle(s, unit=u.deg).deg for s in viewdecl_str])

# Convert ISO time to XMM mission time
# XMM epoch: 1998-01-01T00:00:00 TT
print("Converting timestamps...")
t_ahf_astropy = Time(list(valtime_str), format='isot', scale='utc')
xmm_epoch = Time('1998-01-01T00:00:00', format='isot', scale='tt')
# Approximate: mission time = seconds since epoch
t_ahf_mjd = t_ahf_astropy.mjd
# Actually just use the Time objects for comparison

print(f"AHF: {len(ra_ahf_deg)} rows")
print(f"Time range: {valtime_str[0]} to {valtime_str[-1]}")
print(f"VIEWRA range: [{ra_ahf_deg.min():.6f}, {ra_ahf_deg.max():.6f}] deg")
print(f"VIEWDECL range: [{dec_ahf_deg.min():.6f}, {dec_ahf_deg.max():.6f}] deg")

# ── 2. Read ATTHK ──
with fits.open(str(base / "comet" / "moved_atthk.dat")) as hdul:
    ext = hdul["ATTHK"]
    t_atthk = ext.data["TIME"].astype(np.float64)
    ra_atthk = ext.data["AHFRA"].astype(np.float64)
    dec_atthk = ext.data["AHFDEC"].astype(np.float64)

print(f"\nATTHK: {len(t_atthk)} rows")
print(f"AHFRA range: [{ra_atthk.min():.6f}, {ra_atthk.max():.6f}] deg")
print(f"AHFDEC range: [{dec_atthk.min():.6f}, {dec_atthk.max():.6f}] deg")

# ── 3. Convert ATTHK times to ISO for matching ──
# XMM mission time epoch: 1998-01-01T00:00:00 TT in MJD = 50814.0
# Mission time in seconds since this epoch
xmm_epoch_mjd = 50814.0  # 1998-01-01 in MJD (TT)
t_atthk_mjd = xmm_epoch_mjd + t_atthk / 86400.0
t_atthk_time = Time(t_atthk_mjd, format='mjd', scale='tt')
print(f"ATTHK time range: {t_atthk_time[0].isot} to {t_atthk_time[-1].isot}")

# ── 4. Find AHF entries within the ATTHK time range ──
t_min = t_atthk_time[0]
t_max = t_atthk_time[-1]
mask = (t_ahf_astropy >= t_min) & (t_ahf_astropy <= t_max)
n_overlap = mask.sum()
print(f"\n{'='*60}")
print(f"AHF entries within ATTHK time range: {n_overlap} of {len(ra_ahf_deg)}")
print(f"{'='*60}")

if n_overlap > 0:
    ra_ahf_overlap = ra_ahf_deg[mask]
    dec_ahf_overlap = dec_ahf_deg[mask]
    t_ahf_overlap = t_ahf_astropy[mask]
    
    print(f"AHF (in overlap):  RA range=[{ra_ahf_overlap.min():.6f}, {ra_ahf_overlap.max():.6f}]")
    print(f"                   DEC range=[{dec_ahf_overlap.min():.6f}, {dec_ahf_overlap.max():.6f}]")
    print(f"AHF median:  RA={np.median(ra_ahf_overlap):.6f}, DEC={np.median(dec_ahf_overlap):.6f}")
    print(f"ATTHK median: RA={np.median(ra_atthk):.6f}, DEC={np.median(dec_atthk):.6f}")
    
    # Compute median offset
    dra = (np.median(ra_ahf_overlap) - np.median(ra_atthk)) * 3600 * np.cos(np.radians(np.median(dec_ahf_overlap)))
    ddec = (np.median(dec_ahf_overlap) - np.median(dec_atthk)) * 3600
    print(f"\nMedian pointing offset (AHF - ATTHK):")
    print(f"  dRA*cos(dec) = {dra:.2f}\"")
    print(f"  dDec         = {ddec:.2f}\"")
    print(f"  total        = {np.sqrt(dra**2 + ddec**2):.2f}\"")
    
    # ── 5. Interpolate ATTHK to AHF overlap times and compare per-sample ──
    # Convert AHF overlap times to mission seconds for interpolation
    t_ahf_overlap_mjd = t_ahf_overlap.mjd
    t_atthk_mjd_arr = t_atthk_time.mjd
    
    ra_atthk_interp = np.interp(t_ahf_overlap_mjd, t_atthk_mjd_arr, ra_atthk)
    dec_atthk_interp = np.interp(t_ahf_overlap_mjd, t_atthk_mjd_arr, dec_atthk)
    
    dra_arr = (ra_ahf_overlap - ra_atthk_interp) * 3600 * np.cos(np.radians(dec_ahf_overlap))
    ddec_arr = (dec_ahf_overlap - dec_atthk_interp) * 3600
    sep_arr = np.sqrt(dra_arr**2 + ddec_arr**2)
    
    print(f"\nPer-sample comparison ({n_overlap} samples):")
    print(f"  dRA*cos(dec): mean={dra_arr.mean():.2f}\", median={np.median(dra_arr):.2f}\", std={dra_arr.std():.2f}\"")
    print(f"  dDec:         mean={ddec_arr.mean():.2f}\", median={np.median(ddec_arr):.2f}\", std={ddec_arr.std():.2f}\"")
    print(f"  separation:   mean={sep_arr.mean():.2f}\", median={np.median(sep_arr):.2f}\", std={sep_arr.std():.2f}\"")
    print(f"  separation:   min={sep_arr.min():.2f}\", max={sep_arr.max():.2f}\"")
    
    # Print some samples  
    idx = np.linspace(0, n_overlap-1, min(15, n_overlap), dtype=int)
    print(f"\n{'IDX':>5s} {'TIME':>23s} {'dRA_cos':>10s} {'dDec':>10s} {'sep':>10s}")
    for i in idx:
        print(f"{i:5d} {t_ahf_overlap[i].isot:>23s} {dra_arr[i]:10.2f}\" {ddec_arr[i]:10.2f}\" {sep_arr[i]:10.2f}\"")
    
    # Also check: filter to only "science" pointing (exclude slews)
    # During science, RA should change slowly; during slews, it changes fast
    if n_overlap > 10:
        # Look at where the pointing is near the ATTHK range (comet observation)
        ra_mid = np.median(ra_atthk)
        dec_mid = np.median(dec_atthk)
        science_mask = (np.abs(ra_ahf_overlap - ra_mid) < 0.5) & (np.abs(dec_ahf_overlap - dec_mid) < 0.5)
        n_sci = science_mask.sum()
        print(f"\n{'='*60}")
        print(f"SCIENCE-ONLY (within 0.5° of ATTHK median): {n_sci} of {n_overlap} AHF samples")
        print(f"{'='*60}")
        if n_sci > 0:
            dra_sci = dra_arr[science_mask]
            ddec_sci = ddec_arr[science_mask]
            sep_sci = sep_arr[science_mask]
            print(f"  dRA*cos(dec): mean={dra_sci.mean():.2f}\", median={np.median(dra_sci):.2f}\", std={dra_sci.std():.2f}\"")
            print(f"  dDec:         mean={ddec_sci.mean():.2f}\", median={np.median(ddec_sci):.2f}\", std={ddec_sci.std():.2f}\"")
            print(f"  separation:   mean={sep_sci.mean():.2f}\", median={np.median(sep_sci):.2f}\", std={sep_sci.std():.2f}\"")
            print(f"  separation:   min={sep_sci.min():.2f}\", max={sep_sci.max():.2f}\"")
else:
    print("NO OVERLAP — time ranges don't match!")
    print("This means the AHF and ATTHK cover different time periods")
    print("AHF may need TT/UTC conversion to match")
    
    # Try with TT scale for AHF
    print("\nRetrying with TT scale for AHF timestamps...")
    t_ahf_tt = Time(list(valtime_str), format='isot', scale='tt')
    mask_tt = (t_ahf_tt >= t_min) & (t_ahf_tt <= t_max)
    print(f"With TT: {mask_tt.sum()} AHF rows overlap with ATTHK time range")
