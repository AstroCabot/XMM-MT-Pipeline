#!/usr/bin/env python3
"""Part 6: Check the attitude file. Was the spacecraft tracking the comet?"""
import numpy as np, math, sys
from astropy.io import fits
sys.path.insert(0, "scripts")
from trail_mask_tools import read_track

track = read_track("output/track/comet_track.fits")
track_xmm = track.obs_t0 + (track.mjd - track.mjd[0]) * 86400.0

# Read the attitude file
atthk_path = "output/repro/atthk.dat"
with fits.open(atthk_path) as hdul:
    print("Extensions:", [h.name for h in hdul])
    for ext in hdul:
        if hasattr(ext, 'columns') and ext.columns is not None:
            print(f"\n{ext.name} columns: {ext.columns.names}")
            if ext.data is not None and len(ext.data) > 0:
                print(f"  Rows: {len(ext.data)}")
                for col in ext.columns.names:
                    d = ext.data[col]
                    if np.issubdtype(d.dtype, np.floating):
                        print(f"  {col}: {d.min():.6f} to {d.max():.6f}")

# Load attitude data in one pass
with fits.open(atthk_path) as hdul:
    att_time = hdul[1].data["TIME"].astype(float).copy()
    att_ra = hdul[1].data["AHFRA"].astype(float).copy()
    att_dec = hdul[1].data["AHFDEC"].astype(float).copy()

# The main pointing columns
for ra_col, dec_col in [("AHFRA", "AHFDEC")]:
    try:
        print(f"\n=== Attitude ({ra_col}/{dec_col}) ===")
        print(f"Time: {att_time.min():.1f} to {att_time.max():.1f}")
        print(f"RA:   {att_ra.min():.6f} to {att_ra.max():.6f}")
        print(f"Dec:  {att_dec.min():.6f} to {att_dec.max():.6f}")
        print(f"Total RA sweep:  {(att_ra.max() - att_ra.min()) * 3600:.1f} arcsec")
        print(f"Total Dec sweep: {(att_dec.max() - att_dec.min()) * 3600:.1f} arcsec")
        
        # Check pointing at key times
        for label, t_target in [("first event", 881152770.0), 
                                ("track start", track_xmm[0]),
                                ("obs mid", 0.5*(att_time.min() + att_time.max())),
                                ("track mid", track_xmm[len(track_xmm)//2]),
                                ("last event", 881231472.0),
                                ("track end", track_xmm[-1])]:
            idx = np.argmin(np.abs(att_time - t_target))
            print(f"  At {label} (T={t_target:.0f}): "
                  f"pointing=({att_ra[idx]:.6f}, {att_dec[idx]:.6f})")
        
        # Compare attitude to comet position over time
        print(f"\n=== Pointing vs comet position over time ===")
        # Interpolate comet position to attitude timestamps
        comet_ra_interp = np.interp(att_time, track_xmm, track.ra, 
                                     left=np.nan, right=np.nan)
        comet_dec_interp = np.interp(att_time, track_xmm, track.dec,
                                      left=np.nan, right=np.nan)
        valid = np.isfinite(comet_ra_interp)
        sep_ra = (att_ra[valid] - comet_ra_interp[valid]) * np.cos(np.radians(att_dec[valid])) * 3600
        sep_dec = (att_dec[valid] - comet_dec_interp[valid]) * 3600
        sep_total = np.sqrt(sep_ra**2 + sep_dec**2)
        print(f"  Time overlap: {np.sum(valid)}/{len(att_time)} attitude samples")
        print(f"  RA offset (pointing - comet): {sep_ra.min():.1f} to {sep_ra.max():.1f} arcsec")
        print(f"  Dec offset (pointing - comet): {sep_dec.min():.1f} to {sep_dec.max():.1f} arcsec")
        print(f"  Total sep: min={sep_total.min():.1f}, max={sep_total.max():.1f}, "
              f"median={np.median(sep_total):.1f} arcsec")
        
        # Was the pointing TRACKING (moving with the comet)?
        att_ra_rate = np.diff(att_ra) / np.diff(att_time) * 3600  # arcsec/s
        comet_ra_rate_interp = np.interp(att_time[:-1], track_xmm, 
                                          np.gradient(track.ra, track_xmm - track_xmm[0]) * 3600,
                                          left=np.nan, right=np.nan)
        v = np.isfinite(comet_ra_rate_interp) & np.isfinite(att_ra_rate)
        if np.any(v):
            print(f"\n  Attitude RA rate: {np.nanmedian(att_ra_rate[v]):.4f} arcsec/s")
            print(f"  Comet RA rate:    {np.nanmedian(comet_ra_rate_interp[v]):.4f} arcsec/s")
            print(f"  Ratio (att/comet): {np.nanmedian(att_ra_rate[v] / comet_ra_rate_interp[v]):.3f}")
        
        break
    except KeyError:
        continue
