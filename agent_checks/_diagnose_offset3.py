#!/usr/bin/env python3
"""Part 3: Detailed sub-exposure analysis — where does each land on the mosaic?"""
import numpy as np, math, glob, sys, os
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
sys.path.insert(0, "scripts")
from trail_mask_tools import read_track

track = read_track("output/track/comet_track.fits")
ctr_ra, ctr_dec = 177.94642639160156, 1.323562741279602
cos_dec = math.cos(math.radians(ctr_dec))

def flat_off(ra, dec):
    dx = ((np.asarray(ra) - ctr_ra + 180.) % 360. - 180.) * cos_dec * 3600.
    dy = (np.asarray(dec) - ctr_dec) * 3600.
    return dx, dy

# Check ALL event files in detect, with per-file stats  
evt_files = sorted(glob.glob("output/detect/stack_events*.fits"))
print(f"{'File':<45} {'Tmin':>14} {'Tmax':>14} {'Span_s':>8} "
      f"{'Nevt':>8} {'Xmin':>8} {'Xmax':>8} {'Ymin':>8} {'Ymax':>8}")
print("-" * 155)

all_tmin, all_tmax = [], []
for ef in evt_files:
    name = os.path.basename(ef)
    with fits.open(ef, memmap=True) as hdul:
        evt = hdul["EVENTS"]
        t = evt.data["TIME"]
        
        # Get RA/Dec
        cols_upper = [c.upper() for c in evt.columns.names]
        if "RA" in cols_upper and "DEC" in cols_upper:
            ra_col = evt.columns.names[cols_upper.index("RA")]
            dec_col = evt.columns.names[cols_upper.index("DEC")]
            ra = np.asarray(evt.data[ra_col], dtype=float)
            dec = np.asarray(evt.data[dec_col], dtype=float)
        else:
            ra = dec = None
            
        all_tmin.append(t.min())
        all_tmax.append(t.max())
        
        if ra is not None:
            dx, dy = flat_off(ra, dec)
            print(f"{name:<45} {t.min():>14.1f} {t.max():>14.1f} {t.max()-t.min():>8.0f} "
                  f"{len(t):>8d} {dx.min():>8.1f} {dx.max():>8.1f} {dy.min():>8.1f} {dy.max():>8.1f}")
        else:
            print(f"{name:<45} {t.min():>14.1f} {t.max():>14.1f} {t.max()-t.min():>8.0f} "
                  f"{len(t):>8d} {'no RA':>8} {'':>8} {'':>8} {'':>8}")

obs_tmin = min(all_tmin)
obs_tmax = max(all_tmax)
print(f"\nOverall event time range: {obs_tmin:.1f} to {obs_tmax:.1f} ({obs_tmax-obs_tmin:.0f} s)")

# Track time
track_xmm = track.obs_t0 + (track.mjd - track.mjd[0]) * 86400.0
print(f"Track XMM time range:    {track_xmm[0]:.1f} to {track_xmm[-1]:.1f} ({track_xmm[-1]-track_xmm[0]:.0f} s)")

# What's the comet position at the last event time?
last_evt_idx = np.searchsorted(track_xmm, obs_tmax)
last_evt_idx = min(last_evt_idx, len(track.ra)-1)
dx_last, dy_last = flat_off(track.ra[last_evt_idx], track.dec[last_evt_idx])
print(f"\nComet at last event time (idx {last_evt_idx}):")
print(f"  RA={track.ra[last_evt_idx]:.6f} Dec={track.dec[last_evt_idx]:.6f}")
print(f"  Mosaic offset: ({dx_last:.1f}, {dy_last:.1f}) arcsec")

# What's the comet position at track start (first sample)?
dx_first, dy_first = flat_off(track.ra[0], track.dec[0])
print(f"\nComet at track start (idx 0):")
print(f"  RA={track.ra[0]:.6f} Dec={track.dec[0]:.6f}")
print(f"  Mosaic offset: ({dx_first:.1f}, {dy_first:.1f}) arcsec")

# What about the first event's time in the track?
first_evt_idx = np.searchsorted(track_xmm, obs_tmin)
first_evt_idx = min(first_evt_idx, len(track.ra)-1)
dx_fe, dy_fe = flat_off(track.ra[first_evt_idx], track.dec[first_evt_idx])
print(f"\nComet at first event time (idx {first_evt_idx}):")
print(f"  RA={track.ra[first_evt_idx]:.6f} Dec={track.dec[first_evt_idx]:.6f}")
print(f"  Mosaic offset: ({dx_fe:.1f}, {dy_fe:.1f}) arcsec")

# Now read the soft mosaic and check the actual data profile along the track direction
soft = fits.open("output/detect/EPIC_soft_mosaic.fits")
data = soft[0].data
h = soft[0].header
xmin_as, xmax_as = float(h["XMIN_AS"]), float(h["XMAX_AS"])
ymin_as, ymax_as = float(h["YMIN_AS"]), float(h["YMAX_AS"])
ny, nx = data.shape
bin_as = (xmax_as - xmin_as) / nx

# Check the data sum along columns (x-axis profile)  
col_sums = data.sum(axis=0)
nonzero_cols = np.where(col_sums > 0)[0]
print(f"\n=== Mosaic x profile ===")
print(f"Nonzero columns: {nonzero_cols[0]} to {nonzero_cols[-1]}")
print(f"  = x offsets {xmin_as + nonzero_cols[0]*bin_as:.1f} to {xmin_as + (nonzero_cols[-1]+1)*bin_as:.1f} arcsec")

# Check at x=1862 (track start)
target_col = int((dx_first - xmin_as) / bin_as)
print(f"\nTrack start (x={dx_first:.0f}) is column {target_col}")
if 0 <= target_col < nx:
    print(f"Column sum at track start: {col_sums[target_col]:.0f}")
# Check nearby columns
for offset in [-100, -50, -20, -10, -5, 0, 5]:
    col = target_col + offset
    if 0 <= col < nx:
        x_as = xmin_as + col * bin_as
        print(f"  col {col} (x={x_as:.0f}): sum={col_sums[col]:.0f}")

# Also check: RA_PNT for each sub-exposure from their headers
print(f"\n=== Per-sub-exposure pointing ===")
for ef in evt_files:
    name = os.path.basename(ef)
    with fits.open(ef, memmap=True) as hdul:
        h0 = hdul[0].header
        he = hdul["EVENTS"].header
        ra_pnt = he.get("RA_PNT", h0.get("RA_PNT"))
        dec_pnt = he.get("DEC_PNT", h0.get("DEC_PNT"))
        if ra_pnt is not None:
            pnt_dx, pnt_dy = flat_off(ra_pnt, dec_pnt)
            t = hdul["EVENTS"].data["TIME"]
            print(f"  {name:<45} PNT=({pnt_dx:>8.1f}, {pnt_dy:>8.1f})  "
                  f"T=[{t.min():.0f}, {t.max():.0f}]")
