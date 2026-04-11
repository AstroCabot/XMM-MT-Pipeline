#!/usr/bin/env python3
"""Part 3: Compare track time range to observation events and pointing."""
import numpy as np, math, glob, sys
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
sys.path.insert(0, "scripts")
from trail_mask_tools import read_track

track = read_track("output/track/comet_track.fits")

# Track time range
track_mjd_span = track.mjd[-1] - track.mjd[0]
print(f"Track MJD span: {track_mjd_span:.4f} days = {track_mjd_span*86400:.1f} s")
print(f"Track MJD: {track.mjd[0]:.6f} to {track.mjd[-1]:.6f}")

# Observation actual event coverage across ALL event files
all_tmin, all_tmax = [], []
evt_files = sorted(glob.glob("output/detect/stack_events*.fits"))
for ef in evt_files:
    with fits.open(ef, memmap=True) as hdul:
        t = hdul["EVENTS"].data["TIME"]
        all_tmin.append(t.min())
        all_tmax.append(t.max())

obs_tmin = min(all_tmin)
obs_tmax = max(all_tmax)
print(f"\nActual event time range: {obs_tmin:.1f} to {obs_tmax:.1f}")
print(f"Event time span: {obs_tmax - obs_tmin:.1f} s = {(obs_tmax-obs_tmin)/86400:.4f} days")

# Convert track MJDs to XMM time
# track.obs_t0 is the XMM time of the first track sample
track_xmm = track.obs_t0 + (track.mjd - track.mjd[0]) * 86400.0
print(f"\nTrack XMM time range: {track_xmm[0]:.1f} to {track_xmm[-1]:.1f}")
print(f"Track XMM time span: {track_xmm[-1] - track_xmm[0]:.1f} s")

print(f"\nTrack starts {track_xmm[0] - obs_tmin:.1f} s AFTER first event")
print(f"Track ends   {track_xmm[-1] - obs_tmax:.1f} s AFTER last event")

# Where in the track does the actual observation start/end?
obs_start_idx = np.searchsorted(track_xmm, obs_tmin)
obs_end_idx = np.searchsorted(track_xmm, obs_tmax)
print(f"\nFirst event at track index {obs_start_idx}/{len(track.ra)}")
print(f"Last event at track index {min(obs_end_idx, len(track.ra)-1)}/{len(track.ra)}")

if obs_start_idx < len(track.ra):
    print(f"Comet RA/Dec at first event: {track.ra[obs_start_idx]:.6f}, {track.dec[obs_start_idx]:.6f}")
if obs_end_idx < len(track.ra):
    ei = min(obs_end_idx, len(track.ra)-1)
    print(f"Comet RA/Dec at last event:  {track.ra[ei]:.6f}, {track.dec[ei]:.6f}")

# Flat-sky offsets for these
ctr_ra, ctr_dec = 177.94642639160156, 1.323562741279602
cos_dec = math.cos(math.radians(ctr_dec))
def flat_off(ra, dec):
    dx = ((ra - ctr_ra + 180.) % 360. - 180.) * cos_dec * 3600.
    dy = (dec - ctr_dec) * 3600.
    return dx, dy

if obs_start_idx < len(track.ra):
    dx0, dy0 = flat_off(track.ra[obs_start_idx], track.dec[obs_start_idx])
    print(f"Track at first-event offset: ({dx0:.1f}, {dy0:.1f}) arcsec")
ei = min(obs_end_idx, len(track.ra)-1)
dxe, dye = flat_off(track.ra[ei], track.dec[ei])
print(f"Track at last-event offset:  ({dxe:.1f}, {dye:.1f}) arcsec")

# Pointing
pnt_ra, pnt_dec = 177.547708333333, 1.532
nom_ra, nom_dec = 178.0888755, 1.3128889

def sep(ra1, dec1, ra2, dec2):
    c1 = SkyCoord(ra=ra1*u.deg, dec=dec1*u.deg)
    c2 = SkyCoord(ra=ra2*u.deg, dec=dec2*u.deg)
    return c1.separation(c2).arcsec

mid = len(track.ra) // 2
print(f"\n=== Pointing geometry ===")
print(f"RA_PNT, DEC_PNT = {pnt_ra:.6f}, {pnt_dec:.6f}")
print(f"RA_NOM, DEC_NOM = {nom_ra:.6f}, {nom_dec:.6f}")
print(f"Pointing to track start: {sep(pnt_ra, pnt_dec, track.ra[0], track.dec[0]):.1f} arcsec")
print(f"Pointing to track end:   {sep(pnt_ra, pnt_dec, track.ra[-1], track.dec[-1]):.1f} arcsec")
print(f"Pointing to track mid:   {sep(pnt_ra, pnt_dec, track.ra[mid], track.dec[mid]):.1f} arcsec")
print(f"RA_NOM to track start:   {sep(nom_ra, nom_dec, track.ra[0], track.dec[0]):.1f} arcsec")
print(f"RA_NOM to track end:     {sep(nom_ra, nom_dec, track.ra[-1], track.dec[-1]):.1f} arcsec")
print(f"RA_NOM to track mid:     {sep(nom_ra, nom_dec, track.ra[mid], track.dec[mid]):.1f} arcsec")
pnt_dx, pnt_dy = flat_off(pnt_ra, pnt_dec)
nom_dx, nom_dy = flat_off(nom_ra, nom_dec)
print(f"\nPointing offset from mosaic center: ({pnt_dx:.1f}, {pnt_dy:.1f}) arcsec")
print(f"RA_NOM offset from mosaic center:   ({nom_dx:.1f}, {nom_dy:.1f}) arcsec")
print(f"FOV radius: ~900 arcsec (15 arcmin)")

# Now check: the mosaic should show a circular-ish footprint centered on the pointing
# The data extends from x=-4080 to x=+1374. Let's see if that matches a circle at RA_PNT
print(f"\n=== Data footprint vs pointing circle ===")
print(f"Data x: [-4080, +1374], pointing at x={pnt_dx:.1f}")
print(f"Data y: [-918, +2070],  pointing at y={pnt_dy:.1f}")
print(f"Expected data range (pointing ± 900 arcsec):")
print(f"  x: [{pnt_dx-900:.1f}, {pnt_dx+900:.1f}]")
print(f"  y: [{pnt_dy-900:.1f}, {pnt_dy+900:.1f}]")
