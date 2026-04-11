#!/usr/bin/env python3
"""Diagnose the projection mismatch in the detect QC plot."""
import sys, math
import numpy as np
sys.path.insert(0, "scripts")
from trail_mask_tools import read_track
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

# Read mosaic header + data
soft = fits.open("output/detect/EPIC_soft_mosaic.fits")
h = soft[0].header
data = soft[0].data
for k in ["CTRRA","CTRDEC","XMIN_AS","XMAX_AS","YMIN_AS","YMAX_AS","NAXIS1","NAXIS2"]:
    v = h.get(k)
    print(f"{k} = {v}")
print()

# Read the track
track = read_track("output/track/comet_track.fits")
print(f"Track samples: {len(track.ra)}")
print(f"Track RA range:  {track.ra.min():.6f} to {track.ra.max():.6f}")
print(f"Track Dec range: {track.dec.min():.6f} to {track.dec.max():.6f}")
mid = len(track.ra) // 2
ctr_ra = float(h["CTRRA"])
ctr_dec = float(h["CTRDEC"])
print(f"Mosaic center: RA={ctr_ra:.6f}, Dec={ctr_dec:.6f}")
print(f"Track midpoint: RA={track.ra[mid]:.6f}, Dec={track.dec[mid]:.6f}")
print(f"Track first:    RA={track.ra[0]:.6f}, Dec={track.dec[0]:.6f}")
print(f"Track last:     RA={track.ra[-1]:.6f}, Dec={track.dec[-1]:.6f}")
print()

# Method 1: flat-sky (what the mosaic uses)
def flat_offsets(ra, dec, cra, cdec):
    cos_dec = math.cos(math.radians(cdec))
    dx = ((np.asarray(ra) - cra + 180.) % 360. - 180.) * cos_dec * 3600.
    dy = (np.asarray(dec) - cdec) * 3600.
    return dx, dy

# Method 2: spherical_offsets_to (what checks.py uses)
def sph_offsets(ra, dec, cra, cdec):
    t = SkyCoord(ra=np.asarray(ra)*u.deg, dec=np.asarray(dec)*u.deg, frame="icrs")
    r = SkyCoord(ra=cra*u.deg, dec=cdec*u.deg, frame="icrs")
    dlon, dlat = r.spherical_offsets_to(t)
    return dlon.to_value(u.arcsec), dlat.to_value(u.arcsec)

fx, fy = flat_offsets(track.ra, track.dec, ctr_ra, ctr_dec)
sx, sy = sph_offsets(track.ra, track.dec, ctr_ra, ctr_dec)

print("=== Flat vs Spherical offsets for full track ===")
print(f"Flat  first: ({fx[0]:.2f}, {fy[0]:.2f})")
print(f"Spher first: ({sx[0]:.2f}, {sy[0]:.2f})")
print(f"Delta first: ({sx[0]-fx[0]:.2f}, {sy[0]-fy[0]:.2f})")
print()
print(f"Flat  last:  ({fx[-1]:.2f}, {fy[-1]:.2f})")
print(f"Spher last:  ({sx[-1]:.2f}, {sy[-1]:.2f})")
print(f"Delta last:  ({sx[-1]-fx[-1]:.2f}, {sy[-1]-fy[-1]:.2f})")
print()
print(f"Flat  mid:   ({fx[mid]:.2f}, {fy[mid]:.2f})")
print(f"Spher mid:   ({sx[mid]:.2f}, {sy[mid]:.2f})")
print(f"Delta mid:   ({sx[mid]-fx[mid]:.2f}, {sy[mid]-fy[mid]:.2f})")
print()
print(f"Max |dx|: {np.max(np.abs(sx-fx)):.2f} arcsec")
print(f"Max |dy|: {np.max(np.abs(sy-fy)):.2f} arcsec")
print(f"Max total offset: {np.max(np.sqrt((sx-fx)**2 + (sy-fy)**2)):.2f} arcsec")

# Check image extent
xmin = float(h["XMIN_AS"]); xmax = float(h["XMAX_AS"])
ymin = float(h["YMIN_AS"]); ymax = float(h["YMAX_AS"])
print(f"\nImage extent: x=[{xmin:.1f}, {xmax:.1f}], y=[{ymin:.1f}, {ymax:.1f}]")
print(f"Image size:  {xmax-xmin:.1f} x {ymax-ymin:.1f} arcsec")
print(f"Pixel scale: {(xmax-xmin)/data.shape[1]:.1f} arcsec/pixel")
print(f"Image shape: {data.shape} (ny, nx)")

# Where do the track endpoints fall in the mosaic data?
ny, nx = data.shape
bin_as = (xmax - xmin) / nx  # arcsec per pixel

# Convert flat track offsets to pixel indices
def offset_to_pixel(dx, dy):
    ix = ((dx - xmin) / (xmax - xmin) * nx).astype(int)
    iy = ((dy - ymin) / (ymax - ymin) * ny).astype(int)
    return ix, iy

pix_x, pix_y = offset_to_pixel(fx, fy)

print(f"\n=== Track endpoints in mosaic pixel coords ===")
print(f"Track first pixel: ({pix_x[0]}, {pix_y[0]})")
print(f"Track last pixel:  ({pix_x[-1]}, {pix_y[-1]})")
print(f"Track mid pixel:   ({pix_x[mid]}, {pix_y[mid]})")

# Check actual mosaic counts at track endpoints
for label, i in [("first", 0), ("mid", mid), ("last", -1)]:
    ix, iy = int(pix_x[i]), int(pix_y[i])
    if 0 <= ix < nx and 0 <= iy < ny:
        val = data[iy, ix]
        # Also check a small neighborhood
        r = 3
        iy0, iy1 = max(0,iy-r), min(ny,iy+r+1)
        ix0, ix1 = max(0,ix-r), min(nx,ix+r+1)
        neighborhood = data[iy0:iy1, ix0:ix1]
        print(f"Track {label}: pixel ({ix},{iy}), value={val:.1f}, "
              f"neighborhood max={neighborhood.max():.1f}, "
              f"neighborhood sum={neighborhood.sum():.1f}")
    else:
        print(f"Track {label}: pixel ({ix},{iy}) OUTSIDE IMAGE")

# Find the actual data footprint (where counts > 0)
nonzero = data > 0
rows_with_data = np.any(nonzero, axis=1)
cols_with_data = np.any(nonzero, axis=0)
if np.any(rows_with_data) and np.any(cols_with_data):
    row_range = np.where(rows_with_data)[0]
    col_range = np.where(cols_with_data)[0]
    # Convert to arcsec offsets
    data_xmin = xmin + col_range[0] * bin_as
    data_xmax = xmin + (col_range[-1] + 1) * bin_as
    data_ymin = ymin + row_range[0] * bin_as
    data_ymax = ymin + (row_range[-1] + 1) * bin_as
    print(f"\n=== Actual data footprint (nonzero pixels) ===")
    print(f"Data x range: [{data_xmin:.1f}, {data_xmax:.1f}] arcsec")
    print(f"Data y range: [{data_ymin:.1f}, {data_ymax:.1f}] arcsec")
    print(f"Track first at ({fx[0]:.1f}, {fy[0]:.1f}) - "
          f"{'INSIDE' if data_xmin <= fx[0] <= data_xmax and data_ymin <= fy[0] <= data_ymax else 'OUTSIDE'} data footprint")
    print(f"Track last  at ({fx[-1]:.1f}, {fy[-1]:.1f}) - "
          f"{'INSIDE' if data_xmin <= fx[-1] <= data_xmax and data_ymin <= fy[-1] <= data_ymax else 'OUTSIDE'} data footprint")

    # Check what fraction of the track is on nonzero mosaic pixels
    on_data = 0
    for i in range(len(fx)):
        ix, iy = int(pix_x[i]), int(pix_y[i])
        if 0 <= ix < nx and 0 <= iy < ny and data[iy, ix] > 0:
            on_data += 1
    print(f"\nTrack points on nonzero pixels: {on_data}/{len(fx)} ({100*on_data/len(fx):.1f}%)")

    # Check the nearby column at track first point: is there data above/below?
    print(f"\n=== Mosaic column at track start (x~{fx[0]:.0f} arcsec, col={pix_x[0]}) ===")
    col = data[:, min(pix_x[0], nx-1)]
    nz_rows = np.where(col > 0)[0]
    if len(nz_rows):
        print(f"Nonzero rows: {nz_rows[0]} to {nz_rows[-1]} "
              f"(y ~ {ymin + nz_rows[0]*bin_as:.0f} to {ymin + nz_rows[-1]*bin_as:.0f} arcsec)")
        print(f"Track start row: {pix_y[0]} (y ~ {fy[0]:.0f} arcsec)")
        if pix_y[0] < nz_rows[0]:
            print(f"Track start is {(nz_rows[0] - pix_y[0]) * bin_as:.0f} arcsec BELOW data in this column")
        elif pix_y[0] > nz_rows[-1]:
            print(f"Track start is {(pix_y[0] - nz_rows[-1]) * bin_as:.0f} arcsec ABOVE data in this column")
        else:
            print(f"Track start is within data in this column")
    else:
        print(f"No data in this column!")

# ============================================================
# Part 2: Track time range vs observation / event time range
# ============================================================
print("\n" + "="*60)
print("Part 2: Track vs observation time coverage")
print("="*60)

# Track time range
if hasattr(track, 'mjd') and len(track.mjd):
    print(f"\nTrack MJD range: {track.mjd[0]:.6f} to {track.mjd[-1]:.6f}")
    print(f"Track MJD span:  {track.mjd[-1] - track.mjd[0]:.6f} days = {(track.mjd[-1] - track.mjd[0])*86400:.1f} s")

if hasattr(track, 'obs_t0'):
    print(f"Track obs_t0:    {track.obs_t0}")

# Check track time field names
print(f"\nTrack attributes: {[a for a in dir(track) if not a.startswith('_')]}")

# Read an event file to get the observation time range
import glob
evt_files = sorted(glob.glob("output/detect/stack_events*.fits"))
if not evt_files:
    evt_files = sorted(glob.glob("output/clean/*.fits"))
print(f"\nFound {len(evt_files)} event files")

for ef in evt_files[:3]:
    with fits.open(ef, memmap=True) as hdul:
        for hdu_name in ['EVENTS', 'STDGTI', 'GTI']:
            if hdu_name in hdul:
                ext = hdul[hdu_name]
                if hdu_name == 'EVENTS':
                    t = ext.data['TIME']
                    print(f"\n{ef.split('/')[-1]} EVENTS time range:")
                    print(f"  TIME: {t.min():.1f} to {t.max():.1f} (span {t.max()-t.min():.1f} s)")
                    # Check for TSTART/TSTOP in header
                    ts = ext.header.get('TSTART')
                    te = ext.header.get('TSTOP')
                    if ts is not None:
                        print(f"  TSTART={ts}, TSTOP={te}")
                elif hdu_name in ('STDGTI', 'GTI'):
                    gti = ext.data
                    if gti is not None and len(gti):
                        starts = gti['START'] if 'START' in gti.columns.names else gti['START_TIME']
                        stops = gti['STOP'] if 'STOP' in gti.columns.names else gti['STOP_TIME']
                        print(f"  {hdu_name}: {len(gti)} intervals, "
                              f"first start={starts[0]:.1f}, last stop={stops[-1]:.1f}, "
                              f"total span={stops[-1]-starts[0]:.1f} s")

# Read the attitude file to find pointing center
att_env = "output/attitude.env"
import os
if os.path.exists(att_env):
    print(f"\n=== attitude.env ===")
    with open(att_env) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                print(f"  {line}")

# Check the spacecraft pointing from event header
if evt_files:
    with fits.open(evt_files[0], memmap=True) as hdul:
        h0 = hdul[0].header
        for k in ['RA_PNT', 'DEC_PNT', 'RA_NOM', 'DEC_NOM', 'RA_OBJ', 'DEC_OBJ',
                   'TSTART', 'TSTOP', 'DATE-OBS', 'DATE-END', 'OBJECT']:
            v = h0.get(k)
            if v is not None:
                print(f"  Primary header {k} = {v}")
        if 'EVENTS' in hdul:
            he = hdul['EVENTS'].header
            for k in ['RA_PNT', 'DEC_PNT', 'RA_NOM', 'DEC_NOM', 'TSTART', 'TSTOP']:
                v = he.get(k)
                if v is not None:
                    print(f"  EVENTS header {k} = {v}")

# Compare track time to event times
# First check what time the track uses
print(f"\n=== Track time info ===")
if hasattr(track, 'time') and len(track.time):
    print(f"track.time range: {track.time[0]:.1f} to {track.time[-1]:.1f}")
    print(f"track.time span:  {track.time[-1]-track.time[0]:.1f} s")
elif hasattr(track, 'obs_t0'):
    # Time might be stored as offset from obs_t0
    pass

# Check the source of the track
track_path = "output/track/comet_track.fits"
with fits.open(track_path) as hdul:
    print(f"\nTrack FITS extensions: {[h.name for h in hdul]}")
    for ext in hdul:
        if hasattr(ext, 'columns') and ext.columns is not None:
            print(f"  {ext.name} columns: {ext.columns.names}")
            if ext.data is not None:
                for col in ext.columns.names:
                    d = ext.data[col]
                    if np.issubdtype(d.dtype, np.floating):
                        print(f"    {col}: {d.min():.6f} to {d.max():.6f}")
        h = ext.header
        for k in ['TSTART', 'TSTOP', 'OBS_T0', 'TIMEZERO']:
            v = h.get(k)
            if v is not None:
                print(f"  {ext.name} header {k} = {v}")
