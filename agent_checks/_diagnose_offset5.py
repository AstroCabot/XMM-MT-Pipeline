#!/usr/bin/env python3
"""Part 5: Check the merged clean events — what WCS do they use?
Compare their RA/Dec extent to the mosaic and track."""
import numpy as np, math, glob, sys
from astropy.io import fits
from astropy.wcs import WCS
sys.path.insert(0, "scripts")
from trail_mask_tools import read_track

track = read_track("output/track/comet_track.fits")
ctr_ra, ctr_dec = 177.94642639160156, 1.323562741279602
cos_dec = math.cos(math.radians(ctr_dec))

def flat_off(ra, dec):
    dx = ((np.asarray(ra) - ctr_ra + 180.) % 360. - 180.) * cos_dec * 3600.
    dy = (np.asarray(dec) - ctr_dec) * 3600.
    return dx, dy

# Check the actual merged cleaned event files
for inst in ["PN", "M1", "M2"]:
    path = f"output/clean/{inst}_clean_merged.fits"
    try:
        with fits.open(path, memmap=True) as hdul:
            evt = hdul["EVENTS"]
            h = evt.header
            cols = [c.upper() for c in evt.columns.names]
            print(f"\n=== {inst}: {path} ===")
            print(f"Columns: {evt.columns.names}")
            print(f"Nevents: {len(evt.data)}")
            
            t = evt.data["TIME"]
            print(f"TIME: {t.min():.1f} to {t.max():.1f} ({t.max()-t.min():.0f} s)")
            
            # Get WCS info for X/Y columns
            xcol = ycol = None
            for i in range(1, h.get("TFIELDS", 0) + 1):
                name = h.get(f"TTYPE{i}", "").strip().upper()
                if name == "X": xcol = i
                elif name == "Y": ycol = i
            
            if xcol and ycol:
                crvl_x = h.get(f"TCRVL{xcol}")
                crvl_y = h.get(f"TCRVL{ycol}")
                print(f"TCRVL (WCS ref): ({crvl_x:.6f}, {crvl_y:.6f})")
                
                w = WCS(naxis=2)
                w.wcs.crpix = [h[f"TCRPX{xcol}"], h[f"TCRPX{ycol}"]]
                w.wcs.cdelt = [h[f"TCDLT{xcol}"], h[f"TCDLT{ycol}"]]
                w.wcs.crval = [crvl_x, crvl_y]
                w.wcs.ctype = [h.get(f"TCTYP{xcol}", "RA---TAN"),
                              h.get(f"TCTYP{ycol}", "DEC--TAN")]
                
                x = evt.data["X"].astype(float)
                y = evt.data["Y"].astype(float)
                ra, dec = w.all_pix2world(x, y, 1)
                dx, dy = flat_off(ra, dec)
                
                print(f"X range: {x.min():.0f} to {x.max():.0f}")
                print(f"Y range: {y.min():.0f} to {y.max():.0f}")
                print(f"RA range:  {np.nanmin(ra):.6f} to {np.nanmax(ra):.6f}")
                print(f"Dec range: {np.nanmin(dec):.6f} to {np.nanmax(dec):.6f}")
                print(f"dx range: {np.nanmin(dx):.1f} to {np.nanmax(dx):.1f} arcsec")
                print(f"dy range: {np.nanmin(dy):.1f} to {np.nanmax(dy):.1f} arcsec")
                
            if "RA" in cols and "DEC" in cols:
                ra_col = evt.columns.names[cols.index("RA")]
                dec_col = evt.columns.names[cols.index("DEC")]
                ra_d = evt.data[ra_col].astype(float)
                dec_d = evt.data[dec_col].astype(float)
                dx_d, dy_d = flat_off(ra_d, dec_d)
                print(f"Direct RA range:  {ra_d.min():.6f} to {ra_d.max():.6f}")
                print(f"Direct Dec range: {dec_d.min():.6f} to {dec_d.max():.6f}")
                print(f"Direct dx range: {dx_d.min():.1f} to {dx_d.max():.1f} arcsec")
                print(f"Direct dy range: {dy_d.min():.1f} to {dy_d.max():.1f} arcsec")
                
    except Exception as e:
        print(f"{inst}: {e}")

# Track endpoints for reference
print(f"\n=== Track reference ===")
fx, fy = flat_off(track.ra, track.dec)
print(f"Track start offset: ({fx[0]:.1f}, {fy[0]:.1f}) arcsec")
print(f"Track end offset:   ({fx[-1]:.1f}, {fy[-1]:.1f}) arcsec")
print(f"Track mid offset:   ({fx[len(fx)//2]:.1f}, {fy[len(fy)//2]:.1f}) arcsec")
