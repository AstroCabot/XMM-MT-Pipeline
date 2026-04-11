#!/usr/bin/env python3
"""Part 4: How does the mosaic align photons? Check WCS in event files."""
import numpy as np, math, glob, sys
from astropy.io import fits
from astropy.wcs import WCS
sys.path.insert(0, "scripts")

# Pick one event file and dump its WCS-relevant headers
ef = "output/detect/stack_events0963720201E11PN.fits"
with fits.open(ef, memmap=True) as hdul:
    evt = hdul["EVENTS"]
    h = evt.header
    
    # What columns exist?
    print("=== EVENTS columns ===")
    print(evt.columns.names)
    
    # Find X/Y column indices
    print("\n=== WCS-related header keywords ===")
    for i in range(1, h.get("TFIELDS", 0) + 1):
        name = h.get(f"TTYPE{i}", "").strip()
        if name.upper() in ("X", "Y", "RA", "DEC", "TIME", "RAWX", "RAWY"):
            print(f"\nColumn {i}: {name}")
            for prefix in ["TCRPX", "TCDLT", "TCRVL", "TCTYP", "TCUNI", "TLMIN", "TLMAX"]:
                key = f"{prefix}{i}"
                val = h.get(key)
                if val is not None:
                    print(f"  {key} = {val}")

    # Key header values
    print("\n=== Key header values ===")
    for k in ["RA_PNT", "DEC_PNT", "RA_NOM", "DEC_NOM", "REFXCRVL", "REFYCRVL",
              "REFXCRPX", "REFYCRPX", "REFXCDLT", "REFYCDLT",
              "PA_PNT", "TSTART", "TSTOP", "INSTRUME"]:
        v = h.get(k)
        if v is not None:
            print(f"  {k} = {v}")

    # Show actual X/Y data range
    x = evt.data["X"].astype(float)
    y = evt.data["Y"].astype(float)
    print(f"\n=== X/Y data range ===")
    print(f"X: {x.min():.1f} to {x.max():.1f}")
    print(f"Y: {y.min():.1f} to {y.max():.1f}")
    
    # Build the WCS and convert corners to RA/Dec
    xcol = ycol = None
    for i in range(1, h.get("TFIELDS", 0) + 1):
        name = h.get(f"TTYPE{i}", "").strip().upper()
        if name == "X": xcol = i
        elif name == "Y": ycol = i
    
    w = WCS(naxis=2)
    w.wcs.crpix = [h[f"TCRPX{xcol}"], h[f"TCRPX{ycol}"]]
    w.wcs.cdelt = [h[f"TCDLT{xcol}"], h[f"TCDLT{ycol}"]]
    w.wcs.crval = [h[f"TCRVL{xcol}"], h[f"TCRVL{ycol}"]]
    w.wcs.ctype = [h.get(f"TCTYP{xcol}", "RA---TAN"), h.get(f"TCTYP{ycol}", "DEC--TAN")]
    
    print(f"\n=== Constructed WCS ===")
    print(f"CRPIX = {w.wcs.crpix}")
    print(f"CDELT = {w.wcs.cdelt}")
    print(f"CRVAL = {w.wcs.crval}")
    print(f"CTYPE = {list(w.wcs.ctype)}")
    
    # Convert event X/Y to RA/Dec
    ra, dec = w.all_pix2world(x, y, 1)
    print(f"\n=== Converted RA/Dec range ===")
    print(f"RA:  {np.nanmin(ra):.6f} to {np.nanmax(ra):.6f}")
    print(f"Dec: {np.nanmin(dec):.6f} to {np.nanmax(dec):.6f}")
    
    # Compare to mosaic center
    ctr_ra, ctr_dec = 177.94642639160156, 1.323562741279602
    cos_dec = math.cos(math.radians(ctr_dec))
    dx = ((ra - ctr_ra + 180.) % 360. - 180.) * cos_dec * 3600.
    dy = (dec - ctr_dec) * 3600.
    print(f"\n=== Event offsets from mosaic center ===")
    print(f"dx: {np.nanmin(dx):.1f} to {np.nanmax(dx):.1f} arcsec")
    print(f"dy: {np.nanmin(dy):.1f} to {np.nanmax(dy):.1f} arcsec")

# Now check: is TCRVL the same across all sub-exposures?
print(f"\n=== TCRVL (WCS reference) across all sub-exposures ===")
evt_files = sorted(glob.glob("output/detect/stack_events*.fits"))
for ef in evt_files:
    name = ef.split("/")[-1]
    with fits.open(ef, memmap=True) as hdul:
        h = hdul["EVENTS"].header
        xcol = ycol = None
        for i in range(1, h.get("TFIELDS", 0) + 1):
            nm = h.get(f"TTYPE{i}", "").strip().upper()
            if nm == "X": xcol = i
            elif nm == "Y": ycol = i
        if xcol and ycol:
            crvl_x = h.get(f"TCRVL{xcol}")
            crvl_y = h.get(f"TCRVL{ycol}")
            t = hdul["EVENTS"].data["TIME"]
            print(f"  {name:<45} TCRVL=({crvl_x:.6f}, {crvl_y:.6f})  "
                  f"T=[{t.min():.0f}, {t.max():.0f}]")

# Check: how was this event file created? Look at the pipeline
print("\n=== How does the pipeline create these detect event files? ===")
print("Checking pipeline script for estackrun / stage_detect...")
