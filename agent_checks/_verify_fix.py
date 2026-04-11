#!/usr/bin/env python3
"""Verify the fix: re-query Horizons with quantities='1,20' and compare."""
import sys, math
import numpy as np
from astropy.io import fits
from astropy.time import Time
sys.path.insert(0, "scripts")
from trail_mask_tools import read_track
from astroquery.jplhorizons import Horizons

# Load current (apparent) track
track = read_track("output/track/comet_track.fits")

# Sample epochs
sample_idx = [0, len(track.mjd)//2, -1]
sample_jd = [float(track.mjd[i] + 2400000.5) for i in sample_idx]

# Query with the FIXED quantity (1 = astrometric/ICRS)
obj = Horizons(id="C/2025 N1", id_type="smallbody",
               location="500@399", epochs=sample_jd)
eph = obj.ephemerides(quantities="1,20", extra_precision=True)

print("Fixed query returns columns:", eph.colnames)
ra_fixed = np.array(eph["RA"], float)
dec_fixed = np.array(eph["DEC"], float)

# Check event RA/DEC from clean event files for comparison
from pathlib import Path
from astropy.wcs import WCS
evt_files = sorted(Path("output/clean").glob("*_clean.fits"))

for ef in evt_files[:1]:
    with fits.open(ef, memmap=True) as hdul:
        if "EVENTS" not in hdul:
            continue
        evt = hdul["EVENTS"]
        data = evt.data
        hdr = evt.header
        
        # Get coordinate system info
        radecsys = hdr.get("RADECSYS", hdr.get("RADESYS", "unknown"))
        print(f"\nEvent file: {ef.name}")
        print(f"  Coordinate frame: {radecsys}")
        
        # Get RA/DEC of events
        pi = np.asarray(data["PI"], float)
        soft = (pi >= 200) & (pi <= 1000)
        
        if "RA" in evt.columns.names:
            ra_evt = np.asarray(data["RA"][soft], float)
            dec_evt = np.asarray(data["DEC"][soft], float)
        else:
            # Use X/Y WCS
            xcol = ycol = None
            for i in range(1, hdr.get("TFIELDS", 0) + 1):
                name = hdr.get(f"TTYPE{i}", "").strip().upper()
                if name == "X": xcol = i
                elif name == "Y": ycol = i
            w = WCS(naxis=2)
            w.wcs.crpix = [hdr[f"TCRPX{xcol}"], hdr[f"TCRPX{ycol}"]]
            w.wcs.cdelt = [hdr[f"TCDLT{xcol}"], hdr[f"TCDLT{ycol}"]]
            w.wcs.crval = [hdr[f"TCRVL{xcol}"], hdr[f"TCRVL{ycol}"]]
            w.wcs.ctype = [hdr.get(f"TCTYP{xcol}", "RA---TAN"), 
                          hdr.get(f"TCTYP{ycol}", "DEC--TAN")]
            x = np.asarray(data["X"][soft], float)
            y = np.asarray(data["Y"][soft], float)
            ra_evt, dec_evt = w.all_pix2world(x, y, 1)
        
        print(f"  Soft events: {len(ra_evt)}")
        print(f"  Event RA:  {np.median(ra_evt):.6f} (median)")
        print(f"  Event DEC: {np.median(dec_evt):.6f} (median)")

print("\n=== FIXED track (astrometric/ICRS) vs event coordinates ===")
for i, si in enumerate(sample_idx):
    label = ["start", "mid", "end"][i]
    dra = (ra_fixed[i] - np.median(ra_evt)) * math.cos(math.radians(dec_fixed[i])) * 3600
    ddec = (dec_fixed[i] - np.median(dec_evt)) * 3600
    print(f"  Track {label}: RA={ra_fixed[i]:.6f} DEC={dec_fixed[i]:.6f}")

print("\n=== OLD track (apparent) vs FIXED track (astrometric) ===")
for i, si in enumerate(sample_idx):
    label = ["start", "mid", "end"][i]
    old_ra = track.ra[si]
    old_dec = track.dec[si]
    dra = (old_ra - ra_fixed[i]) * math.cos(math.radians(old_dec)) * 3600
    ddec = (old_dec - dec_fixed[i]) * 3600
    sep = math.sqrt(dra**2 + ddec**2)
    print(f"  {label}: old=({old_ra:.6f},{old_dec:.6f}) fixed=({ra_fixed[i]:.6f},{dec_fixed[i]:.6f})")
    print(f"     offset: dRA*cos(d)={dra:+.1f}\" dDEC={ddec:+.1f}\" sep={sep:.1f}\"")
