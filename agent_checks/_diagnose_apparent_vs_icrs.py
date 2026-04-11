#!/usr/bin/env python3
"""Diagnose the comet track offset: apparent vs ICRS coordinates.

The Horizons query uses quantities="2,20" which returns apparent RA/DEC
(aberration + precession/nutation to true equator/equinox of date).
XMM event RA/DEC are in ICRS. This script quantifies the mismatch.
"""
import sys, math
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, "scripts")
from trail_mask_tools import read_track

# ── 1. Load the stored comet track ──────────────────────────────────────
track = read_track("output/track/comet_track.fits")
print(f"Track has {len(track.mjd)} points")
print(f"Track MJD range: {track.mjd[0]:.6f} to {track.mjd[-1]:.6f}")
print(f"Track RA range:  {track.ra[0]:.6f} to {track.ra[-1]:.6f}")
print(f"Track DEC range: {track.dec[0]:.6f} to {track.dec[-1]:.6f}")
print(f"MJDREF={track.mjdref}, TIMESYS={track.timesys}")

# ── 2. Re-query Horizons with BOTH astrometric and apparent ────────────
from astroquery.jplhorizons import Horizons

# Use a few sample epochs from the track
sample_idx = [0, len(track.mjd)//4, len(track.mjd)//2, 3*len(track.mjd)//4, -1]
sample_jd = [float(track.mjd[i] + 2400000.5) for i in sample_idx]

print(f"\nSample JD epochs: {sample_jd}")

# Astrometric (quantity 1) = ICRF, no aberration
obj_astro = Horizons(id="C/2025 N1", id_type="smallbody",
                     location="500@399", epochs=sample_jd)
eph_astro = obj_astro.ephemerides(quantities="1,20", extra_precision=True)

print("\n=== Astrometric (quantity 1) columns ===")
print(eph_astro.colnames)

# Apparent (quantity 2) = true equator/equinox of date, with aberration
obj_app = Horizons(id="C/2025 N1", id_type="smallbody",
                   location="500@399", epochs=sample_jd)
eph_app = obj_app.ephemerides(quantities="2,20", extra_precision=True)

print("\n=== Apparent (quantity 2) columns ===")
print(eph_app.colnames)

# Extract coordinates
# Astrometric: RA_ICRF, DEC_ICRF (or similar)
ra_col_a = [c for c in eph_astro.colnames if 'RA' in c.upper() and 'ICRF' in c.upper()]
dec_col_a = [c for c in eph_astro.colnames if 'DEC' in c.upper() and 'ICRF' in c.upper()]
print(f"\nAstrometric RA col: {ra_col_a}")
print(f"Astrometric DEC col: {dec_col_a}")

# Apparent: RA_app, DEC_app (or similar)
ra_col_app = [c for c in eph_app.colnames if 'RA' in c.upper() and 'app' in c.lower()]
dec_col_app = [c for c in eph_app.colnames if 'DEC' in c.upper() and 'app' in c.lower()]
print(f"Apparent RA col: {ra_col_app}")
print(f"Apparent DEC col: {dec_col_app}")

# Get the actual values
if ra_col_a:
    ra_astro = np.array(eph_astro[ra_col_a[0]], float)
    dec_astro = np.array(eph_astro[dec_col_a[0]], float)
else:
    # Try just "RA"
    ra_astro = np.array(eph_astro["RA"], float)
    dec_astro = np.array(eph_astro["DEC"], float)

if ra_col_app:
    ra_app = np.array(eph_app[ra_col_app[0]], float)
    dec_app = np.array(eph_app[dec_col_app[0]], float)
else:
    ra_app = np.array(eph_app["RA"], float)
    dec_app = np.array(eph_app["DEC"], float)

print("\n=== Comparison: Apparent vs Astrometric ===")
print(f"{'Epoch JD':>16s} {'RA_astro':>12s} {'RA_app':>12s} {'DEC_astro':>12s} {'DEC_app':>12s}  {'dRA*cos(d)':>12s}  {'dDEC':>12s}  {'sep':>8s}")
print(" " * 70 + "(arcsec)    (arcsec)   (arcsec)")

for i in range(len(sample_jd)):
    dra = (ra_app[i] - ra_astro[i]) * math.cos(math.radians(dec_astro[i])) * 3600.0
    ddec = (dec_app[i] - dec_astro[i]) * 3600.0
    sep = math.sqrt(dra**2 + ddec**2)
    print(f"{sample_jd[i]:16.6f} {ra_astro[i]:12.6f} {ra_app[i]:12.6f} {dec_astro[i]:12.6f} {dec_app[i]:12.6f}  {dra:+12.2f}  {ddec:+12.2f}  {sep:8.2f}")

# ── 3. Compare stored track to astrometric values ──────────────────────
print("\n=== Stored track vs astrometric at sample points ===")
for i, si in enumerate(sample_idx):
    stored_ra = track.ra[si]
    stored_dec = track.dec[si]
    dra = (stored_ra - ra_astro[i]) * math.cos(math.radians(dec_astro[i])) * 3600.0
    ddec = (stored_dec - dec_astro[i]) * 3600.0
    sep_stored_astro = math.sqrt(dra**2 + ddec**2)
    
    dra2 = (stored_ra - ra_app[i]) * math.cos(math.radians(dec_app[i])) * 3600.0
    ddec2 = (stored_dec - dec_app[i]) * 3600.0
    sep_stored_app = math.sqrt(dra2**2 + ddec2**2)
    
    print(f"  Point {i}: stored RA={stored_ra:.6f} DEC={stored_dec:.6f}")
    print(f"    vs astrometric: dRA*cos(d)={dra:+.2f}\", dDEC={ddec:+.2f}\", sep={sep_stored_astro:.2f}\"")
    print(f"    vs apparent:    dRA*cos(d)={dra2:+.2f}\", dDEC={ddec2:+.2f}\", sep={sep_stored_app:.2f}\"")

# ── 4. Also check what the mosaic center is ────────────────────────────
print("\n=== Detect mosaic properties ===")
for band in ("soft", "hard"):
    mpath = Path("output/detect") / f"EPIC_{band}_mosaic.fits"
    if mpath.exists():
        with fits.open(mpath) as hdul:
            hdr = hdul[0].header
            print(f"\n{band} mosaic:")
            print(f"  CTRRA  = {hdr.get('CTRRA', 'missing')}")
            print(f"  CTRDEC = {hdr.get('CTRDEC', 'missing')}")
            print(f"  XMIN_AS = {hdr.get('XMIN_AS')}, XMAX_AS = {hdr.get('XMAX_AS')}")
            print(f"  YMIN_AS = {hdr.get('YMIN_AS')}, YMAX_AS = {hdr.get('YMAX_AS')}")

# ── 5. Check event file RA/DEC to confirm ICRS ────────────────────────
print("\n=== Event file coordinate frame check ===")
evt_files = list(Path("output/clean").glob("*_clean.fits"))
if not evt_files:
    evt_files = list(Path("output/repro").glob("*.fits"))
for ef in evt_files[:1]:
    with fits.open(ef, memmap=True) as hdul:
        if "EVENTS" in hdul:
            hdr = hdul["EVENTS"].header
            print(f"\n{ef.name}:")
            for kw in ("RADECSYS", "RADESYS", "EQUINOX", "REFXCTYP", "REFYCTYP",
                        "TCTYP11", "TCTYP12", "TCTYP8", "TCTYP9"):
                v = hdr.get(kw)
                if v is not None:
                    print(f"  {kw} = {v}")
            # Check X/Y WCS columns
            for i in range(1, hdr.get("TFIELDS", 0) + 1):
                name = hdr.get(f"TTYPE{i}", "").strip().upper()
                if name in ("X", "Y", "RA", "DEC"):
                    ctype = hdr.get(f"TCTYP{i}", "")
                    cunit = hdr.get(f"TCUNI{i}", "")
                    print(f"  Column {i} ({name}): CTYPE={ctype}, UNIT={cunit}")
