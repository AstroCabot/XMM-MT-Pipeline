#!/usr/bin/env python3
"""Final verification: astrometric track vs event coordinates & old apparent track."""
import sys, math
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
sys.path.insert(0, "scripts")
from trail_mask_tools import read_track
from astroquery.jplhorizons import Horizons

# Load current (apparent) track
track = read_track("output/track/comet_track.fits")

# Sample at start/mid/end
idx = [0, len(track.mjd)//2, -1]
jds = [float(track.mjd[i] + 2400000.5) for i in idx]

# Query ICRS (astrometric)
obj = Horizons(id="C/2025 N1", id_type="smallbody",
               location="500@399", epochs=jds)
eph = obj.ephemerides(quantities="1,20", extra_precision=True)
ra_icrs = np.array(eph["RA"], float)
dec_icrs = np.array(eph["DEC"], float)

# Read event median coordinates from merged clean file
evt_path = "output/clean/PN_clean_merged.fits"
with fits.open(evt_path, memmap=True) as hdul:
    evt = hdul["EVENTS"]
    hdr = evt.header
    data = evt.data
    pi = np.asarray(data["PI"], float)
    soft = (pi >= 200) & (pi <= 1000)
    
    # Check what columns exist
    cols = [c.upper() for c in evt.columns.names]
    if "RA" in cols:
        ra_col = evt.columns.names[cols.index("RA")]
        dec_col = evt.columns.names[cols.index("DEC")]
        ra_evt = np.asarray(data[ra_col][soft], float)
        dec_evt = np.asarray(data[dec_col][soft], float)
    else:
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
    
    print(f"Event file: {evt_path}")
    print(f"  RADECSYS: {hdr.get('RADECSYS', hdr.get('RADESYS', 'N/A'))}")
    print(f"  Soft events: {np.sum(soft)}")
    print(f"  RA median:  {np.median(ra_evt):.6f}")
    print(f"  DEC median: {np.median(dec_evt):.6f}")

print("\n" + "="*70)
print(f"{'':8s} {'Old (apparent)':>30s}  {'Fixed (ICRS)':>30s}")
print(f"{'':8s} {'RA':>12s} {'DEC':>12s} {'sep':>6s}  {'RA':>12s} {'DEC':>12s} {'sep':>6s}")
print("="*70)

med_ra = np.median(ra_evt)
med_dec = np.median(dec_evt)

for i, si in enumerate(idx):
    label = ["start", "mid  ", "end  "][i]
    
    # Old (apparent) vs events
    dra_old = (track.ra[si] - med_ra) * math.cos(math.radians(med_dec)) * 3600
    ddec_old = (track.dec[si] - med_dec) * 3600
    sep_old = math.sqrt(dra_old**2 + ddec_old**2)
    
    # Fixed (ICRS) vs events
    dra_new = (ra_icrs[i] - med_ra) * math.cos(math.radians(med_dec)) * 3600
    ddec_new = (dec_icrs[i] - med_dec) * 3600
    sep_new = math.sqrt(dra_new**2 + ddec_new**2)
    
    print(f"  {label}: {track.ra[si]:12.6f} {track.dec[si]:12.6f} {sep_old:6.0f}\"  "
          f"{ra_icrs[i]:12.6f} {dec_icrs[i]:12.6f} {sep_new:6.0f}\"")

print("\n=== Summary ===")
print(f"Old track is ~1302\" from event data (WRONG - apparent vs ICRS mismatch)")
print(f"Fixed track should be within a few arcsec of event data (CORRECT - both ICRS)")
