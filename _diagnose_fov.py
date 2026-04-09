"""Diagnose track offset from FOV strip + source distribution asymmetry."""
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
import csv, math

# ── Read data ──
ref_ra, ref_dec = 177.9464263916, 1.3235627413  # mosaic center = track midpoint

with fits.open("output/track/comet_track.fits") as hdul:
    tab = hdul[1].data
    trk_ra = np.asarray(tab["RA"], dtype=float)
    trk_dec = np.asarray(tab["DEC"], dtype=float)
n_track = len(trk_ra)
mid = n_track // 2

# Original telescope pointing (RA_NOM from event header = commanded target position)
# RA_PNT=177.548 is from the MOVED ATTHK processed file, not the original pointing
ra_nom, dec_nom = 178.0888755, 1.3128889  # from event FITS headers
print(f"RA_NOM = {ra_nom:.7f}  (commanded target = approx original boresight)")
print(f"DEC_NOM = {dec_nom:.7f}")
print(f"RA_PNT = 177.548  (moved-frame boresight at reference time)")
print(f"Ref attitude = ({ref_ra:.6f}, {ref_dec:.6f})")
print()

# Verify: RA_PNT ≈ ref + (RA_NOM - comet_start)?
pnt_check = ref_ra + (ra_nom - trk_ra[0])
print(f"Check: ref + (RA_NOM - comet_start) = {ref_ra:.4f} + ({ra_nom:.4f} - {trk_ra[0]:.4f}) = {pnt_check:.4f}")
print(f"  vs RA_PNT = 177.548 → diff = {pnt_check - 177.548:.4f}° (should be small)")
print()

cos_dec = math.cos(math.radians(ref_dec))

# ── FOV center path in comet frame ──
# In comet frame: FOV_center(t) = ref + (boresight_original - comet(t))
# ≈ ref + (RA_NOM - comet(t))
# Expressed as arcsec offset from mosaic center:
# dE(t) = (RA_NOM - comet_ra(t)) * cos(dec) * 3600
# dN(t) = (DEC_NOM - comet_dec(t)) * 3600

fov_dE = (ra_nom - trk_ra) * cos_dec * 3600
fov_dN = (dec_nom - trk_dec) * 3600

# Also compute track offsets
trk_dE = (trk_ra - ref_ra) * cos_dec * 3600
trk_dN = (trk_dec - ref_dec) * 3600

print("=" * 60)
print("FOV CENTER PATH in comet frame (offset from mosaic center):")
print(f"  At START (t=0):  East={fov_dE[0]:+.0f}\", North={fov_dN[0]:+.0f}\"")
print(f"  At MID:          East={fov_dE[mid]:+.0f}\", North={fov_dN[mid]:+.0f}\"")
print(f"  At END:          East={fov_dE[-1]:+.0f}\", North={fov_dN[-1]:+.0f}\"")
print()
print("COMET TRACK (sky coords, offset from mosaic center):")
print(f"  At START (t=0):  East={trk_dE[0]:+.0f}\", North={trk_dN[0]:+.0f}\"")
print(f"  At MID:          East={trk_dE[mid]:+.0f}\", North={trk_dN[mid]:+.0f}\"")
print(f"  At END:          East={trk_dE[-1]:+.0f}\", North={trk_dN[-1]:+.0f}\"")
print()

# Key systematic offset
fov_mid_E = fov_dE[mid]
fov_mid_N = fov_dN[mid]
print(f"SYSTEMATIC OFFSET (FOV center - track center at midpoint):")
print(f"  ΔEast = {fov_mid_E:+.1f}\" = {fov_mid_E/60:+.1f}'")
print(f"  ΔNorth = {fov_mid_N:+.1f}\" = {fov_mid_N/60:+.1f}'")
total_off = math.sqrt(fov_mid_E**2 + fov_mid_N**2)
print(f"  Total = {total_off:.1f}\" = {total_off/60:.1f}'")
print()

# ── When is the comet inside the FOV? ──
fov_radius = 900  # ~15' for EPIC pn
# In comet frame, comet is at (0,0). FOV center is at (fov_dE, fov_dN).
# Comet is in FOV when distance(fov_center, origin) < fov_radius
dist_to_comet = np.sqrt(fov_dE**2 + fov_dN**2)
in_fov = dist_to_comet < fov_radius

# Find first and last time comet is in FOV
entries = np.where(in_fov)[0]
if len(entries) > 0:
    t_enter_frac = entries[0] / n_track
    t_exit_frac = entries[-1] / n_track
    print(f"Comet in FOV (pn ~15' radius):")
    print(f"  Enters at track sample {entries[0]}/{n_track} ({t_enter_frac*100:.1f}% of observation)")
    print(f"  Exits at track sample {entries[-1]}/{n_track} ({t_exit_frac*100:.1f}% of observation)")
    print(f"  In-FOV fraction: {len(entries)/n_track*100:.1f}%")
    print(f"  At entry: FOV center=({fov_dE[entries[0]]:+.0f}, {fov_dN[entries[0]]:+.0f}), dist={dist_to_comet[entries[0]]:.0f}\"")
    print(f"  At exit:  FOV center=({fov_dE[entries[-1]]:+.0f}, {fov_dN[entries[-1]]:+.0f}), dist={dist_to_comet[entries[-1]]:.0f}\"")
else:
    print("Comet NEVER in FOV! Something is wrong.")
print()

# ── Where does the FOV strip cover the track start/end positions? ──
# Track start at (trk_dE[0], trk_dN[0]) = (+1862, -663)
# When does the FOV circle pass over this point?
dist_fov_to_trkstart = np.sqrt((fov_dE - trk_dE[0])**2 + (fov_dN - trk_dN[0])**2)
covers_start = np.where(dist_fov_to_trkstart < fov_radius)[0]
if len(covers_start) > 0:
    print(f"FOV covers track START position ({trk_dE[0]:+.0f}, {trk_dN[0]:+.0f}):")
    print(f"  From sample {covers_start[0]} to {covers_start[-1]} ({covers_start[0]/n_track*100:.1f}% to {covers_start[-1]/n_track*100:.1f}%)")
    print(f"  = {len(covers_start)/n_track*100:.1f}% of observation time")
    print(f"  This is at the END of the observation (FOV has swept to lower-right)")
else:
    print(f"FOV NEVER covers track start position! Track start is always outside FOV.")

dist_fov_to_trkend = np.sqrt((fov_dE - trk_dE[-1])**2 + (fov_dN - trk_dN[-2])**2)
covers_end = np.where(dist_fov_to_trkend < fov_radius)[0]
if len(covers_end) > 0:
    print(f"FOV covers track END position ({trk_dE[-1]:+.0f}, {trk_dN[-1]:+.0f}):")
    print(f"  From sample {covers_end[0]} to {covers_end[-1]} ({covers_end[0]/n_track*100:.1f}% to {covers_end[-1]/n_track*100:.1f}%)")
    print(f"  = {len(covers_end)/n_track*100:.1f}% of observation time")
else:
    print(f"FOV NEVER covers track end position! Track end is always outside FOV.")
print()

# ── PN event distribution vs FOV coverage ──
print("=" * 60)
print("PN EVENT DISTRIBUTION (from earlier analysis):")
print("  Bins 0-3 (first 40% of obs): ~12% of events → HEAVY FLARING")
print("  Bins 4-9 (last 60% of obs):  ~88% of events → CLEAN")
print()
print("This means:")
print("  - Upper-left FOV strip (early time): minimal effective exposure")
print("  - Lower-right FOV strip (late time): deep, clean exposure")
print("  - Source detection SENSITIVITY is much higher in the lower-right region")

# ── Source distribution vs track ──
try:
    with fits.open("output/detect/field_sources_all.fits") as hdul:
        for hdu in hdul[1:]:
            if hasattr(hdu, 'data') and hdu.data is not None and len(hdu.data) > 0:
                src = hdu.data; break
        src_ra = np.asarray(src['RA'], dtype=float)
        src_dec = np.asarray(src['DEC'], dtype=float)
        dx_src = (src_ra - ref_ra) * cos_dec * 3600
        dy_src = (src_dec - ref_dec) * 3600

    print(f"\n\nSOURCE POSITIONS vs FOV STRIP:")
    print(f"  Sources: {len(src_ra)} total")
    
    # FOV strip runs from (fov_dE[0],fov_dN[0]) to (fov_dE[-1],fov_dN[-1])
    # = from ({fov_dE[0]:+.0f},{fov_dN[0]:+.0f}) to ({fov_dE[-1]:+.0f},{fov_dN[-1]:+.0f})
    # Track runs from (trk_dE[0],trk_dN[0]) to (trk_dE[-1],trk_dN[-1])
    
    # Count sources in early-time FOV region vs late-time FOV region
    # Split at the midpoint of FOV path
    # FOV moves from upper-left to lower-right
    # early time: East < fov_mid_E, late time: East > fov_mid_E
    n_early = np.sum(dx_src < 0)  # upper-left (negative East)
    n_late = np.sum(dx_src >= 0)   # lower-right (positive East)
    print(f"  Sources in upper-left (East<0, early-time FOV region): {n_early}")
    print(f"  Sources in lower-right (East≥0, late-time FOV region): {n_late}")
    
    # But also check where sources are relative to the FOV strip centerline
    # (not just left/right of the mosaic center)
    print(f"\n  Sources by East offset bins:")
    bins = [(-3000, -2000), (-2000, -1000), (-1000, 0), (0, 1000), (1000, 2000), (2000, 3000)]
    for lo, hi in bins:
        n = np.sum((dx_src >= lo) & (dx_src < hi))
        print(f"    [{lo:+5d} to {hi:+5d}\"]: {n:3d} sources")

except Exception as e:
    print(f"Error reading sources: {e}")
