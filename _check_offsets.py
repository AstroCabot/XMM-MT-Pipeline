#!/usr/bin/env python3
"""Quick diagnostic: compare track positions with mosaic center."""
import numpy as np
import csv
import math
from astropy.io import fits

# --- Read motion segments ---
segs = []
with open('output/track/motion_segments.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        segs.append(row)

print(f"Number of 30-s segments: {len(segs)}")
print(f"First segment UTC: {segs[0]['utc_start']}")
print(f"Last segment UTC:  {segs[-1]['utc_stop']}")

# --- Read track FITS ---
with fits.open('output/track/comet_track.fits') as hdul:
    for hdu in hdul[1:]:
        if hasattr(hdu, 'data') and hdu.data is not None and len(hdu.data) > 0:
            tab = hdu.data
            break
    track_ra = np.asarray(tab['RA'], dtype=float)
    track_dec = np.asarray(tab['DEC'], dtype=float)
    track_t = np.asarray(tab['TIME_XMM'], dtype=float) if 'TIME_XMM' in tab.dtype.names else None
    print(f"\nTrack FITS: {len(track_ra)} samples")
    print(f"  Columns: {tab.dtype.names}")
    print(f"  RA range:  {track_ra[0]:.6f} -> {track_ra[-1]:.6f}")
    print(f"  Dec range: {track_dec[0]:.6f} -> {track_dec[-1]:.6f}")
    mid = len(track_ra) // 2
    print(f"  Midpoint [index {mid}]: RA={track_ra[mid]:.6f}, Dec={track_dec[mid]:.6f}")
    if track_t is not None:
        print(f"  Time range: {track_t[0]:.3f} -> {track_t[-1]:.3f}")
        print(f"  Duration: {track_t[-1]-track_t[0]:.1f} s = {(track_t[-1]-track_t[0])/3600:.2f} hr")
        dt = np.diff(track_t)
        print(f"  Step: median={np.median(dt):.3f} s, min={np.min(dt):.3f}, max={np.max(dt):.3f}")

# --- Mosaic center ---
ctr_ra = 177.9464263916
ctr_dec = 1.3235627413
print(f"\nMosaic center (ref_attitude.env): RA={ctr_ra:.10f}, Dec={ctr_dec:.10f}")

# --- Compute offsets ---
cos_dec = math.cos(math.radians(ctr_dec))

def offsets(ra, dec):
    dx = (ra - ctr_ra) * cos_dec * 3600
    dy = (dec - ctr_dec) * 3600
    return dx, dy

dx_track, dy_track = offsets(track_ra, track_dec)
print(f"\nTrack offsets from mosaic center (arcsec):")
print(f"  Start:  dE={dx_track[0]:+.1f}, dN={dy_track[0]:+.1f}")
print(f"  Mid:    dE={dx_track[mid]:+.1f}, dN={dy_track[mid]:+.1f}")
print(f"  End:    dE={dx_track[-1]:+.1f}, dN={dy_track[-1]:+.1f}")
print(f"  E range: {dx_track.min():.1f} to {dx_track.max():.1f}")
print(f"  N range: {dy_track.min():.1f} to {dy_track.max():.1f}")

# --- Check segment RA/Dec span ---
seg_ra_first = float(segs[0]['ra_mid_deg'])
seg_dec_first = float(segs[0]['dec_mid_deg'])
seg_ra_last = float(segs[-1]['ra_mid_deg'])
seg_dec_last = float(segs[-1]['dec_mid_deg'])
dx_sf, dy_sf = offsets(seg_ra_first, seg_dec_first)
dx_sl, dy_sl = offsets(seg_ra_last, seg_dec_last)
print(f"\nSegment midpoints offset from mosaic center:")
print(f"  First: dE={dx_sf:+.1f}, dN={dy_sf:+.1f}")
print(f"  Last:  dE={dx_sl:+.1f}, dN={dy_sl:+.1f}")

# --- Read source list ---
try:
    with fits.open('output/detect/field_sources_all.fits') as hdul:
        for hdu in hdul[1:]:
            if hasattr(hdu, 'data') and hdu.data is not None and len(hdu.data) > 0:
                src = hdu.data
                break
        src_ra = np.asarray(src['RA'], dtype=float)
        src_dec = np.asarray(src['DEC'], dtype=float)
        dx_src, dy_src = offsets(src_ra, src_dec)
        print(f"\nAll detected sources: {len(src_ra)}")
        print(f"  E range: {dx_src.min():.1f} to {dx_src.max():.1f}")
        print(f"  N range: {dy_src.min():.1f} to {dy_src.max():.1f}")
        # How many sources left vs right of center
        n_left = np.sum(dx_src < 0)
        n_right = np.sum(dx_src >= 0)
        print(f"  Sources left (E<0): {n_left}")
        print(f"  Sources right (E>=0): {n_right}")
        # How many sources left vs right of track midpoint
        dx_track_mid = dx_track[mid]
        n_left_trk = np.sum(dx_src < dx_track_mid)
        n_right_trk = np.sum(dx_src >= dx_track_mid)
        print(f"  Sources left of track mid ({dx_track_mid:.0f}\"): {n_left_trk}")
        print(f"  Sources right of track mid: {n_right_trk}")
except Exception as e:
    print(f"\nCould not read source list: {e}")

# --- Check how choose_center works ---
# The make_still_sky_mosaic uses choose_center(), which if no --center-ra/dec
# is given but --track is provided, uses track midpoint:
#   mid = len(track.ra) // 2; return (track.ra[mid], track.dec[mid])
# But the pipeline calls it WITHOUT --center-ra/--center-dec args
# Let's check what the pipeline actually passes
print(f"\n--- KEY FINDING ---")
print(f"ref_attitude.env center: RA={ctr_ra:.10f}, Dec={ctr_dec:.10f}")
print(f"Track midpoint [{mid}]:    RA={track_ra[mid]:.10f}, Dec={track_dec[mid]:.10f}")
sep_ra = (track_ra[mid] - ctr_ra) * cos_dec * 3600
sep_dec = (track_dec[mid] - ctr_dec) * 3600
print(f"Difference: dRA*cos(dec)={sep_ra:.2f}\", dDec={sep_dec:.2f}\"")
print(f"These should be ~identical if the mosaic uses the track midpoint as center")

# --- Check GTI / exposure patterns ---
print(f"\n\n=== GTI / EXPOSURE ANALYSIS ===")
for inst in ["PN", "M1", "M2"]:
    for fname in [f"output/clean/{inst}_clean_merged.fits", f"output/clean/{inst}_clean.fits"]:
        try:
            with fits.open(fname) as hdul:
                if "EVENTS" in hdul:
                    evt = hdul["EVENTS"].data
                    times = np.asarray(evt["TIME"], dtype=float)
                    print(f"\n{inst} ({fname}):")
                    print(f"  Events: {len(times)}")
                    print(f"  Time: {times.min():.3f} to {times.max():.3f}")
                    print(f"  Duration: {(times.max() - times.min())/3600:.2f} hr")
                    
                    # Check for big gaps (soft proton flaring removal)
                    sorted_t = np.sort(times)
                    dt = np.diff(sorted_t)
                    big_gaps = np.where(dt > 300)[0]  # gaps > 5 min
                    if len(big_gaps) > 0:
                        print(f"  *** {len(big_gaps)} gaps > 300s detected:")
                        for gi in big_gaps[:10]:
                            gap_start = sorted_t[gi]
                            gap_dur = dt[gi]
                            print(f"      t={gap_start:.0f}, dur={gap_dur:.0f}s = {gap_dur/3600:.2f}hr")
                    
                    # Time histogram to see exposure distribution
                    t0 = float(segs[0]["tstart"])
                    t1 = float(segs[-1]["tstop"])
                    bins = np.linspace(t0, t1, 11)  # 10 time bins
                    hist, _ = np.histogram(times, bins=bins)
                    bin_dur = (t1 - t0) / 10
                    print(f"  Event counts in 10 equal time bins ({bin_dur/3600:.2f} hr each):")
                    for i, h in enumerate(hist):
                        pct = h / len(times) * 100
                        bar = "#" * int(pct)
                        print(f"    bin {i}: {h:7d} events ({pct:5.1f}%) {bar}")

                # Check for GTI extension
                for ext_name in hdul:
                    if hasattr(ext_name, 'name') and 'GTI' in str(ext_name.name).upper():
                        gti = ext_name.data
                        if gti is not None and len(gti) > 0:
                            cols = gti.dtype.names
                            start_col = [c for c in cols if 'START' in c.upper()][0] if any('START' in c.upper() for c in cols) else cols[0]
                            stop_col = [c for c in cols if 'STOP' in c.upper()][0] if any('STOP' in c.upper() for c in cols) else cols[1]
                            total_gti = sum(gti[stop_col] - gti[start_col])
                            print(f"  GTI ({ext_name.name}): {len(gti)} intervals, total={total_gti:.0f}s = {total_gti/3600:.2f}hr")
                            if len(gti) <= 20:
                                for j, row in enumerate(gti):
                                    print(f"    [{j}] {row[start_col]:.0f} -- {row[stop_col]:.0f} ({(row[stop_col]-row[start_col])/3600:.2f}hr)")
                break
        except FileNotFoundError:
            continue

# --- Check if FOV pattern explains source asymmetry ---
print(f"\n\n=== SOURCE DISTRIBUTION vs EXPOSURE ===")
try:
    with fits.open('output/detect/field_sources_all.fits') as hdul:
        for hdu in hdul[1:]:
            if hasattr(hdu, 'data') and hdu.data is not None and len(hdu.data) > 0:
                src = hdu.data
                break
        src_ra = np.asarray(src['RA'], dtype=float)
        src_dec = np.asarray(src['DEC'], dtype=float)
        dx_src = (src_ra - ctr_ra) * cos_dec * 3600
        dy_src = (src_dec - ctr_dec) * 3600
        
        # Divide FOV into quadrants relative to track
        # Track runs roughly East-West (dx ~ +1862 to -1880)
        # "Right side of trail" in image = EAST = positive dx (early obs)
        # "Left side of trail" = WEST = negative dx (late obs)
        
        # Count sources in East vs West halves
        # But relative to where the FOV has exposure
        n_east_1000 = np.sum(dx_src > 1000)
        n_east_500 = np.sum((dx_src > 500) & (dx_src <= 1000))
        n_center = np.sum((dx_src >= -500) & (dx_src <= 500))
        n_west_500 = np.sum((dx_src < -500) & (dx_src >= -1000))
        n_west_1000 = np.sum(dx_src < -1000)
        
        print(f"Source distribution along E-W axis:")
        print(f"  East > 1000\":  {n_east_1000}")
        print(f"  500-1000\":     {n_east_500}")
        print(f"  -500 to 500\":  {n_center}")
        print(f"  -1000 to -500\": {n_west_500}")
        print(f"  West < -1000\": {n_west_1000}")
except Exception as e:
    print(f"Error: {e}")

# --- Check mosaic FITS headers for WCS details ---
print(f"\n\n=== MOSAIC FITS HEADERS ===")
import glob
for pattern in ["output/detect_mosaic/*.fits", "output/detect/*.fits"]:
    files = glob.glob(pattern)
    for f in sorted(files):
        if "mosaic" in f.lower():
            hdr = fits.getheader(f)
            print(f"\n{f}:")
            for key in ["NAXIS1", "NAXIS2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", 
                        "CDELT1", "CDELT2", "CTYPE1", "CTYPE2", 
                        "CTRRA", "CTRDEC", "XMIN_AS", "XMAX_AS", "YMIN_AS", "YMAX_AS"]:
                if key in hdr:
                    print(f"  {key} = {hdr[key]}")

# --- Check telescope pointing from event headers ---
print(f"\n\n=== TELESCOPE POINTING ===")
for inst in ["PN", "M1", "M2"]:
    for fname in [f"output/clean/{inst}_clean_merged.fits", f"output/clean/{inst}_clean.fits"]:
        try:
            hdr = fits.getheader(fname, "EVENTS")
            pnt_ra = hdr.get("RA_PNT", hdr.get("RA_NOM"))
            pnt_dec = hdr.get("DEC_PNT", hdr.get("DEC_NOM"))
            pa = hdr.get("PA_PNT", hdr.get("PA_NOM"))
            print(f"{inst}: RA_PNT={pnt_ra}, DEC_PNT={pnt_dec}, PA_PNT={pa}")
            for k in sorted(hdr.keys()):
                if any(x in k.upper() for x in ["RA_", "DEC_", "PA_", "TSTART", "TSTOP", "EXPID", "OBS_ID"]):
                    print(f"  {k} = {hdr[k]}")
            break
        except:
            continue

# --- Compute FOV center position in comet frame ---
print(f"\n\n=== FOV GEOMETRY IN COMET FRAME ===")
# Read original attitude pointing from first event file
try:
    with fits.open("output/clean/PN_clean_merged.fits") as hdul:
        hdr = hdul["EVENTS"].header
        pnt_ra = float(hdr.get("RA_PNT", hdr.get("RA_NOM", 0)))
        pnt_dec = float(hdr.get("DEC_PNT", hdr.get("DEC_NOM", 0)))
    
    print(f"Telescope pointing (approx fixed): RA={pnt_ra:.6f}, Dec={pnt_dec:.6f}")
    
    # In comet frame, FOV center at time t = ref + (pnt - comet_t)
    # At track start: comet at (178.462, 1.140)
    # At track end: comet at (177.424, 1.510)
    
    comet_start_ra, comet_start_dec = float(trk_ra[0]), float(trk_dec[0])
    comet_end_ra, comet_end_dec = float(trk_ra[-1]), float(trk_dec[-1])
    
    # FOV center in comet frame at start
    fov_start_dra = (pnt_ra - comet_start_ra) * cos_dec * 3600
    fov_start_ddec = (pnt_dec - comet_start_dec) * 3600
    
    # FOV center in comet frame at end
    fov_end_dra = (pnt_ra - comet_end_ra) * cos_dec * 3600
    fov_end_ddec = (pnt_dec - comet_end_dec) * 3600
    
    # FOV center in comet frame at midpoint
    comet_mid_ra, comet_mid_dec = float(trk_ra[mid]), float(trk_dec[mid])
    fov_mid_dra = (pnt_ra - comet_mid_ra) * cos_dec * 3600
    fov_mid_ddec = (pnt_dec - comet_mid_dec) * 3600
    
    print(f"FOV center (comet frame offsets from mosaic center, arcsec):")
    print(f"  At obs START: East={fov_start_dra:+.1f}, North={fov_start_ddec:+.1f}")
    print(f"  At obs MID:   East={fov_mid_dra:+.1f}, North={fov_mid_ddec:+.1f}")
    print(f"  At obs END:   East={fov_end_dra:+.1f}, North={fov_end_ddec:+.1f}")
    
    # XMM FOV radius ~15 arcmin = 900 arcsec
    fov_r = 900
    print(f"\n  XMM FOV radius: ~{fov_r}\"")
    print(f"  FOV right edge at START: East={fov_start_dra + fov_r:.1f}\"")
    print(f"  FOV left edge at END:    East={fov_end_dra - fov_r:.1f}\"")
    print(f"  Track START position:    East={dx_track[0]:+.1f}\"")
    print(f"  Track END position:      East={dx_track[-1]:+.1f}\"")
    
    print(f"\n  Gap: Track start ({dx_track[0]:+.1f}\") vs FOV right edge at start ({fov_start_dra + fov_r:.1f}\")")
    print(f"  Gap: Track end ({dx_track[-1]:+.1f}\") vs FOV left edge at end ({fov_end_dra - fov_r:.1f}\")")
    
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
