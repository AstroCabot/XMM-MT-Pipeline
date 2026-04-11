#!/usr/bin/env python3
"""Predict the offset between comet emission and track overlay in the QC plot."""
import math

# The mosaic center (from FITS header) = apparent coords from track midpoint
ctr_ra = 177.94642639160156
ctr_dec = 1.323562741279602

# Astrometric position of comet at midpoint (from diagnostic)
astro_ra = 177.614551
astro_dec = 1.467528

# How events appear: events are in ICRS. A photon FROM the comet at mid-obs
# has ICRS RA/DEC ~ astrometric position. The mosaic offsets it from the
# apparent-coordinate center:
dx_event = (astro_ra - ctr_ra) * math.cos(math.radians(ctr_dec)) * 3600
dy_event = (astro_dec - ctr_dec) * 3600

# Where does the track midpoint plot? It's (apparent - apparent) = 0
dx_track = 0.0
dy_track = 0.0

print(f"Mosaic center (apparent):       RA={ctr_ra:.6f}, DEC={ctr_dec:.6f}")
print(f"Comet mid-obs (astrometric):    RA={astro_ra:.6f}, DEC={astro_dec:.6f}")
print()
print(f"Event (actual comet) offset in mosaic: dx={dx_event:+.1f}, dy={dy_event:+.1f} arcsec")
print(f"Track midpoint offset in mosaic:       dx={dx_track:.1f}, dy={dy_track:.1f} arcsec")
print(f"Track appears {-dx_event:.0f} arcsec EAST and {-dy_event:.0f} arcsec SOUTH of actual emission")
print(f"Total separation: {math.sqrt(dx_event**2 + dy_event**2):.0f} arcsec = "
      f"{math.sqrt(dx_event**2 + dy_event**2)/60:.1f} arcmin")
