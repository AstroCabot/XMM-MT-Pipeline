from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.signal import fftconvolve

ANNULUS_CLOSURE_RTOL = 2e-5


def clean_wcs_header(header):
    header = header.copy()
    if "RADECSYS" in header:
        header["RADESYS"] = header.pop("RADECSYS")
    return header


def image_radii(shape, header, center):
    matrix = np.asarray(WCS(header, fix=False).pixel_scale_matrix, float)
    pixel_area_arcmin2 = abs(float(np.linalg.det(matrix))) * 3600
    rows, columns = np.indices(shape)
    offset = matrix @ np.vstack(
        ((columns - center[0]).ravel(), (rows - center[1]).ravel())
    )
    radius = np.hypot(offset[0], offset[1]).reshape(shape) * 60
    return radius, pixel_area_arcmin2


def annulus_exposure(image, header, annulus_edges_arcsec):
    header = clean_wcs_header(header)
    center = WCS(header, fix=False).celestial.wcs.crpix - 1
    radius, pixel_area = image_radii(image.shape, header, center)
    edges = np.asarray(annulus_edges_arcsec, float) / 60
    ring = np.searchsorted(edges, radius, side="right") - 1
    values = np.asarray(image, float)
    use = (ring >= 0) & (ring < len(edges) - 1) & (values > 0)
    return np.bincount(ring[use], values[use] * pixel_area, minlength=len(edges) - 1)


def deposit(shape, x, y, values):
    """Bilinearly deposit point values onto an image."""
    image = np.zeros(shape, float)
    ix, iy = np.floor(x).astype(int), np.floor(y).astype(int)
    fx, fy = x - ix, y - iy
    for dx, dy, weight in (
        (0, 0, (1 - fx) * (1 - fy)),
        (1, 0, fx * (1 - fy)),
        (0, 1, (1 - fx) * fy),
        (1, 1, fx * fy),
    ):
        xx, yy = ix + dx, iy + dy
        use = (xx >= 0) & (xx < shape[1]) & (yy >= 0) & (yy < shape[0])
        np.add.at(image, (yy[use], xx[use]), values[use] * weight[use])
    return image


def comoving_exposure(
    image,
    source_header,
    track_ra_deg,
    track_dec_deg,
    sample_weights,
    output_header,
    output_shape,
):
    """Reproject one exposure through the comet-motion kernel."""
    source_wcs = WCS(clean_wcs_header(source_header), fix=False).celestial
    output_wcs = WCS(clean_wcs_header(output_header), fix=False).celestial
    image = np.asarray(image, float)
    keep = np.isfinite(image) & (image != 0)
    y, x = np.nonzero(keep)
    ra, dec = source_wcs.pixel_to_world_values(x, y)
    sx, sy = output_wcs.world_to_pixel_values(ra, dec)
    cx, cy = output_wcs.world_to_pixel_values(track_ra_deg, track_dec_deg)
    weights = np.asarray(sample_weights, float)
    base_x, base_y = np.average(cx, weights=weights), np.average(cy, weights=weights)
    center = (np.asarray(output_shape[::-1]) - 1) / 2
    factor = abs(np.linalg.det(source_wcs.pixel_scale_matrix)) / abs(
        np.linalg.det(output_wcs.pixel_scale_matrix)
    )
    base = deposit(
        output_shape,
        sx - base_x + center[0],
        sy - base_y + center[1],
        image[keep] * factor,
    )
    dx, dy = -(cx - base_x), -(cy - base_y)
    radius = int(np.ceil(max(np.max(abs(dx)), np.max(abs(dy))))) + 2
    kernel = deposit(
        (2 * radius + 1,) * 2, radius + dx, radius + dy, weights / weights.sum()
    )
    output = fftconvolve(base, kernel, mode="same")
    expected = float(np.nansum(image) * factor)
    if abs(output.sum() - expected) > 1e-5 * max(abs(expected), 1):
        raise RuntimeError("exposure integral changed in comet reprojection")
    return output


def radial_exposure(image, header, maximum_radius_arcmin):
    header = clean_wcs_header(header)
    wcs = WCS(header, fix=False)
    radius, pixel_area = image_radii(image.shape, header, wcs.wcs.crpix - 1)
    step = np.sqrt(pixel_area)
    count = int(np.ceil(maximum_radius_arcmin / step)) + 1
    index = np.floor(radius / step).astype(int)
    use = index < count
    pixels = np.bincount(index[use], minlength=count)
    radius_sum = np.bincount(index[use], radius[use], minlength=count)
    exposure_sum = np.bincount(
        index[use], np.asarray(image, float)[use], minlength=count
    )
    valid = pixels > 0
    return radius_sum[valid] / pixels[valid], exposure_sum[valid] / pixels[valid]


def task04_plane(path, soft_band_ev):
    image, header = fits.getdata(Path(path), "EXPOSURE", header=True)
    image = np.asarray(image, float)
    numerical_floor = -1e-10 * max(float(np.max(image)), 1)
    if (
        header.get("DETECTOR") != "EPIC"
        or header.get("BAND") != "soft"
        or [header.get("E_MIN"), header.get("E_MAX")] != list(soft_band_ev)
        or not header.get("COMOVING", False)
        or image.ndim != 2
        or np.any(~np.isfinite(image))
        or np.any(image < numerical_floor)
        or not np.any(image > 0)
    ):
        raise RuntimeError("invalid Task 04 MOS2-equivalent exposure plane")
    image = image.copy()
    image[image < 0] = 0
    return image, header


def fractional_closure(actual, expected, label):
    actual, expected = np.asarray(actual, float), np.asarray(expected, float)
    if (
        actual.shape != expected.shape
        or np.any(~np.isfinite(actual))
        or np.any(~np.isfinite(expected))
        or np.any(expected <= 0)
    ):
        raise RuntimeError(f"invalid {label} annular exposure")
    fraction = actual / expected - 1
    if np.max(np.abs(fraction)) > ANNULUS_CLOSURE_RTOL:
        raise RuntimeError(f"{label} annular exposure did not close")
    return fraction


def composite_density(frame_density, frame_weights):
    density, weights = np.asarray(frame_density), np.asarray(frame_weights)
    if (
        density.ndim != 2
        or weights.ndim != 2
        or density.shape[0] != weights.shape[0]
        or np.any(~np.isfinite(density))
        or np.any(density < 0)
        or np.any(~np.isfinite(weights))
        or np.any(weights < 0)
        or np.any(weights.sum(axis=0) <= 0)
    ):
        raise ValueError("invalid frame PSF weights")
    return np.asarray(
        [
            np.average(density, axis=0, weights=weights[:, index])
            for index in range(weights.shape[1])
        ]
    )


def annulus_quadrature(
    annulus_edges_arcsec, radial_order, exposure_radius_arcmin=None, exposure_mean=None
):
    if (exposure_radius_arcmin is None) != (exposure_mean is None):
        raise ValueError("exposure radius and mean must be supplied together")
    edges = np.asarray(annulus_edges_arcsec, float) / 60
    nodes, legendre = np.polynomial.legendre.leggauss(radial_order)
    output = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        radius = 0.5 * (hi + lo) + 0.5 * (hi - lo) * nodes
        weight = legendre * radius
        if exposure_radius_arcmin is not None:
            weight *= np.interp(radius, exposure_radius_arcmin, exposure_mean)
        if np.any(~np.isfinite(weight)) or np.any(weight < 0) or weight.sum() <= 0:
            raise ValueError("invalid annular exposure quadrature")
        output.append((radius, weight / weight.sum()))
    return output


def radial_operator(
    psf_radius_arcsec,
    density_arcsec2,
    source_radius_arcmin,
    annulus_edges_arcsec,
    radial_order,
    azimuth_order,
    exposure_radius_arcmin=None,
    exposure_mean=None,
):
    """Map intrinsic radial brightness to exposure-weighted annulus means."""
    source = np.asarray(source_radius_arcmin, float)
    step = float(np.median(np.diff(source)))
    psf_r = np.asarray(psf_radius_arcsec, float) / 60
    density = np.asarray(density_arcsec2, float).copy()
    if density.ndim == 1:
        density = np.tile(density, (len(annulus_edges_arcsec) - 1, 1))
    if density.shape != (len(annulus_edges_arcsec) - 1, len(psf_r)):
        raise ValueError("PSF density does not match the annuli and radial grid")
    density *= 3600
    angle = 2 * np.pi * (np.arange(azimuth_order) + 0.5) / azimuth_order
    rows = []
    quadrature = annulus_quadrature(
        annulus_edges_arcsec, radial_order, exposure_radius_arcmin, exposure_mean
    )
    for row_index, (rho, radial_weight) in enumerate(quadrature):
        row = np.zeros_like(source)
        for value, weight in zip(rho, radial_weight):
            distance = np.sqrt(
                value**2
                + source[:, None] ** 2
                - 2 * value * source[:, None] * np.cos(angle)
            )
            kernel = np.interp(
                distance, psf_r, density[row_index], left=density[row_index, 0], right=0
            )
            row += weight * source * step * kernel.mean(axis=1) * 2 * np.pi
        rows.append(row)
    matrix = np.asarray(rows)
    if (
        np.any(~np.isfinite(matrix))
        or np.any(matrix < 0)
        or np.any(matrix.sum(axis=1) <= 0)
    ):
        raise ValueError("invalid radial PSF operator")
    return matrix / matrix.sum(axis=1)[:, None]
