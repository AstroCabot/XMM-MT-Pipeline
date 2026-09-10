import sys

import numpy as np

from _common import DATA, REPRO, read_json, stage_dir, write_json

sys.path.insert(0, str(REPRO))

from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

from analysis_02_pps_luminosity import run as pps

PPS = DATA / "pps"


def main():
    output = stage_dir(2)
    gtis = pps.read_gtis(PPS / pps.GTI_IMAGE)
    sources, ephemeris, event_header = pps.motion_inputs(
        event=PPS / "events-header-only.fits",
        sources_path=PPS / "sources.tsv",
        ephemeris_path=PPS / "ephemeris.tsv",
    )
    times, weights = pps.time_samples(gtis)
    epoch = Time(
        np.average(pps.event_mjd(times, event_header), weights=weights),
        format="mjd",
        scale="utc",
    )
    if abs((epoch - pps.EXPECTED_EPOCH).sec) > 0.01:
        raise ValueError("unexpected GTI-weighted epoch")

    maps = [np.asarray(fits.getdata(PPS / name), float) for name in pps.EXPMAPS]
    coverage = pps.exposure_coverage(
        gtis,
        event_header,
        sources,
        ephemeris,
        maps[0].shape,
        pps.PARAMS["source_mask_radius_arcsec"],
    )
    header = fits.getheader(PPS / pps.EXPMAPS[0])
    header["RADESYS"] = header.pop("RADECSYS")
    center = WCS(header, fix=False).all_world2pix([pps.REFERENCE], 0)[0]
    if np.max(np.abs(center - pps.MAP_CENTER)) > 1e-4:
        raise ValueError("unexpected PPS map center")
    yy, xx = np.indices(maps[0].shape)
    radius = np.hypot(xx - pps.MAP_CENTER[0], yy - pps.MAP_CENTER[1]) * pps.MAP_PIXEL
    source_region = radius <= pps.APERTURE
    background_region = (radius >= pps.BACKGROUND[0]) & (radius < pps.BACKGROUND[1])
    pixel_area = (pps.MAP_PIXEL / 60) ** 2

    archived = read_json(PPS / "archived_result.json")
    rows, checks = [], []
    for (lo, hi), image, stored in zip(pps.BANDS, maps, archived["bands"], strict=True):
        if stored["energy_eV"] != [lo, hi]:
            raise ValueError("archived band order changed")
        exposed = image * coverage
        source_area = float(exposed[source_region].sum() * pixel_area)
        background_area = float(exposed[background_region].sum() * pixel_area)
        self_subtraction = pps.correction(
            exposed, radius, source_region, background_region
        )
        mask_fill = float(
            np.sum(image[source_region] / radius[source_region])
            / np.sum(exposed[source_region] / radius[source_region])
        )
        checks.append(
            {
                "energy_eV": [lo, hi],
                "recomputed_vs_archived": {
                    "source_exposure_s_arcmin2": [
                        source_area,
                        stored["source_exposure_s_arcmin2"],
                    ],
                    "annulus_self_subtraction_factor_1_over_r": [
                        self_subtraction,
                        stored["annulus_self_subtraction_factor_1_over_r"],
                    ],
                    "source_mask_fill_factor_1_over_r": [
                        mask_fill,
                        stored["source_mask_fill_factor_1_over_r"],
                    ],
                },
            }
        )
        rows.append(
            {
                **stored,
                "source_exposure_s_arcmin2": source_area,
                "background_exposure_s_arcmin2": background_area,
                "annulus_self_subtraction_factor_1_over_r": self_subtraction,
                "source_mask_fill_factor_1_over_r": mask_fill,
            }
        )

    worst = max(
        abs(new / old - 1)
        for check in checks
        for new, old in check["recomputed_vs_archived"].values()
    )
    if worst > 1e-6:
        raise RuntimeError(f"recomputed PPS geometry deviates from archive: {worst}")

    result = pps.make_result(epoch, rows)
    write_json(output / "result.json", result)
    luminosity = result["flux_luminosity"]["full_aperture"]["luminosity_erg_s"]
    print(f"stage 02: epoch {epoch.isot}, L_X(full aperture) = {luminosity:.3e} erg/s")
    print(f"stage 02: recomputed geometry matches archive to {worst:.2e}")


if __name__ == "__main__":
    main()
