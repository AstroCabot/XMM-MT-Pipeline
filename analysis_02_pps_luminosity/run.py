import json
import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

from .selection import (
    EVENT_CENTER,
    MAP_CENTER,
    EVENT_PIXEL,
    MAP_PIXEL,
    correction,
    event_mjd,
    exposure_coverage,
    gti_expression,
    read_gtis,
    read_table,
    source_keep,
    time_samples,
)

HERE = Path(__file__).resolve().parent
REPRO = HERE.parent
REDUCTION = REPRO / "analysis_01_reduction"
PARAMS = json.loads((REPRO / "parameters.json").read_text())
SETTINGS = json.loads((REDUCTION / "settings.json").read_text())
MT_PPS = Path(SETTINGS["mt_pps"])
OUTPUT = HERE / "outputs"
APERTURE = float(PARAMS["science_aperture_radius_arcsec"])
BACKGROUND = tuple(map(float, PARAMS["pps_background_annulus_arcsec"]))
SOFT_BAND = tuple(map(int, PARAMS["soft_band_ev"]))

EVENT = "P0963720201PNS003PIEVLI0000.FTZ"
GTI_IMAGE = "P0963720201PNS003IMAGE_2000.FTZ"
EXPMAPS = (
    "P0963720201PNS003EXPMAP1000.FTZ",
    "P0963720201PNS003EXPMAP2000.FTZ",
    "P0963720201PNS003EXPMAP3000.FTZ",
)
RMF = HERE / "calibration/epn_e4_ff20_sdY9_v22.0.rmf"
SOURCES = REDUCTION / "outputs/04_sources/sources.tsv"
EPHEMERIS = REDUCTION / "outputs/07_comove/ephemeris.tsv"
CCF = REDUCTION / "outputs/01_init/ccf.cif"

BANDS = ((SOFT_BAND[0], 500), (500, 1000), (1000, SOFT_BAND[1]))
FINAL_PRODUCTS = (
    "result.json",
    "response.arf",
    *(f"{stem}-{lo}-{hi}.pi" for lo, hi in BANDS for stem in ("source", "background")),
)
REFERENCE = (177.705, 1.43536)
AU_CM = 1.495978707e13
DISTANCE_AU = 1.8815323353
EXPECTED_EPOCH = Time("2025-12-03T21:25:46.398", scale="utc")


def _inputs():
    paths = [MT_PPS / EVENT, MT_PPS / GTI_IMAGE]
    paths += [MT_PPS / name for name in EXPMAPS]
    summaries = list((REDUCTION / "outputs/01_init/odf").glob("*SUM.SAS"))
    paths += [RMF, SOURCES, EPHEMERIS, CCF, *summaries]
    missing = [str(path) for path in paths if not path.is_file()]
    if len(summaries) != 1:
        raise FileNotFoundError("expected one private ODF summary")
    if missing:
        raise FileNotFoundError("missing inputs:\n" + "\n".join(missing))
    return paths


def motion_inputs(event=MT_PPS / EVENT, sources_path=SOURCES, ephemeris_path=EPHEMERIS):
    source_rows = read_table(sources_path)
    sources = [(float(row["ra"]), float(row["dec"])) for row in source_rows]
    eph_rows = read_table(ephemeris_path)
    ephemeris = np.array(
        [[float(row[key]) for key in ("mjd", "ra_deg", "dec_deg")] for row in eph_rows]
    )
    if (
        len(sources) != 15
        or len(set(sources)) != 15
        or ephemeris.shape[0] < 2
        or np.any(~np.isfinite(ephemeris))
        or np.any(np.diff(ephemeris[:, 0]) <= 0)
    ):
        raise ValueError("unexpected source catalog or ephemeris")
    with fits.open(event, memmap=True) as hdul:
        if not hdul[0].header.get("MTFLAG"):
            raise ValueError("PPS event list is not in the moving-target frame")
        header = hdul["EVENTS"].header.copy()
    return sources, ephemeris, header


def integrate_flux(path, lo, hi):
    with fits.open(path, memmap=True) as hdul:
        hdu = next(
            h
            for h in hdul
            if getattr(h.data, "dtype", None) is not None
            and h.data.dtype.names
            and {"ENERGY", "ENERGY_BIN", "FLUX"} <= set(h.data.dtype.names)
        )
        energy = np.asarray(hdu.data["ENERGY"], float)
        width = np.asarray(hdu.data["ENERGY_BIN"], float)
        flux = np.asarray(hdu.data["FLUX"], float)
    step = np.nanmedian(np.diff(energy))
    half = (
        width / 2
        if abs(np.nanmedian(width) - step) <= abs(2 * np.nanmedian(width) - step)
        else width
    )
    overlap = np.maximum(
        0, np.minimum(energy + half, hi) - np.maximum(energy - half, lo)
    )
    return float(np.sum(flux * overlap))


def summarize_flux(rows):
    raw = sum(row["raw_masked_flux_erg_cm-2_s-1"] for row in rows)
    corrected = sum(
        row["raw_masked_flux_erg_cm-2_s-1"]
        * row["annulus_self_subtraction_factor_1_over_r"]
        for row in rows
    )
    full = sum(
        row["raw_masked_flux_erg_cm-2_s-1"]
        * row["annulus_self_subtraction_factor_1_over_r"]
        * row["source_mask_fill_factor_1_over_r"]
        for row in rows
    )
    if raw <= 0 or corrected <= 0:
        raise ValueError("nonpositive integrated flux")
    luminosity_factor = 4 * math.pi * (DISTANCE_AU * AU_CM) ** 2
    fluxes = {
        "raw_masked": raw,
        "annulus_self_subtraction_corrected_masked": corrected,
        "full_aperture": full,
    }
    return {
        "flux_luminosity": {
            name: {
                "flux_erg_cm-2_s-1": flux,
                "luminosity_erg_s": flux * luminosity_factor,
            }
            for name, flux in fluxes.items()
        },
        "correction_factors": {
            "annulus_self_subtraction_exposure_weighted_1_over_r": corrected / raw,
            "source_mask_fill_1_over_r": full / corrected,
            "raw_masked_to_full_aperture_1_over_r": full / raw,
        },
    }


def make_result(epoch, rows):
    return {
        "epoch_utc": epoch.isot,
        "distance_au": DISTANCE_AU,
        "selection": {
            "energy_eV": list(SOFT_BAND),
            "aperture_arcsec": APERTURE,
            "background_arcsec": BACKGROUND,
            "event_filter": "exact CCD GTIs; RAWY>12; PATTERN<=4; FLAG==0",
            "source_mask_radius_arcsec": PARAMS["source_mask_radius_arcsec"],
        },
        "response": {
            "arf": "response.arf",
            "arf_model": "flat extended source",
            "rmf": f"calibration/{RMF.name}",
        },
        "bands": rows,
        **summarize_flux(rows),
        "units": {
            "counts": "count",
            "background_scale": "dimensionless",
            "flux": "erg cm-2 s-1",
            "luminosity": "erg s-1",
        },
    }


def _environment(work):
    command = (
        f'export HEADAS="{SETTINGS["headas"]}" && '
        'source "$HEADAS/headas-init.sh" && '
        f'source "{SETTINGS["sas"]}" && '
        f'export LD_LIBRARY_PATH="{SETTINGS["python_lib"]}:${{LD_LIBRARY_PATH:-}}" && '
        "env -0"
    )
    result = subprocess.run(["bash", "-lc", command], capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stderr.decode(errors="replace"))
    env = dict(
        item.split("=", 1) for item in result.stdout.decode().split("\0") if "=" in item
    )
    summary = next((REDUCTION / "outputs/01_init/odf").glob("*SUM.SAS"))
    env.update(
        SAS_CCFPATH=SETTINGS["ccf"],
        SAS_CCF=str(CCF),
        SAS_ODF=str(summary),
        TMPDIR=str(work),
        XDG_CACHE_HOME=str(work),
        XDG_CONFIG_HOME=str(work),
    )
    pfiles = work / "pfiles"
    pfiles.mkdir()
    env["PFILES"] = f"{pfiles};{env.get('PFILES', '').split(';')[-1]}"
    return env


def _sas(tool, work, env, *parameters):
    result = subprocess.run(
        [tool, *parameters], cwd=work, env=env, text=True, capture_output=True
    )
    text = result.stdout + result.stderr
    errors = [
        line
        for line in text.splitlines()
        if line.startswith("** ") and ": error" in line.lower()
    ]
    if result.returncode or errors:
        raise RuntimeError(f"{tool} failed:\n" + "\n".join(text.splitlines()[-20:]))


def _set_pha(path, counts, backscal):
    with fits.open(path, mode="update", memmap=False) as hdul:
        hdu = hdul["SPECTRUM"]
        if len(counts) != len(hdu.data):
            raise ValueError("count array does not match the PHA channel grid")
        hdu.data["COUNTS"] = np.asarray(counts, hdu.data["COUNTS"].dtype)
        hdu.header.update(
            POISSERR=True,
            BACKSCAL=float(backscal),
            AREASCAL=1.0,
            RESPFILE=f"../calibration/{RMF.name}",
            ANCRFILE="response.arf",
            BACKFILE="none",
        )


def run():
    OUTPUT.mkdir(exist_ok=True)
    for name in FINAL_PRODUCTS:
        (OUTPUT / name).unlink(missing_ok=True)
    inputs = _inputs()
    witness = [(path.stat().st_size, path.stat().st_mtime_ns) for path in inputs]
    scratch = HERE / "work"
    scratch.mkdir(exist_ok=True)
    gtis = read_gtis(MT_PPS / GTI_IMAGE)
    sources, ephemeris, event_header = motion_inputs()
    times, weights = time_samples(gtis)
    epoch = Time(
        np.average(event_mjd(times, event_header), weights=weights),
        format="mjd",
        scale="utc",
    )
    if abs((epoch - EXPECTED_EPOCH).sec) > 0.01:
        raise ValueError("unexpected GTI-weighted epoch")

    with tempfile.TemporaryDirectory(dir=scratch) as name:
        work = Path(name)
        for short, path in (
            ("events.FTZ", MT_PPS / EVENT),
            ("gti.FTZ", MT_PPS / GTI_IMAGE),
            (RMF.name, RMF),
        ):
            (work / short).symlink_to(path)
        env = _environment(work)
        quality = (
            f"{gti_expression()}&&(RAWY>12)&&(PATTERN<=4)&&(FLAG==0)"
            "&&(PI>=0)&&(PI<=20479)"
        )
        _sas(
            "evselect",
            work,
            env,
            "table=events.FTZ:EVENTS",
            "withfilteredset=yes",
            "filteredset=selected.fits",
            "keepfilteroutput=yes",
            "destruct=yes",
            "filterexposure=yes",
            "updateexposure=yes",
            "writedss=yes",
            f"expression={quality}",
        )
        source_region = (
            f"circle({EVENT_CENTER[0]},{EVENT_CENTER[1]},{APERTURE / EVENT_PIXEL})"
        )
        background_region = (
            f"annulus({EVENT_CENTER[0]},{EVENT_CENTER[1]},"
            f"{BACKGROUND[0] / EVENT_PIXEL},{BACKGROUND[1] / EVENT_PIXEL})"
        )
        for stem, region in (
            ("source", source_region),
            ("background", background_region),
        ):
            _sas(
                "evselect",
                work,
                env,
                "table=selected.fits:EVENTS",
                "withspectrumset=yes",
                f"spectrumset={stem}.pi",
                "energycolumn=PI",
                "spectralbinsize=5",
                "specchannelmin=0",
                "specchannelmax=20479",
                "withspecranges=yes",
                "destruct=yes",
                "filterexposure=yes",
                "updateexposure=yes",
                "writedss=yes",
                f"expression=((X,Y) IN {region})",
            )
        _sas(
            "arfgen",
            work,
            env,
            "spectrumset=source.pi",
            "arfset=response.arf",
            "withrmfset=no",
            "detmaptype=flat",
            "extendedsource=yes",
            "withbadpixcorr=yes",
            "badpixlocation=events.FTZ",
            "setbackscale=no",
            "filterdss=yes",
            "modelootcorr=no",
        )

        channels = len(fits.getdata(work / "source.pi", "SPECTRUM"))
        if (
            len(fits.getdata(work / "background.pi", "SPECTRUM")) != channels
            or len(fits.getdata(RMF, "EBOUNDS")) != channels
        ):
            raise ValueError("PHA and RMF channel grids differ")

        with fits.open(work / "selected.fits", memmap=True) as hdul:
            events, header = hdul["EVENTS"].data, hdul["EVENTS"].header.copy()
            keep = source_keep(
                events["X"],
                events["Y"],
                event_mjd(events["TIME"], header),
                sources,
                ephemeris,
                PARAMS["source_mask_radius_arcsec"],
            )
            radius_event = (
                np.hypot(events["X"] - EVENT_CENTER[0], events["Y"] - EVENT_CENTER[1])
                * EVENT_PIXEL
            )
            event_pi = np.asarray(events["PI"], int)
            event_channel = event_pi // 5
            source_events = keep & (radius_event <= APERTURE)
            background_events = (
                keep & (radius_event >= BACKGROUND[0]) & (radius_event < BACKGROUND[1])
            )
            source_counts = np.bincount(
                event_channel[source_events], minlength=channels
            )
            background_counts = np.bincount(
                event_channel[background_events], minlength=channels
            )
            if len(source_counts) != channels or len(background_counts) != channels:
                raise ValueError("event channel exceeds the PHA grid")

        maps = [np.asarray(fits.getdata(MT_PPS / name), float) for name in EXPMAPS]
        if len({image.shape for image in maps}) != 1:
            raise ValueError("PPS exposure-map shapes differ")
        coverage = exposure_coverage(
            gtis,
            event_header,
            sources,
            ephemeris,
            maps[0].shape,
            PARAMS["source_mask_radius_arcsec"],
        )
        header = fits.getheader(MT_PPS / EXPMAPS[0])
        header["RADESYS"] = header.pop("RADECSYS")
        center = WCS(header, fix=False).all_world2pix([REFERENCE], 0)[0]
        if np.max(np.abs(center - MAP_CENTER)) > 1e-4:
            raise ValueError("unexpected PPS map center")
        yy, xx = np.indices(maps[0].shape)
        radius = np.hypot(xx - MAP_CENTER[0], yy - MAP_CENTER[1]) * MAP_PIXEL
        source_region = radius <= APERTURE
        background_region = (radius >= BACKGROUND[0]) & (radius < BACKGROUND[1])
        pixel_area = (MAP_PIXEL / 60) ** 2
        source_exposure = float(
            fits.getheader(work / "source.pi", "SPECTRUM")["EXPOSURE"]
        )
        background_exposure = float(
            fits.getheader(work / "background.pi", "SPECTRUM")["EXPOSURE"]
        )

        rows = []
        for index, ((lo, hi), image) in enumerate(zip(BANDS, maps, strict=True)):
            exposed = image * coverage
            source_area = float(exposed[source_region].sum() * pixel_area)
            background_area = float(exposed[background_region].sum() * pixel_area)
            source_backscal = source_area / source_exposure
            background_backscal = background_area / background_exposure
            alpha = source_backscal / background_backscal
            source_name, background_name = (
                f"source-{lo}-{hi}.pi",
                f"background-{lo}-{hi}.pi",
            )
            shutil.copyfile(work / "source.pi", work / source_name)
            shutil.copyfile(work / "background.pi", work / background_name)
            _set_pha(work / source_name, source_counts, source_backscal)
            _set_pha(work / background_name, background_counts, background_backscal)
            fluxed = work / f"flux-{lo}-{hi}.fits"
            _sas(
                "efluxer",
                work,
                env,
                f"spectrumset={source_name}",
                f"backgndset={background_name}",
                "withbgdset=yes",
                "arfset=response.arf",
                f"rmfset={RMF.name}",
                f"fluxedset={fluxed.name}",
            )
            net_flux = integrate_flux(fluxed, lo / 1000, hi / 1000)
            background_correction = correction(
                exposed, radius, source_region, background_region
            )
            mask_fill = np.sum(image[source_region] / radius[source_region]) / np.sum(
                exposed[source_region] / radius[source_region]
            )
            band = (event_pi >= lo) & (
                (event_pi <= hi) if index == len(BANDS) - 1 else (event_pi < hi)
            )
            rows.append(
                {
                    "energy_eV": [lo, hi],
                    "source_counts": int(np.count_nonzero(source_events & band)),
                    "background_counts": int(
                        np.count_nonzero(background_events & band)
                    ),
                    "background_scale": alpha,
                    "source_exposure_s_arcmin2": source_area,
                    "background_exposure_s_arcmin2": background_area,
                    "raw_masked_flux_erg_cm-2_s-1": net_flux,
                    "annulus_self_subtraction_factor_1_over_r": background_correction,
                    "source_mask_fill_factor_1_over_r": mask_fill,
                }
            )
        if witness != [
            (path.stat().st_size, path.stat().st_mtime_ns) for path in inputs
        ]:
            raise RuntimeError("an input changed during the run")
        for lo, hi in BANDS:
            for stem in ("source", "background"):
                name = f"{stem}-{lo}-{hi}.pi"
                shutil.copyfile(work / name, OUTPUT / name)
        shutil.copyfile(work / "response.arf", OUTPUT / "response.arf")

    result = make_result(epoch, rows)
    (OUTPUT / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    run()
