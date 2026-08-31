import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits

HERE = Path(__file__).resolve().parent
REPRO = HERE.parent
SETTINGS = json.loads((REPRO / "analysis_01_reduction/settings.json").read_text())
BACKGROUND = REPRO / "analysis_01_reduction/outputs/05_background"
PARAMS = json.loads((REPRO / "parameters.json").read_text())
SIZE = 900  # binadapt reads a fixed 900x900 buffer
BINFACTOR = 2
ARCMIN2_PER_DEG2 = 3600.0
SOFT_BAND = tuple(map(int, PARAMS["soft_band_ev"]))
SOFT_TAG = f"{SOFT_BAND[0]}-{SOFT_BAND[1]}"


def environment(work):
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
    pfiles = work / "pfiles"
    pfiles.mkdir()
    env.update(
        PFILES=f"{pfiles};{env.get('PFILES', '').split(';')[-1]}",
        TMPDIR=str(work),
        XDG_CACHE_HOME=str(work),
        XDG_CONFIG_HOME=str(work),
        SAS_CCFPATH=SETTINGS["ccf"],
        SAS_CCF=str(REPRO / "analysis_01_reduction/outputs/01_init/ccf.cif"),
    )
    return env


def run(command, work):
    result = subprocess.run(
        command, cwd=work, env=environment(work), text=True, capture_output=True
    )
    if result.returncode:
        signal = (
            f" (killed by signal {-result.returncode})" if result.returncode < 0 else ""
        )
        raise RuntimeError(
            f"{command[0]} exited {result.returncode}{signal}\n"
            + "\n".join((result.stdout + result.stderr).splitlines()[-40:])
        )


def write_plane(path, data, template, header, exposure=None):
    """Write a padded binadapt input with comet-frame WCS."""
    rows, columns = data.shape
    if rows > SIZE or columns > SIZE:
        raise ValueError(f"binadapt input {rows}x{columns} exceeds {SIZE}")
    x0, y0 = (SIZE - columns) // 2, (SIZE - rows) // 2
    padded = np.zeros((SIZE, SIZE), data.dtype)
    padded[y0 : y0 + rows, x0 : x0 + columns] = data
    with fits.open(template, memmap=False) as source:
        hdul = fits.HDUList([hdu.copy() for hdu in source])
    hdul[0].data = padded
    hdul[0].header["CRPIX1"] = float(header["CRPIX1"]) + x0
    hdul[0].header["CRPIX2"] = float(header["CRPIX2"]) + y0
    for key in ("CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "CTYPE1", "CTYPE2"):
        hdul[0].header[key] = header[key]
    if exposure is not None:
        hdul[0].header["EXPOSURE"] = exposure
        hdul[0].header["ONTIME"] = exposure
    hdul.writeto(path)
    return x0, y0


def adaptive(events, qpb, sp, oot, exposure, header, mask=None):
    """Run binadapt with 50-count kernels on a padded 900-pixel input."""
    work_root = REPRO / "work"
    work_root.mkdir(exist_ok=True)
    template = next(
        iter(sorted(BACKGROUND.glob(f"mos1-*/soft/mos1S001-fovimsky-{SOFT_TAG}.fits"))),
        None,
    )
    if template is None:
        raise FileNotFoundError("no ESAS image template for binadapt")
    planes = {
        "fovimsky": np.asarray(events, np.float32),
        "bkgimsky": np.asarray(qpb + PARAMS["pn_oot_fraction"] * oot, np.float32),
        "protimsky": np.asarray(sp, np.float32),
        "expimsky": np.asarray(exposure, np.float32),
    }
    if mask is not None:
        planes["maskimsky"] = np.asarray(mask, np.int32)
    with tempfile.TemporaryDirectory(prefix="binadapt-", dir=work_root) as name:
        work = Path(name)
        mask_options = (
            ["withmask=no"]
            if mask is None
            else [
                "withmask=yes",
                f"maskfile={work / f'comb-maskimsky-{SOFT_TAG}.fits'}",
                "maskthresh=0.5",
            ]
        )
        offsets = {
            write_plane(
                work / f"comb-{stem}-{SOFT_TAG}.fits",
                data,
                template,
                header,
                float(exposure.max()) if stem == "expimsky" else None,
            )
            for stem, data in planes.items()
        }
        ((x0, y0),) = offsets
        run(
            [
                "binadapt",
                "prefix=comb",
                f"elow={SOFT_BAND[0]}",
                f"ehigh={SOFT_BAND[1]}",
                "withpartbkg=yes",
                "withspbkg=yes",
                "withswcxbkg=no",
                *mask_options,
                "withbinning=yes",
                f"binfactor={BINFACTOR}",
                "withsmoothing=yes",
                "smoothcounts=50",
            ],
            work,
        )
        blocked = np.asarray(
            fits.getdata(work / f"comb-adaptimsky-{SOFT_TAG}.fits"), float
        )
    output = surface_brightness(blocked, events.shape, x0, y0)
    if mask is not None:
        output[~mask] = np.nan
    return output


def surface_brightness(blocked, shape, x0, y0):
    """Return binadapt output on the input grid in count/s/arcmin2."""
    padded = np.repeat(np.repeat(blocked, BINFACTOR, 0), BINFACTOR, 1)
    return padded[y0 : y0 + shape[0], x0 : x0 + shape[1]] / ARCMIN2_PER_DEG2
