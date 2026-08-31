import contextlib
import os
from pathlib import Path

import numpy as np
from astropy.io import fits

import pha
from background import mask_sky_map
from common import CFG, directory, finish, read_tsv, sas, write_tsv
from flare import replace_hard_map

STAGE = "06_soft_protons"
BREAK_KEV = 3.0
NH_1E22 = 0.020
CXB_INDEX = 1.46
CXB_NORM_ARCMIN2 = 8.88e-7
SP_NORM_FLOOR_ARCMIN2 = 3e-6
HARD_EV = tuple(CFG["bands"]["hard"]["pi"])
LINES = {
    "mos1": (1.486, 1.740),
    "mos2": (1.486, 1.740),
    "pn": (1.486, 8.041),
}


def xspec_module():
    try:
        import xspec
    except ImportError as error:
        raise RuntimeError("run after sourcing the configured HEASoft setup") from error
    xspec.Xset.chatter = 0
    xspec.Xset.logChatter = 0
    xspec.Xset.abund = "wilm"
    xspec.Xset.xsect = "vern"
    return xspec


@contextlib.contextmanager
def xspec_session(path):
    xspec, cwd = xspec_module(), Path.cwd()
    xspec.AllData.clear()
    xspec.AllModels.clear()
    try:
        os.chdir(Path(path).resolve().parent)
        yield xspec
    finally:
        xspec.AllData.clear()
        xspec.AllModels.clear()
        os.chdir(cwd)


def source_path(row):
    if row["det"] != "pn":
        return Path(row["source_pha"])
    prefix = Path(row["source_pha"]).name.removesuffix("-fovt.pi")
    path = Path(row["dir"]) / f"{prefix}-fovtootsub.pi"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def frame_spectrum(row):
    return pha.net(source_path(row), row["qpb_pha"])


def sky_model(xspec, detector, area, low, high):
    energies = [energy for energy in LINES[detector] if low <= energy <= high]
    if low < 1:
        energies += [0.561, 0.654]
    expression = "apec+tbabs*(apec+powerlaw)" + "+gaussian" * len(energies)
    model = xspec.Model(expression, "sky", 1)
    model.apec.kT.values = [0.10, -1]
    model.apec.Abundanc.values = [1, -1]
    model.apec.Redshift.values = [0, -1]
    model.TBabs.nH.values = [NH_1E22, -1]
    model.apec_3.kT.values = [0.25, -1]
    model.apec_3.Abundanc.values = [1, -1]
    model.apec_3.Redshift.values = [0, -1]
    for component in (model.apec, model.apec_3):
        component.norm.values = [1e-4, 0.01, 0, 0, 1, 1] if low < 1 else [0, -1]
    model.powerlaw.PhoIndex.values = [CXB_INDEX, -1]
    model.powerlaw.norm.values = [CXB_NORM_ARCMIN2 * area, -1]
    gaussians = [
        getattr(model, name)
        for name in model.componentNames
        if name.startswith("gaussian")
    ]
    for component, energy in zip(gaussians, energies):
        component.LineE.values = [energy, -1]
        component.Sigma.values = [0, -1]
        component.norm.values = [1e-5, 0.01, 0, 0, 10, 10]
    return model


def fit_master(path, rmf, arf, diagonal, detector, band, area):
    path, rmf, arf, diagonal = map(
        lambda item: Path(item).resolve(), (path, rmf, arf, diagonal)
    )
    low, high = (
        (2.5, 8.5) if band == "hard" else (0.5, 5.5) if detector == "pn" else (0.5, 9.0)
    )
    with xspec_session(path) as xspec:
        spectrum = xspec.Spectrum(path.name)
        spectrum.multiresponse[0] = str(rmf)
        spectrum.multiresponse[0].arf = str(arf)
        spectrum.multiresponse[1] = str(diagonal)
        xspec.AllData.ignore("bad")
        spectrum.ignore(f"**-{low} {high}-**")
        sky_model(xspec, detector, area, low, high)
        proton = xspec.Model("bknpow", "proton", 2)
        proton.bknpower.PhoIndx1.values = [0.8, 0.01, 0.1, 0.1, 1.4, 1.4]
        proton.bknpower.BreakE.values = [BREAK_KEV, -1]
        proton.bknpower.PhoIndx2.values = [0.8, 0.01, 0.1, 0.1, 5.0, 5.0]
        proton.bknpower.norm.values = [1e-3, 0.01, 0, 0, 1e4, 1e4]
        xspec.Fit.statMethod = "chi"
        xspec.Fit.query = "no"
        xspec.Fit.nIterations = 200
        xspec.Fit.renorm()
        xspec.Fit.perform()
        values = {
            "fit_low_keV": low,
            "fit_high_keV": high,
            "gamma_low": float(proton.bknpower.PhoIndx1.values[0]),
            "gamma_high": float(proton.bknpower.PhoIndx2.values[0]),
            "break_keV": BREAK_KEV,
            "norm": float(proton.bknpower.norm.values[0]),
            "chi2": float(xspec.Fit.statistic),
            "dof": int(xspec.Fit.dof),
        }
        if (
            not all(np.isfinite(value) for value in values.values())
            or values["dof"] <= 0
            or values["chi2"] / values["dof"] > 2
        ):
            raise RuntimeError(f"invalid soft-proton fit for {detector}/{band}")
        if values["norm"] < SP_NORM_FLOOR_ARCMIN2 * area:
            values["norm"] = 0.0
        return values


def make_master(detector, band, rows):
    work = directory(STAGE) / "master" / detector / band
    work.mkdir(parents=True, exist_ok=True)
    spectra = [frame_spectrum(row) for row in rows]
    master = pha.stack(spectra)
    source_counts = pha.stack(pha.read(source_path(row)) for row in rows).counts
    net_path = work / "net.pi"
    pha.write(
        net_path, source_path(rows[0]), master, rate=True, group_counts=source_counts
    )
    weights = [item.exposure * item.backscal for item in spectra]  # Diffuse counts.
    rmf, arf = work / "response.rmf", work / "response.arf"
    pha.add_response([row["rmf"] for row in rows], weights, rmf, "addrmf")
    pha.add_response([row["arf"] for row in rows], weights, arf, "addarf")
    diagonal = Path(CFG["esas"]) / f"{detector}-diag.rsp"
    if not diagonal.exists():
        raise FileNotFoundError(diagonal)
    result = fit_master(
        net_path, rmf, arf, diagonal, detector, band, pha.area_arcmin2(master)
    )
    hard_rate = pha.band_counts(master, *HARD_EV) / master.exposure
    if hard_rate <= 0:
        raise RuntimeError(f"non-positive hard residual for {detector}/{band}")
    pattern = CFG["bands"][band]["pn" if detector == "pn" else "mos"]["pattern"]
    result.update(
        det=detector,
        band=band,
        pattern=pattern,
        n_frames=len(rows),
        hard_net_rate=hard_rate,
        diagonal=diagonal,
    )
    return result


def safe_index(value):
    # proton is singular at index 1.
    return 1.001 if abs(value - 1) < 1e-3 else value


def fake_spectrum(path, diagonal, fit, norm, template):
    path, diagonal = Path(path).resolve(), Path(diagonal).resolve()
    with xspec_session(path) as xspec:
        model = xspec.Model("bknpow")
        model.bknpower.PhoIndx1.values = [fit["gamma_low"], -1]
        model.bknpower.BreakE.values = [BREAK_KEV, -1]
        model.bknpower.PhoIndx2.values = [fit["gamma_high"], -1]
        model.bknpower.norm.values = [max(norm, 1e-30), -1]
        settings = xspec.FakeitSettings(
            response=str(diagonal),
            arf="",
            exposure=template.exposure,
            fileName=path.name,
        )
        xspec.AllData.fakeit(1, [settings], applyStats=False, noWrite=False)
        generated = pha.read(path)
        if generated.channel.size != template.channel.size:
            raise ValueError("soft-proton and source channel counts differ")
        counts = pha.response_counts(generated.counts, diagonal, template)
        generated = generated._replace(
            channel=template.channel,
            counts=counts,
            channel_ev=template.channel_ev,
            variance=np.zeros_like(counts),
            backscal=template.backscal,
        )
        pha.write(path, path, generated)


def frame_norm(spectrum, fit):
    # FOV spectra retain hard PI channels in soft-band jobs.
    rate = pha.band_counts(spectrum, *HARD_EV) / spectrum.exposure
    scale = max(rate, 0) / fit["hard_net_rate"] if fit["hard_net_rate"] > 0 else 0
    return scale, fit["norm"] * scale


def make_frame(row, fit, sources):
    work = directory(STAGE) / row["frame"] / row["band"]
    work.mkdir(parents=True, exist_ok=True)
    spectrum = frame_spectrum(row)
    scale, norm = frame_norm(spectrum, fit)
    low, high = CFG["bands"][row["band"]]["pi"]
    prefix = Path(row["source_pha"]).name.removesuffix("-fovt.pi")
    detector_image = Path(row["dir"]) / f"{prefix}-fovimdet-{low}-{high}.fits"
    input_spectrum = source_path(row)
    for source in (detector_image, input_spectrum):
        link = work / source.name
        if not link.exists():
            link.symlink_to(source)
    detector_map = work / f"{row['frame']}-{row['band']}-sp-det.fits"
    sky_map = work / f"{row['frame']}-{row['band']}-sp-sky.fits"
    detector_map.unlink(missing_ok=True)
    sky_map.unlink(missing_ok=True)
    proton_parameters = {
        "imagefile": detector_image.name,
        "specfile": input_spectrum.name,
        "ccds": row["chips"],
        "speccontrol": 2,
        "bindl": f"{safe_index(fit['gamma_low']):.6f}",
        "bindh": f"{safe_index(fit['gamma_high']):.6f}",
        "bbreak": BREAK_KEV,
        "bnorm": f"{max(norm, 1e-30):.8e}",
    }
    sas(
        "proton", work, spmapdet=detector_map, elow=low, ehigh=high, **proton_parameters
    )
    if row["band"] == "hard":
        if (low, high) != (2500, 8500):
            raise RuntimeError("hard FLARE correction requires the 2500-8500 eV band")
        support = None
        if "F" in row["chips"].split():
            support = work / f".{row['frame']}-hard-sp-support.fits"
            support.unlink(missing_ok=True)
            try:
                sas(
                    "proton",
                    work,
                    spmapdet=support,
                    elow=300,
                    ehigh=12000,
                    **proton_parameters,
                )
                replace_hard_map(detector_map, row["det"], row["chips"], support)
            finally:
                support.unlink(missing_ok=True)
        else:
            replace_hard_map(detector_map, row["det"], row["chips"])
    sas(
        "rotdet2sky",
        work,
        intemplate=row["events_image"],
        inimage=detector_map,
        outimage=sky_map,
    )
    detector_counts = float(np.nansum(fits.getdata(detector_map)))
    sky_counts = float(np.nansum(fits.getdata(sky_map)))
    if not np.isclose(sky_counts, detector_counts, rtol=0.01, atol=0.05):
        raise RuntimeError(
            f"soft-proton reprojection mismatch for {row['frame']}/{row['band']}"
        )
    sp_pha = work / f"{row['frame']}-{row['band']}-sp.pi"
    sp_pha.unlink(missing_ok=True)
    fake_spectrum(sp_pha, fit["diagonal"], fit, norm, spectrum)
    mask_sky_map(sky_map, sources)
    return {
        "frame": row["frame"],
        "det": row["det"],
        "band": row["band"],
        "sp_pha": sp_pha,
        "sp_map": sky_map,
        "gamma_low": f"{fit['gamma_low']:.5f}",
        "gamma_high": f"{fit['gamma_high']:.5f}",
        "break_keV": f"{BREAK_KEV:.1f}",
        "norm": f"{norm:.8e}",
        "hard_rate_scale": f"{scale:.6f}",
    }


def run():
    os.environ.setdefault("ATOMDB", CFG["atomdb"])
    directory(STAGE)
    xspec_module()
    rows = read_tsv(directory("05_background") / "background.tsv")
    grouped = {}
    for row in rows:
        grouped.setdefault((row["det"], row["band"]), []).append(row)
    models = {key: make_master(*key, value) for key, value in sorted(grouped.items())}
    fit_columns = (
        "det",
        "band",
        "pattern",
        "n_frames",
        "fit_low_keV",
        "fit_high_keV",
        "gamma_low",
        "gamma_high",
        "break_keV",
        "norm",
        "chi2",
        "dof",
        "hard_net_rate",
    )
    write_tsv(
        directory(STAGE) / "fits.tsv",
        fit_columns,
        [{key: fit[key] for key in fit_columns} for fit in models.values()],
    )
    sources = read_tsv(directory("04_sources") / "sources.tsv")
    output = [
        make_frame(row, models[(row["det"], row["band"])], sources) for row in rows
    ]
    write_tsv(directory(STAGE) / "soft_protons.tsv", tuple(output[0]), output)
    finish(STAGE)
