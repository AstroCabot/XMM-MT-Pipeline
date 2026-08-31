import contextlib
import math
import os

import numpy as np

import spectrum
from model import (
    DETECTORS,
    EXPRESSION,
    components,
    configure,
    free_parameters,
    global_index,
    link,
)

AU_CM = 1.495978707e13
KEV_ERG = 1.602176634e-9


def minimise(xspec):
    previous = math.inf
    for _ in range(6):
        xspec.Fit.renorm()
        xspec.Fit.perform()
        value = float(xspec.Fit.statistic)
        if previous - value < 0.01:
            break
        previous = value
    return value


def profile_errors(xspec, targets):
    models = [
        xspec.AllModels(index, "joint") for index in range(1, xspec.AllData.nGroups + 1)
    ]
    free_pars = free_parameters(models)
    saved = [float(parameter.values[0]) for parameter in free_pars]
    centers = {
        name: float(parameter.values[0]) for name, (_, parameter) in targets.items()
    }
    try:
        command = "maximum 3.0 1.0 " + " ".join(
            global_index(model, parameter) for model, parameter in targets.values()
        )
        for attempt in range(2):
            try:
                xspec.Fit.error(command)
                break
            except Exception as error:
                for parameter, value in zip(free_pars, saved):
                    parameter.values = value
                minimise(xspec)
                if attempt:
                    raise RuntimeError("XSPEC profile errors failed") from error
        intervals = {
            name: [float(parameter.error[0]), float(parameter.error[1])]
            for name, (_, parameter) in targets.items()
        }
        for name, (low, high) in intervals.items():
            if not np.isfinite([low, high]).all() or not low <= centers[name] <= high:
                raise RuntimeError(f"invalid XSPEC profile interval: {name}")
    finally:
        for parameter, value in zip(free_pars, saved):
            parameter.values = value
        minimise(xspec)
    return intervals


@contextlib.contextmanager
def isolated(models, region_models, keep):
    first = components(models[0])
    norms = {"lhb": first[2].norm, "halo": first[4].norm, "cxb": first[5].norm}
    norms.update(
        {
            f"cx_{name}": components(model)[6].norm
            for name, model in region_models.items()
        }
    )
    saved = {name: list(parameter.values[:2]) for name, parameter in norms.items()}
    for name, parameter in norms.items():
        on = name == keep or (keep == "cx" and name.startswith("cx_"))
        parameter.values = [saved[name][0] if on else 0, -1]
    try:
        yield
    finally:
        for name, parameter in norms.items():
            parameter.values = saved[name]


def band_text(band):
    low, high = map(float, band)
    if not 0 < low < high:
        raise ValueError("invalid fit band")
    return f"{low:g} {high:g}"


def fluxes(xspec, count, band):
    xspec.AllModels.calcFlux(band_text(band))
    return [
        (float(xspec.AllData(index).flux[0]), float(xspec.AllData(index).flux[3]))
        for index in range(1, count + 1)
    ]


def mean_energy_eV(energy_flux, photon_flux):
    if energy_flux <= 0 or photon_flux <= 0:
        raise ValueError("non-positive intrinsic VACX2 flux")
    return energy_flux / photon_flux / KEV_ERG * 1000


def profiled_fluxes(xspec, models, region_models, groups, intervals, band):
    states = [
        (parameter, list(parameter.values), bool(parameter.frozen))
        for parameter in free_parameters(models)
    ]
    output = {}
    try:
        for name, model in region_models.items():
            parameter = components(model)[6].norm
            values = []
            index = next(
                i
                for i, group in enumerate(groups)
                if group["region"] == name and group["det"] == "pn"
            )
            for endpoint in intervals[f"cx_norm_{name}"]:
                for item, saved, frozen in states:
                    item.values, item.frozen = saved, frozen
                parameter.values, parameter.frozen = max(float(endpoint), 0), True
                minimise(xspec)
                with isolated(models, region_models, "cx"):
                    intrinsic = fluxes(xspec, len(groups), band)
                surface = intrinsic[index][0] / groups[index]["omega_arcmin2"]
                values.append(surface * groups[index]["full_area_arcmin2"])
            output[name] = values
    finally:
        for parameter, values, frozen in states:
            parameter.values, parameter.frozen = values, frozen
        minimise(xspec)
    return output


def plot_arrays(xspec, count):
    xspec.Plot.device, xspec.Plot.xAxis = "/null", "keV"
    xspec.Plot.background, xspec.Plot.add = False, False
    xspec.Plot("data")
    arrays = [
        {
            "energy_keV": np.asarray(xspec.Plot.x(i), float),
            "half_width_keV": np.asarray(xspec.Plot.xErr(i), float),
            "rate_ct_s_keV": np.asarray(xspec.Plot.y(i), float),
            "rate_error_ct_s_keV": np.asarray(xspec.Plot.yErr(i), float),
            "model_ct_s_keV": np.asarray(xspec.Plot.model(i), float),
        }
        for i in range(1, count + 1)
    ]
    xspec.Plot("delchi")
    for i, entry in enumerate(arrays, 1):
        entry["residual_sigma"] = np.asarray(xspec.Plot.y(i), float)
    return arrays


def component_arrays(xspec, models, region_models, count, arrays):
    for name in ("lhb", "halo", "cxb", "cx"):
        with isolated(models, region_models, name):
            xspec.Plot("data")
            for index, entry in enumerate(arrays, 1):
                entry[f"{name}_ct_s_keV"] = np.asarray(xspec.Plot.model(index), float)


def conversion(xspec, models, region_models, groups, band):
    with isolated(models, region_models, "cx"):
        intrinsic = fluxes(xspec, len(groups), band)
        selected = []
        try:
            for index, group in enumerate(groups, 1):
                low, high = spectrum.grouped_range(
                    group["source_pha"], band[0] * 1000, band[1] * 1000
                )
                data = xspec.AllData(index)
                data.notice("all")
                data.ignore(f"**-{low - 1} {high + 1}-**")
                selected.append(high - low + 1)
            folded = plot_arrays(xspec, len(groups))
            if any(
                len(entry["energy_keV"]) != count
                for entry, count in zip(folded, selected)
            ):
                raise RuntimeError("XSPEC image-band selection differs from the PHA")
        finally:
            xspec.AllData.notice("all")
            xspec.AllData.ignore("bad")
    factors = {}
    for index, group in enumerate(groups):
        rate = float(
            np.sum(
                folded[index]["model_ct_s_keV"] * 2 * folded[index]["half_width_keV"]
            )
        )
        camera = float(components(models[index])[1].factor.values[0])
        factors.setdefault(group["region"], {})[group["det"]] = (
            rate * camera / intrinsic[index][0]
        )
    energy_flux, photon_flux = intrinsic[
        next(i for i, group in enumerate(groups) if group["region"] == "inner")
    ]
    return factors, mean_energy_eV(energy_flux, photon_flux), intrinsic


def profile_targets(models, detector_models, region_models):
    base, cx = components(models[0]), components(models[0])[6]
    targets = {
        "lhb_norm": (models[0], base[2].norm),
        "halo_norm": (models[0], base[4].norm),
        "cxb_norm": (models[0], base[5].norm),
        "cx_temperature_keV": (models[0], cx.temperature),
        "cx_collision_keV_u": (models[0], cx.collnpar),
    }
    for detector in ("mos1", "mos2"):
        if detector in detector_models:
            fitted = detector_models[detector]
            targets[f"crosscal_{detector}"] = (fitted, components(fitted)[1].factor)
    targets.update(
        {
            f"cx_norm_{name}": (fitted, components(fitted)[6].norm)
            for name, fitted in region_models.items()
        }
    )
    return targets


def fit(groups, settings, distance_au, atomdb, band, profile=False):
    os.environ["ATOMDB"] = atomdb
    import xspec

    __import__("acx2_xspec")
    xspec.AllData.clear()
    xspec.AllModels.clear()
    xspec.Xset.chatter = xspec.Xset.logChatter = 0
    xspec.Xset.abund, xspec.Xset.xsect = "wilm", "vern"
    xspec.Fit.statMethod, xspec.Fit.query, xspec.Fit.nIterations = "chi", "yes", 300
    groups = sorted(
        groups,
        key=lambda group: (
            ("inner", "outer").index(group["region"]),
            DETECTORS.index(group["det"]),
        ),
    )
    cwd = os.getcwd()
    os.chdir(groups[0]["source_pha"].parent)
    try:
        xspec.AllData(
            " ".join(
                f"{i}:{i} {group['source_pha'].name}"
                for i, group in enumerate(groups, 1)
            )
        )
        for index, group in enumerate(groups, 1):
            data = xspec.AllData(index)
            data.background = group["background_pha"].name
            data.response = group["rmf"].name
            data.response.arf = group["arf"].name
        xspec.AllData.ignore("bad")
        xspec.Model(EXPRESSION, "joint", 1)
        models = [
            xspec.AllModels(index, "joint") for index in range(1, len(groups) + 1)
        ]
        for group, model in zip(groups, models):
            configure(model, group["det"], group["omega_arcmin2"], settings)
        detector_models, region_models = link(groups, models)
        first = components(models[0])[6]
        best = minimise(xspec)
        saved = [float(parameter.values[0]) for parameter in free_parameters(models)]
        first.collnpar.values = settings["cx_collision_starts_keV_u"][1]
        alternative = minimise(xspec)
        if alternative >= best:
            for parameter, value in zip(free_parameters(models), saved):
                parameter.values = value
            best = minimise(xspec)
        targets = profile_targets(models, detector_models, region_models)
        intervals = profile_errors(xspec, targets) if profile else {}
        crosscal = {
            detector: float(components(model)[1].factor.values[0])
            for detector, model in detector_models.items()
        }
        region_factors, mean_energy, intrinsic = conversion(
            xspec, models, region_models, groups, band
        )
        flux_intervals = (
            profiled_fluxes(xspec, models, region_models, groups, intervals, band)
            if profile
            else {}
        )
        factors = region_factors["inner"]
        regions = {}
        for name, region_model in region_models.items():
            members = [i for i, group in enumerate(groups) if group["region"] == name]
            i = next(i for i in members if groups[i]["det"] == "pn")
            surface = intrinsic[i][0] / groups[i]["omega_arcmin2"]
            full_area = groups[members[0]]["full_area_arcmin2"]
            live_area = sum(
                groups[i]["exposure_s"] * groups[i]["omega_arcmin2"] for i in members
            ) / sum(groups[i]["exposure_s"] for i in members)
            flux = surface * full_area
            norm = float(components(region_model)[6].norm.values[0])
            regions[name] = {
                "surface_brightness_erg_cm2_s_arcmin2": surface,
                "flux_erg_cm2_s": flux,
                "luminosity_erg_s": 4 * np.pi * (distance_au * AU_CM) ** 2 * flux,
                "full_solid_angle_arcmin2": full_area,
                "live_solid_angle_arcmin2": live_area,
                "cx_norm_per_arcmin2": norm,
                "flux_1sigma_erg_cm2_s": flux_intervals.get(name),
                "luminosity_1sigma_erg_s": [
                    value * 4 * np.pi * (distance_au * AU_CM) ** 2
                    for value in flux_intervals.get(name, [])
                ],
            }
        result = {
            "statistic": float(xspec.Fit.statistic),
            "dof": int(xspec.Fit.dof),
            "crosscal": crosscal,
            "n_free_parameters": len(free_parameters(models)),
            "parameters": {
                "lhb_norm": float(components(models[0])[2].norm.values[0]),
                "halo_norm": float(components(models[0])[4].norm.values[0]),
                "cxb_norm": float(components(models[0])[5].norm.values[0]),
                "cx_temperature_keV": float(first.temperature.values[0]),
                "cx_collision_keV_u": float(first.collnpar.values[0]),
            },
            "profile_intervals_1sigma": intervals,
            "regions": regions,
            "detector_K": factors,
            "region_detector_K": region_factors,
            "mean_photon_energy_eV": mean_energy,
        }
        if profile:
            arrays = plot_arrays(xspec, len(groups))
            component_arrays(xspec, models, region_models, len(groups), arrays)
            result["arrays"] = arrays
            result["groups"] = groups
        return result
    finally:
        xspec.AllData.clear()
        xspec.AllModels.clear()
        os.chdir(cwd)
