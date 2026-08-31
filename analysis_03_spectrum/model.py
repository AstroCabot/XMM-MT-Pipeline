EXPRESSION = "constant*constant*(apec+tbabs*(apec+powerlaw)+vacx2)"
DETECTORS = ("pn", "mos1", "mos2")


def components(model):
    return [getattr(model, name) for name in model.componentNames]


def free(parameter, start, bounds):
    low, high = bounds
    parameter.values = [start, 0.01 * max(abs(start), high - low), low, low, high, high]


def fixed(parameter, value):
    parameter.values = [value, -1]


def choice(parameter, value):
    parameter.values = value
    parameter.frozen = True


def configure(model, detector, area, settings):
    scale, camera, lhb, absorption, halo, cxb, cx = components(model)
    fixed(scale.factor, area)
    if detector == "pn":
        fixed(camera.factor, 1)
    else:
        free(camera.factor, 1, settings["crosscal_bounds"])
    for plasma, temperature, norm in (
        (lhb, settings["lhb_kT_keV"], settings["lhb_norm_start"]),
        (halo, settings["halo_kT_keV"], settings["halo_norm_start"]),
    ):
        fixed(plasma.kT, temperature)
        fixed(plasma.Abundanc, 1)
        fixed(plasma.Redshift, 0)
        free(plasma.norm, norm, (0, 0.01))
    fixed(absorption.nH, settings["nh_1e22_cm2"])
    fixed(cxb.PhoIndex, settings["cxb_index"])
    free(cxb.norm, settings["cxb_norm_start"], (0, 0.01))
    free(
        cx.temperature,
        settings["cx_temperature_start_keV"],
        settings["cx_temperature_bounds_keV"],
    )
    free(
        cx.collnpar,
        settings["cx_collision_starts_keV_u"][0],
        settings["cx_collision_bounds_keV_u"],
    )
    choice(cx.collntype, 1)
    choice(cx.acxmodel, 8)
    choice(cx.recombtype, 1)
    fixed(cx.Hefrac, settings["cx_donor_helium_fraction"])
    for name in (
        "H",
        "He",
        "C",
        "N",
        "O",
        "Ne",
        "Mg",
        "Al",
        "Si",
        "S",
        "Ar",
        "Ca",
        "Fe",
        "Ni",
    ):
        fixed(getattr(cx, name), 1)
    for name in ("vbroad", "tbroad", "Redshift"):
        fixed(getattr(cx, name), 0)
    free(cx.norm, 0.01, (0, 100))


def link(groups, models):
    reference = components(models[0])
    detector_reference, region_reference = {}, {}
    for group, fitted in zip(groups, models):
        detector_reference.setdefault(group["det"], fitted)
        region_reference.setdefault(group["region"], fitted)
    for group, fitted in zip(groups, models):
        _, camera, lhb, _, halo, cxb, cx = components(fitted)
        if fitted is not models[0]:
            for item, source in ((lhb, reference[2]), (halo, reference[4])):
                item.norm.link = source.norm
            cxb.norm.link = reference[5].norm
            cx.temperature.link, cx.collnpar.link = (
                reference[6].temperature,
                reference[6].collnpar,
            )
        region_model = region_reference[group["region"]]
        if fitted is not region_model:
            cx.norm.link = components(region_model)[6].norm
        detector_model = detector_reference[group["det"]]
        if group["det"] == "pn":
            fixed(camera.factor, 1)
        elif fitted is not detector_model:
            camera.factor.link = components(detector_model)[1].factor
    return detector_reference, region_reference


def free_parameters(models):
    return [
        fitted(index)
        for fitted in models
        for index in range(1, fitted.nParameters + 1)
        if not fitted(index).frozen and not fitted(index).link
    ]


def global_index(fitted, parameter):
    return f"joint:{fitted.startParIndex + parameter.index - 1}"
