import hashlib
import json
import os
from urllib.request import urlretrieve

import numpy as np
from cdflib import CDF
from scipy.optimize import brentq
from scipy.special import ndtr

from lightcurve import PARAMS, REPRO, SETTINGS
from solar_wind import CACHE

CONTEXT = PARAMS["ion_fraction_context"]
SPECTRUM = REPRO / "analysis_03_spectrum"
ACE_URL = (
    "https://spdf.gsfc.nasa.gov/pub/data/ace/swics/level_2_cdaweb/sw2_h3/2025/"
    "ac_h3_sw2_{}_v04.cdf"
)
HIS_URL = (
    "https://soar.esac.esa.int/soar-sl-tap/data?retrieval_type=LAST_PRODUCT"
    "&product_type=SCIENCE&data_item_id=solo_L3_swa-his-comp-10min_{}"
)
ACE_FIELDS = (
    "vHe2",
    "C6to5",
    "C6to5_qual",
    "O7to6",
    "O7to6_qual",
    "O8to6",
    "O8to6_qual",
)
HIS_FIELDS = (
    "QUALITY_FLAG",
    "HE2_SPEED",
    "C_O_ABUN",
    "MG_O_ABUN",
    "SI_O_ABUN",
    "FE_O_ABUN",
    "O_CHARGE_DIST",
    "C_CHARGE_DIST",
    "MG_CHARGE_DIST",
    "O7_O6_RATIO",
)
FILL_LIMIT = -1e30
C65_LIMITS = (0.1104, 1.425)
O76_FLOOR = 0.052272
O86_FLOOR = 0.00080294
ACX_ELEMENTS = [1, 2, 6, 7, 8, 10, 12, 13, 14, 16, 18, 20, 26, 28]
ACX_STEP_EV = 0.5


def days(first, stop):
    dates = np.arange(first, stop, dtype="datetime64[D]")
    return [str(date).replace("-", "") for date in dates]


ACE_DAYS = days("2025-01-01", "2025-03-17")
HIS_DAYS = days("2023-04-01", "2023-04-08") + days("2023-04-12", "2023-04-28")


def retrieve(url, label, dates):
    paths = []
    for date in dates:
        path = CACHE / "composition" / f"{label}_{date}.cdf"
        if not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_suffix(".tmp")
            urlretrieve(url.format(date), temporary)
            temporary.replace(path)
        paths.append(path)
    return paths


def read(paths, fields):
    rows = {field: [] for field in fields}
    for path in paths:
        with CDF(path) as cdf:
            for field in fields:
                rows[field].append(cdf.varget(field))
    return {field: np.concatenate(values) for field, values in rows.items()}


def c6_o(c65):
    return np.where(
        c65 >= 0.31818,
        0.90788 - 0.45079 * np.exp(-0.081265 * c65) - 0.43521 * np.exp(-0.64476 * c65),
        10.0 ** (-0.42166 + 1.0475 * np.log10(np.maximum(c65, 1e-12))),
    )


def c5_o(c65):
    return 0.026529 + 0.22179 * np.exp(-0.16456 * c65) + 0.16355 * np.exp(-1.2468 * c65)


def o8_o(o86):
    return np.where(
        o86 >= 0.00284,
        0.62398 - 0.56883 * np.exp(-0.60366 * o86) - 0.054758 * np.exp(-7.82875 * o86),
        10.0 ** (-0.05933 + 0.99345 * np.log10(np.maximum(o86, 1e-12))),
    )


def o7_o(o76):
    return np.where(
        o76 >= 0.0475,
        0.85269 - 0.43445 * np.exp(-0.22965 * o76) - 0.41785 * np.exp(-2.0811 * o76),
        10.0 ** (-0.063384 + 0.97403 * np.log10(np.maximum(o76, 1e-12))),
    )


def ne9_o(o76):
    return 0.69869 - 0.69827 * np.exp(-0.033402 * o76)


def mg11_o(o76):
    return 0.00037163 + 0.0074849 * o76**1.9113


def nitrogen_o(speed):
    return np.where(speed < 500, 0.064, 0.011)


def o_h(speed):
    return -0.0028735 + 0.0012157 * np.log10(speed)


def finite(data, fields):
    return np.logical_and.reduce(
        [np.isfinite(data[field]) & (data[field] > FILL_LIMIT) for field in fields]
    )


def ace_fractions(data):
    keep = finite(data, ("vHe2", "C6to5", "O7to6", "O8to6"))
    for field in ("C6to5_qual", "O7to6_qual", "O8to6_qual"):
        keep &= data[field] == 0
    keep &= (
        (data["C6to5"] >= C65_LIMITS[0])
        & (data["C6to5"] <= C65_LIMITS[1])
        & (data["O7to6"] > O76_FLOOR)
        & (data["O8to6"] > O86_FLOOR)
    )
    speed, c65, o76, o86 = (
        data[field][keep] for field in ("vHe2", "C6to5", "O7to6", "O8to6")
    )
    return (
        c5_o(c65)
        + c6_o(c65)
        + o7_o(o76)
        + o8_o(o86)
        + ne9_o(o76)
        + mg11_o(o76)
        + nitrogen_o(speed)
    ) * o_h(speed)


def his_fractions(data):
    keep = finite(
        data, ("HE2_SPEED", "C_O_ABUN", "MG_O_ABUN", "SI_O_ABUN", "FE_O_ABUN")
    )
    keep &= data["QUALITY_FLAG"] == 0
    for field in ("O_CHARGE_DIST", "C_CHARGE_DIST", "MG_CHARGE_DIST"):
        keep &= np.all(np.isfinite(data[field]) & (data[field] >= 0), axis=1)
    keep &= np.isfinite(data["O7_O6_RATIO"]) & (data["O7_O6_RATIO"] > 0)
    speed = data["HE2_SPEED"][keep]
    oxygen, carbon, magnesium = (
        data[field][keep]
        for field in ("O_CHARGE_DIST", "C_CHARGE_DIST", "MG_CHARGE_DIST")
    )
    return (
        data["C_O_ABUN"][keep] * (carbon[:, 1] + carbon[:, 2])
        + oxygen[:, 2]
        + oxygen[:, 3]
        + data["MG_O_ABUN"][keep] * magnesium[:, 5]
        + ne9_o(data["O7_O6_RATIO"][keep])
        + nitrogen_o(speed)
    ) * o_h(speed)


def percentiles(base):
    sigma = CONTEXT["o_h_scatter_fwhm_dex"] / (2 * np.sqrt(2 * np.log(2))) * np.log(10)
    locations = np.log(base)
    limits = locations.min() - 8 * sigma, locations.max() + 8 * sigma

    def quantile(level):
        cdf = lambda value: np.mean(ndtr((value - locations) / sigma)) - level
        return float(np.exp(brentq(cdf, *limits)))

    return [quantile(level) for level in (0.16, 0.50, 0.84)]


def tabulate(paths, fields, fractions):
    data = read(paths, fields)
    values = fractions(data)
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.read_bytes())
    return {
        "records": {"all": len(data[fields[0]]), "retained": int(values.size)},
        "input_sha256": digest.hexdigest(),
        "f_16_50_84": percentiles(values),
    }


def band_energy(model, collision, centers):
    spectrum = model.calc_spectrum(collision, Tbroaden=0, vbroaden=0)
    return float(np.sum(spectrum * centers))


def selected_emission():
    os.environ["ATOMDB"] = SETTINGS["atomdb"]
    import acx2

    result = json.loads((SPECTRUM / "outputs/result.json").read_text())
    settings = json.loads((SPECTRUM / "parameters.json").read_text())["model"]
    temperature = result["fit"]["parameters"]["cx_temperature_keV"]
    collision = result["fit"]["parameters"]["cx_collision_keV_u"]
    helium = settings["cx_donor_helium_fraction"]
    low, high = PARAMS["soft_band_ev"]
    edges = np.arange(low, high + ACX_STEP_EV / 2, ACX_STEP_EV) / 1000
    centers = (edges[1:] + edges[:-1]) / 2
    model = acx2.ACXModel()
    for donor in ("H", "He"):
        model.add_donor(donor, elements=ACX_ELEMENTS)
    model.set_abund(np.ones(len(ACX_ELEMENTS)), elements=ACX_ELEMENTS)
    model.set_ebins(edges)
    model.set_temperature(temperature)
    model.set_acxmodel(result["model"]["acxmodel"])
    model.set_recombtype(result["model"]["recombtype"])
    model.set_donorabund(["H", "He"], [1 - helium, helium])
    model.set_collisiontype(result["model"]["collntype"], "kev/amu")
    total = band_energy(model, collision, centers)
    kept = {element: np.zeros_like(value) for element, value in model.ionfrac.items()}
    for element, charges in CONTEXT["selected_incident_charges_by_Z"].items():
        kept[int(element)][charges] = model.ionfrac[int(element)][charges]
    model.set_ionfrac(kept)
    return {
        "band_eV": [low, high],
        "energy_step_eV": ACX_STEP_EV,
        "cx_temperature_keV": temperature,
        "cx_collision_keV_u": collision,
        "cx_donor_helium_fraction": helium,
        "acxmodel": result["model"]["acxmodel"],
        "recombtype": result["model"]["recombtype"],
        "collntype": result["model"]["collntype"],
        "energy_flux_fraction": band_energy(model, collision, centers) / total,
    }


def calculate():
    return {
        "ace_swics": tabulate(
            retrieve(ACE_URL, "ace", ACE_DAYS), ACE_FIELDS, ace_fractions
        ),
        "solar_orbiter_his": tabulate(
            retrieve(HIS_URL, "his", HIS_DAYS), HIS_FIELDS, his_fractions
        ),
        "selected_emission": selected_emission(),
    }
