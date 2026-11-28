import hashlib
import json
import tarfile
from urllib.request import urlretrieve

import numpy as np
from astropy.io import fits

from .data import HERE, SHARED, pattern_max

URL = "https://heasarc.gsfc.nasa.gov/FTP/xmm/data/responses/xmm_ao26_calfiles.tar.gz"
CCF_URL = "https://heasarc.gsfc.nasa.gov/FTP/xmm/data/CCF/"
SINGLE_RMF = "epn_e4_ff20_sY9_v22.0.rmf"
QE_FILE = "EPN_QUANTUMEF_0018.CCF"


def channel_weights(channels, low, high):
    # PPS PI is integer eV; RMF channel n contains PI=5n,...,5n+4.
    start = 5 * np.asarray(channels, float)
    return np.clip(
        (np.minimum(start + 5, high + 1) - np.maximum(start, low + 1)) / 5, 0, 1
    )


def fold_singles(photons, cache):
    arf = fits.getdata(cache / "pn-thin-5.arf", "SPECRESP")
    with fits.open(cache / QE_FILE) as ccf:
        table = ccf["FRACTION_ENERGY"].data
        row = table[
            (table["MODE_ID"] == 0)
            & (table["THRESH"] == 20)
            & (table["REGION"] == "((RAWX,RAWY) IN RECTANGLE(1,181,64,200,0))")
        ][0]
        energy = ccf["EBINS_FRACTION"].data[0]["ENERGY"]
        middle = 500 * (arf["ENERG_LO"] + arf["ENERG_HI"])
        area = arf["SPECRESP"] * np.interp(
            middle, energy, row["FRAC_S"] / row["FRAC_SD"]
        )
    with fits.open(cache / SINGLE_RMF) as hdul:
        matrix = hdul["MATRIX"].data
        assert np.allclose(matrix["ENERG_LO"], arf["ENERG_LO"])
        assert np.allclose(matrix["ENERG_HI"], arf["ENERG_HI"])
        counts = np.zeros(len(hdul["EBOUNDS"].data))
        for scale, row in zip(photons * area, matrix):
            offset = 0
            for start, length in zip(
                row["F_CHAN"][: row["N_GRP"]], row["N_CHAN"][: row["N_GRP"]]
            ):
                counts[start : start + length] += (
                    scale * row["MATRIX"][offset : offset + length]
                )
                offset += length
    return counts


def main():
    import xspec

    cache = HERE / "cache"
    cache.mkdir(exist_ok=True)
    archive = cache / URL.rsplit("/", 1)[1]
    if not archive.exists():
        urlretrieve(URL, archive)
    paths = [cache / f"pn-thin-5.{suffix}" for suffix in ("rmf", "arf")]
    with tarfile.open(archive) as bundle:
        for path in paths:
            path.write_bytes(
                bundle.extractfile(f"xmm_ao26_calfiles/{path.name}").read()
            )
    for name, url in [
        (SINGLE_RMF, CCF_URL + "extras/responses/PN/" + SINGLE_RMF),
        (QE_FILE, CCF_URL + QE_FILE),
    ]:
        path = cache / name
        if not path.exists():
            urlretrieve(url, path)
        paths.append(path)
    xspec.Xset.chatter = xspec.Xset.logChatter = 0
    xspec.Xset.abund, xspec.Xset.xsect = "angr", "bcmc"
    xspec.AllData.clear()
    xspec.AllModels.clear()
    model = xspec.Model("phabs*powerlaw", setPars={1: 0.03, 2: 1.7, 3: 1.0})
    xspec.AllData.fakeit(
        1,
        xspec.FakeitSettings(response=str(paths[0]), arf=str(paths[1]), exposure=1),
        applyStats=False,
        noWrite=True,
    )
    counts = {
        4: np.asarray(model.folded(1)),
        0: fold_singles(np.asarray(model.values(1)), cache),
    }
    channels = fits.getdata(paths[0], "EBOUNDS")["CHANNEL"]

    parameters = json.loads((HERE / "parameters.json").read_text())
    requests = parameters["bands_eV"] + [SHARED["soft_band_ev"]]
    bands = {}
    for low, high in requests:
        xspec.AllModels.calcFlux(f"{low / 1000} {high / 1000}")
        energy_flux = xspec.AllData(1).flux[0]
        rate = float(
            counts[pattern_max((low, high))] @ channel_weights(channels, low, high)
        )
        bands[f"{low}-{high}"] = {
            "pattern_max": pattern_max((low, high)),
            "model_rate_ct_s": rate,
            "model_energy_flux_erg_cm2_s": energy_flux,
            "ecf_counts_cm2_per_erg": rate / energy_flux,
        }
    result = {
        "method": "predicted band count rate / model band energy flux, directly from standard responses",
        "model": "phabs*powerlaw",
        "photon_index": 1.7,
        "nh_cm-2": 3e20,
        "abundances": "angr",
        "cross_sections": "bcmc",
        "flux_convention": "absorbed; no deabsorption correction",
        "response_url": URL,
        "singles_rmf_url": CCF_URL + "extras/responses/PN/" + SINGLE_RMF,
        "pattern_fraction_url": CCF_URL + QE_FILE,
        "singles_response": "epoch e4, full frame, Y9; standard ARF multiplied by CCF FRAC_S/FRAC_SD at the same position",
        "response_files": [
            {"name": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for path in paths
        ],
        "xspec_version": xspec.Xset.version,
        "approximation": "standard point-source responses, including their aperture losses, used directly for preliminary luminosities",
        "bands": bands,
    }
    (HERE / "ecfs.json").write_text(json.dumps(result, indent=2) + "\n")
    for name, row in bands.items():
        print(f"{name} eV: K={row['ecf_counts_cm2_per_erg']:.6g} counts cm2/erg")


if __name__ == "__main__":
    main()
