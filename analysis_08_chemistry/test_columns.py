import numpy as np

import run
from columns import FITS, GROUPS, load, shares

REFERENCE = {
    "parents_shallow": {
        "number": (0.0, 0.1225161336, 0.7127088157),
        "sigma_weighted": (0.0, 0.1213488901, 0.7059186506),
    },
    "parents_steep": {
        "number": (0.0, 0.2870055757, 0.5021255903),
        "sigma_weighted": (0.0, 0.2840938753, 0.4970314757),
    },
    "radiolytic_shallow": {
        "number": (0.6984476103, 0.0765687877, 0.1348567145),
        "sigma_weighted": (0.6669219352, 0.0829415518, 0.1460807401),
    },
    "radiolytic_steep": {
        "number": (0.4978551412, 0.0903968305, 0.3043494334),
        "sigma_weighted": (0.4641934147, 0.0956155223, 0.3219198051),
    },
}

inputs = run.load_inputs()
operator = run.load_psf(run.OUTPUT / "psf.npz", inputs, run.psf_dependencies())
for name, (folder, *model) in FITS.items():
    posterior, area, best, samples = load(inputs, operator, folder, *model)
    for theta in (*samples[::12000], best):
        result = shares(theta, posterior, area)
        for values in result.values():
            assert tuple(values) == GROUPS
            assert abs(sum(values.values()) - 1) < 1e-12, (name, values)
    for label, expected in REFERENCE[name].items():
        found = [result[label][group] for group in GROUPS[:3]]
        assert np.allclose(found, expected, rtol=0, atol=1e-9), (name, label, found)
