import json
import math
from pathlib import Path

import numpy as np
from thermal import distance, simulate

HERE = Path(__file__).resolve().parent
MODEL = json.loads((HERE / "parameters.json").read_text())["model"]
MODEL = {**MODEL, "density_kg_m3": 1000.0, "depth_m": 10.0}
ORBIT = (1.357, 6.1437)
PEAKS_K = {
    1e-4: (254.0, 92.8),
    1e-3: (253.7, 169.9, 123.2),
    1e-2: (252.9, 220.4, 193.6),
}
DEPTH_150K_M = {1e-4: 0.041, 1e-3: 0.130, 1e-2: 0.400}
DEPTH_80K_M = {1e-3: (0.169, 0.267), 1e-2: (0.49, 0.81)}


def close(value, expected, absolute):
    assert math.isclose(value, expected, abs_tol=absolute), (value, expected)


def test_orbit():
    close(distance(*ORBIT, 0.0), 1.357, 1e-9)
    close(distance(*ORBIT, -120.0), 4.50, 0.01)
    close(distance(*ORBIT, 873.0), 30.0, 0.05)


def test_yaginuma():
    for k, peaks in PEAKS_K.items():
        out = simulate(k, MODEL, *ORBIT, 700.0, (80.0, 150.0), depths=(0.0, 0.1, 0.2))
        for value, expected in zip(out["T_K"].max(axis=0), peaks):
            close(value, expected, 1.5)
        close(out["deepest_m"][:, 1].max(), DEPTH_150K_M[k], 0.025)
        if k in DEPTH_80K_M:
            days = out["days"]
            for day, expected in zip((-120.0, 0.0), DEPTH_80K_M[k]):
                close(out["deepest_m"][np.searchsorted(days, day), 0], expected, 0.015)


def main():
    test_orbit()
    test_yaginuma()


if __name__ == "__main__":
    main()
