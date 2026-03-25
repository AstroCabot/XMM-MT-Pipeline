from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

try:  # pragma: no cover - exercised only when astropy is installed
    from astropy.io import fits
    from astropy.time import Time
    import astropy.units as u
    HAVE_ASTROPY = True
except Exception:  # pragma: no cover
    fits = None
    Time = None
    u = None
    HAVE_ASTROPY = False


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = PACKAGE_ROOT / "scripts"


def run_script(script: str, *args: str) -> subprocess.CompletedProcess[str]:
    cmd = [sys.executable, str(SCRIPTS / script), *map(str, args)]
    env = os.environ.copy()
    env.setdefault("PYTHONPATH", str(SCRIPTS))
    return subprocess.run(cmd, text=True, capture_output=True, env=env, check=False)


def write_track(path: Path, *, ra_start: float = 10.0, ra_stop: float = 10.01, dec: float = 20.0) -> None:
    mjdref = 59000.0
    times = np.asarray([0.0, 50.0, 100.0], dtype=float)
    mjd = (Time(mjdref, format="mjd", scale="utc") + times * u.s).mjd
    cols = fits.ColDefs(
        [
            fits.Column(name="MJD", format="D", array=mjd.astype(np.float64)),
            fits.Column(name="RA", format="D", array=np.asarray([ra_start, 0.5 * (ra_start + ra_stop), ra_stop], dtype=np.float64)),
            fits.Column(name="DEC", format="D", array=np.asarray([dec, dec, dec], dtype=np.float64)),
        ]
    )
    hdu = fits.BinTableHDU.from_columns(cols, name="TRACK")
    hdu.header["MJDREF"] = mjdref
    hdu.header["TIMESYS"] = "UTC"
    hdu.header["OBS_T0"] = 0.0
    hdu.header["OBS_T1"] = 100.0
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path, overwrite=True)


def write_srclist(path: Path, ras: list[float], decs: list[float]) -> None:
    n = len(ras)
    cols = fits.ColDefs(
        [
            fits.Column(name="SRCID", format="J", array=np.arange(1, n + 1, dtype=np.int32)),
            fits.Column(name="RA", format="D", array=np.asarray(ras, dtype=np.float64)),
            fits.Column(name="DEC", format="D", array=np.asarray(decs, dtype=np.float64)),
            fits.Column(name="DET_ML", format="E", array=np.full(n, 20.0, dtype=np.float32)),
            fits.Column(name="EXTENT", format="E", array=np.zeros(n, dtype=np.float32)),
        ]
    )
    hdu = fits.BinTableHDU.from_columns(cols, name="SRCLIST")
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path, overwrite=True)


def write_empty_detect_srclist(path: Path) -> None:
    cols = fits.ColDefs(
        [
            fits.Column(name="SRCID", format="J", array=np.asarray([], dtype=np.int32)),
            fits.Column(name="RA", format="D", array=np.asarray([], dtype=np.float64)),
            fits.Column(name="DEC", format="D", array=np.asarray([], dtype=np.float64)),
            fits.Column(name="DET_ML", format="E", array=np.asarray([], dtype=np.float32)),
            fits.Column(name="EXTENT", format="E", array=np.asarray([], dtype=np.float32)),
        ]
    )
    hdu = fits.BinTableHDU.from_columns(cols, name="SRCLIST")
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path, overwrite=True)


def gti_total(path: Path) -> float:
    with fits.open(path) as hdul:
        tab = hdul[1].data
        if tab is None or len(tab) == 0:
            return 0.0
        return float(np.sum(np.asarray(tab["STOP"], dtype=float) - np.asarray(tab["START"], dtype=float)))


@unittest.skipUnless(HAVE_ASTROPY, "astropy is required for helper integration tests")
class RefactorHelperTests(unittest.TestCase):
    def test_curate_emllist_handles_empty_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            src_in = tmpdir / "empty_stack_srclist.fits"
            src_out = tmpdir / "curated.fits"
            ds9 = tmpdir / "curated.reg"
            write_empty_detect_srclist(src_in)

            proc = run_script(
                "curate_emllist_for_region.py",
                "--input",
                src_in,
                "--output",
                src_out,
                "--ds9",
                ds9,
                "--min-det-ml",
                "10",
            )
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            self.assertTrue(src_out.exists())
            self.assertTrue(ds9.exists())
            with fits.open(src_out) as hdul:
                self.assertEqual(len(hdul[1].data), 0)

    def test_contamination_builder_handles_empty_source_list(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            track = tmpdir / "track.fits"
            srclist = tmpdir / "sources_empty.fits"
            outdir = tmpdir / "contam"
            write_track(track)
            write_srclist(srclist, [], [])

            proc = run_script(
                "build_contamination_products.py",
                "--track",
                track,
                "--srclist",
                srclist,
                "--output-dir",
                outdir,
                "--science-policy",
                "full",
                "--sample-step-s",
                "10",
                "--mask-radius-arcsec",
                "5",
                "--src-dx-arcsec",
                "0",
                "--src-dy-arcsec",
                "0",
                "--src-r-arcsec",
                "5",
                "--bkg-shape",
                "annulus",
                "--bkg-dx-arcsec",
                "0",
                "--bkg-dy-arcsec",
                "0",
                "--bkg-rin-arcsec",
                "20",
                "--bkg-rout-arcsec",
                "30",
            )
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            self.assertAlmostEqual(gti_total(outdir / "gti_full.fits"), 100.0, places=6)
            self.assertAlmostEqual(gti_total(outdir / "science_gti.fits"), 100.0, places=6)
            summary = json.loads((outdir / "contamination_summary.json").read_text())
            self.assertEqual(summary["n_sources"], 0)

    def test_policy_based_contamination_products_reduce_exposure_when_source_crosses_aperture(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            track = tmpdir / "track.fits"
            srclist = tmpdir / "sources.fits"
            outdir = tmpdir / "contam"
            write_track(track)
            # Middle source crosses the source aperture in the comet frame.
            write_srclist(srclist, [10.005], [20.0])

            proc = run_script(
                "build_contamination_products.py",
                "--track",
                track,
                "--srclist",
                srclist,
                "--output-dir",
                outdir,
                "--science-policy",
                "strict",
                "--sample-step-s",
                "5",
                "--mask-radius-arcsec",
                "5",
                "--src-dx-arcsec",
                "0",
                "--src-dy-arcsec",
                "0",
                "--src-r-arcsec",
                "5",
                "--bkg-shape",
                "annulus",
                "--bkg-dx-arcsec",
                "0",
                "--bkg-dy-arcsec",
                "0",
                "--bkg-rin-arcsec",
                "40",
                "--bkg-rout-arcsec",
                "50",
            )
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)

            full_time = gti_total(outdir / "gti_full.fits")
            src_time = gti_total(outdir / "gti_src.fits")
            bkg_time = gti_total(outdir / "gti_bkg.fits")
            strict_time = gti_total(outdir / "gti_strict.fits")
            science_time = gti_total(outdir / "science_gti.fits")

            self.assertGreater(full_time, strict_time)
            self.assertAlmostEqual(src_time, strict_time, places=6)
            self.assertAlmostEqual(bkg_time, full_time, places=6)
            self.assertAlmostEqual(science_time, strict_time, places=6)

    def test_trail_mask_and_event_filtering(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            track = tmpdir / "track.fits"
            srclist = tmpdir / "sources.fits"
            grid = tmpdir / "grid.json"
            mask = tmpdir / "trail_mask.fits"
            event = tmpdir / "events.fits"
            filtered = tmpdir / "events_filtered.fits"
            report = tmpdir / "events_filtered.json"

            write_track(track)
            write_srclist(srclist, [10.005], [20.0])
            grid.write_text(
                json.dumps(
                    {
                        "nx": 200,
                        "ny": 200,
                        "center_x_phys": 100.0,
                        "center_y_phys": 100.0,
                        "x_min_phys": 0.0,
                        "x_max_phys": 200.0,
                        "y_min_phys": 0.0,
                        "y_max_phys": 200.0,
                        "bin_phys": 1.0,
                        "x_arcsec_per_phys": 1.0,
                        "y_arcsec_per_phys": 1.0,
                        "scale_arcsec_per_phys": 1.0,
                    }
                )
            )

            proc = run_script(
                "make_trail_mask.py",
                "--grid-json",
                grid,
                "--track",
                track,
                "--srclist",
                srclist,
                "--mask-radius-arcsec",
                "4",
                "--out-mask",
                mask,
            )
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)

            # Three events on the trail, one safely away from it.
            cols = fits.ColDefs(
                [
                    fits.Column(name="X", format="E", array=np.asarray([83.5, 100.0, 116.5, 100.0], dtype=np.float32)),
                    fits.Column(name="Y", format="E", array=np.asarray([100.0, 100.0, 100.0, 140.0], dtype=np.float32)),
                    fits.Column(name="TIME", format="D", array=np.asarray([10.0, 50.0, 90.0, 50.0], dtype=np.float64)),
                ]
            )
            evt_hdu = fits.BinTableHDU.from_columns(cols, name="EVENTS")
            fits.HDUList([fits.PrimaryHDU(), evt_hdu]).writeto(event, overwrite=True)

            proc = run_script(
                "apply_spatial_trail_mask.py",
                "--event",
                event,
                "--mask",
                mask,
                "--grid-json",
                grid,
                "--output",
                filtered,
                "--report-json",
                report,
            )
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)

            with fits.open(filtered) as hdul:
                kept = hdul["EVENTS"].data
                self.assertEqual(len(kept), 1)
                self.assertAlmostEqual(float(kept["X"][0]), 100.0, places=4)
                self.assertAlmostEqual(float(kept["Y"][0]), 140.0, places=4)

            payload = json.loads(report.read_text())
            self.assertEqual(payload["removed_rows"], 3)
            self.assertEqual(payload["kept_rows"], 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
