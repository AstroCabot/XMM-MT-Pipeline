#!/usr/bin/env python3
"""End-to-end QC runner for the comet pipeline.

Replaces the old xmm_comet_spotcheck_frames.sh and xmm_qc_build_mosaics.sh
shell helpers with a single Python entry-point.  Orchestrates:
  1. Stage manifest (xmm_qc_manifest.py)
  2. Quick comet-frame diagnostic (xmm_comet_quick_diag.py)
  3. Full QC checks (xmm_comet_quick_checks.py)
  4. QC mosaic building (SAS evselect + eexpmap, then combine + postproc)
  5. Segment spot-checks
  6. Final HTML/Markdown index

Usage:  xmm_qc_runner.py <config.env>
"""
from __future__ import annotations
import argparse, json, os, shlex, shutil, subprocess, sys
from pathlib import Path


def load_env(*sources: str) -> dict[str, str]:
    src = " ".join(
        f"source {shlex.quote(s)} >/dev/null 2>&1;"
        for s in sources
        if s and Path(s).is_file()
    )
    cmd = ["bash", "-lc", f"set -a; {src} env -0"]
    proc = subprocess.run(cmd, capture_output=True, check=True)
    env: dict[str, str] = {}
    for chunk in proc.stdout.split(b"\x00"):
        if chunk and b"=" in chunk:
            k, v = chunk.split(b"=", 1)
            env[k.decode()] = v.decode()
    return env


def run(cmd: list[str], env: dict[str, str] | None = None) -> None:
    subprocess.run(cmd, check=True, env=env)


def need(*cmds: str) -> None:
    missing = [cmd for cmd in cmds if shutil.which(cmd) is None]
    if missing:
        raise RuntimeError(f"Missing required command(s): {', '.join(missing)}")


def parse_bands(spec: str) -> list[tuple[str, int, int]]:
    out: list[tuple[str, int, int]] = []
    for entry in spec.split(";"):
        if not entry.strip():
            continue
        label, pimin, pimax = [x.strip() for x in entry.split(":", 2)]
        out.append((label, int(pimin), int(pimax)))
    return out


def truthy(value: str | None) -> bool:
    return str(value or "").strip().lower() in {"1", "y", "yes", "true"}


def skip(path: Path, force: bool) -> bool:
    return path.exists() and path.stat().st_size > 0 and not force


def symlink_force(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    dst.symlink_to(src.resolve())


def first_existing(*paths: Path) -> Path | None:
    return next((p for p in paths if p.is_file()), None)


def grid_from_json(path: Path) -> dict[str, float | int]:
    return json.loads(path.read_text(encoding="utf-8"))


def sas_env(cfg: dict[str, str], workdir: Path) -> dict[str, str]:
    env = os.environ.copy()
    env.update(cfg)
    env["SAS_CCF"] = str(workdir / "init" / "ccf.cif")
    sumsas = next(iter(sorted((workdir / "init").glob("*SUM.SAS"))), None)
    if sumsas is not None:
        env["SAS_ODF"] = str(sumsas)
    return env


def ensure_grid(
    cfg: dict[str, str], workdir: Path, scripts: Path
) -> dict[str, float | int]:
    grid_json = workdir / "images" / "grid.json"
    if not grid_json.is_file():
        ref_evt = first_existing(
            workdir / "comet" / "PN_comet.fits",
            workdir / "comet" / "M1_comet.fits",
            workdir / "comet" / "M2_comet.fits",
        )
        if ref_evt is None:
            raise RuntimeError("No comet-frame event list found")
        run(
            [
                sys.executable,
                str(scripts / "xmm_event_geometry.py"),
                "image-grid",
                "--event",
                str(ref_evt),
                "--radius-arcsec",
                cfg.get("IMAGE_RADIUS_ARCSEC", "300.0"),
                "--bin-arcsec",
                cfg.get("IMAGE_BIN_ARCSEC", "2.0"),
                "--out-json",
                str(grid_json),
            ]
        )
    return grid_from_json(grid_json)


def qc_band(cfg: dict[str, str], requested: str | None) -> tuple[str, int, int]:
    bands = parse_bands(cfg.get("IMAGE_BANDS", "broad:300:2000"))
    if requested:
        for band in bands:
            if band[0] == requested:
                return band
        raise RuntimeError(
            f"Requested band label {requested!r} not found in IMAGE_BANDS"
        )
    return bands[0]


def segment_rows(
    csv_path: Path, manual: str | None
) -> list[tuple[int, str, str, str, str]]:
    import csv

    rows = list(csv.DictReader(csv_path.open("r", newline="", encoding="utf-8")))
    if not rows:
        raise RuntimeError("motion_segments.csv is empty")
    if manual and manual.strip():
        wanted = {int(x.strip()) for x in manual.split(",") if x.strip()}
        keep = [r for r in rows if int(float(r["segid"])) in wanted]
    else:
        idxs = sorted({0, len(rows) // 2, len(rows) - 1})
        keep = [rows[i] for i in idxs]
    return [
        (
            int(float(r["segid"])),
            r["tstart"],
            r["tstop"],
            r.get("utc_start", ""),
            r.get("utc_stop", ""),
        )
        for r in keep
    ]


def build_mosaics(cfg: dict[str, str], scripts: Path, force: bool) -> None:
    need("evselect", "eexpmap")
    workdir = Path(cfg["WORKDIR"])
    env = sas_env(load_env(cfg["CONFIG_FILE"], str(workdir / "sas_setup.env")), workdir)
    if not env.get("SAS_CCF") or not env.get("SAS_ODF"):
        raise RuntimeError("Missing SAS setup")
    grid = grid_from_json(workdir / "images" / "grid.json")
    orig_atthk = workdir / "repro" / "atthk.dat"
    if not orig_atthk.is_file():
        raise RuntimeError(f"Missing original attitude: {orig_atthk}")
    stillsky_root = workdir / "qc" / "mosaics" / "stillsky"
    comet_root = workdir / "qc" / "mosaics" / "comet"
    only_insts = os.environ.get("ONLY_INSTS", "PN M1 M2").split()
    for label, pimin, pimax in parse_bands(
        cfg.get("DETECT_QC_BANDS", "soft:200:1000;hard:1000:12000")
    ):
        banddir = stillsky_root / label
        banddir.mkdir(parents=True, exist_ok=True)
        counts_list: list[str] = []
        expos_list: list[str] = []
        for inst in only_insts:
            evt = workdir / "clean" / f"{inst}_clean_merged.fits"
            if not evt.is_file():
                continue
            img = banddir / f"{inst}_counts.fits"
            exp = banddir / f"{inst}_exp.fits"
            if not skip(img, force):
                run(
                    [
                        "evselect",
                        f"table={evt}:EVENTS",
                        "withimageset=yes",
                        f"imageset={img}",
                        "xcolumn=X",
                        "ycolumn=Y",
                        "imagebinning=binSize",
                        f'ximagebinsize={grid["bin_phys"]}',
                        f'yimagebinsize={grid["bin_phys"]}',
                        "withxranges=yes",
                        f'ximagemin={grid["x_min_phys"]}',
                        f'ximagemax={grid["x_max_phys"]}',
                        "withyranges=yes",
                        f'yimagemin={grid["y_min_phys"]}',
                        f'yimagemax={grid["y_max_phys"]}',
                        "writedss=yes",
                        f"expression=PI in [{pimin}:{pimax}]",
                    ],
                    env,
                )
            if not skip(exp, force):
                run(
                    [
                        "eexpmap",
                        f"imageset={img}",
                        f"attitudeset={orig_atthk}",
                        f"eventset={evt}",
                        f"expimageset={exp}",
                        f"pimin={pimin}",
                        f"pimax={pimax}",
                        f'attrebin={cfg.get("QC_EEXPMAP_ATTREBIN", cfg.get("SCIENCE_EEXPMAP_ATTREBIN", "0.020626481"))}',
                    ],
                    env,
                )
            counts_list.append(str(img))
            expos_list.append(str(exp))
        if counts_list and not skip(banddir / "EPIC_counts.fits", force):
            run(
                [
                    sys.executable,
                    str(scripts / "combine_epic_images.py"),
                    "--counts",
                    *counts_list,
                    "--exposure",
                    *expos_list,
                    "--out-counts",
                    str(banddir / "EPIC_counts.fits"),
                    "--out-exposure",
                    str(banddir / "EPIC_exp.fits"),
                    "--out-rate",
                    str(banddir / "EPIC_rate.fits"),
                    "--out-json",
                    str(banddir / "EPIC_summary.json"),
                ]
            )
    mask_srclist = first_existing(
        workdir / "detect" / "field_sources_all.fits",
        workdir / "detect" / "field_sources_curated.fits",
    )
    for label, pimin, pimax in parse_bands(
        cfg.get("IMAGE_BANDS", "soft:200:1000;broad:300:2000;hard:1000:12000")
    ):
        src_banddir = workdir / "images" / label
        dstdir = comet_root / label
        dstdir.mkdir(parents=True, exist_ok=True)
        counts_list = []
        expos_list = []
        for inst in only_insts:
            cimg = src_banddir / f"{inst}_counts.fits"
            eimg = src_banddir / f"{inst}_exp.fits"
            if not (cimg.is_file() and eimg.is_file()):
                continue
            csrc, esrc = dstdir / f"{inst}_counts.fits", dstdir / f"{inst}_exp.fits"
            symlink_force(cimg, csrc)
            symlink_force(eimg, esrc)
            counts_list.append(str(csrc))
            expos_list.append(str(esrc))
        out_rate = dstdir / "EPIC_rate.fits"
        if counts_list and not skip(out_rate, force):
            cmd = [
                sys.executable,
                str(scripts / "combine_epic_images.py"),
                "--counts",
                *counts_list,
                "--exposure",
                *expos_list,
                "--out-counts",
                str(dstdir / "EPIC_counts.fits"),
                "--out-exposure",
                str(dstdir / "EPIC_exp.fits"),
                "--out-rate",
                str(out_rate),
                "--out-json",
                str(dstdir / "EPIC_summary.json"),
            ]
            if (
                mask_srclist
                and (workdir / "images" / "grid.json").is_file()
                and (workdir / "track" / "comet_track.fits").is_file()
            ):
                cmd += [
                    "--grid-json",
                    str(workdir / "images" / "grid.json"),
                    "--track",
                    str(workdir / "track" / "comet_track.fits"),
                    "--srclist",
                    str(mask_srclist),
                    "--mask-radius-arcsec",
                    cfg.get(
                        "IMAGE_MASK_RADIUS_ARCSEC",
                        cfg.get("FIELD_SOURCE_MASK_R_ARCSEC", "20.0"),
                    ),
                    "--out-mask",
                    str(dstdir / "EPIC_trail_mask.fits"),
                    "--out-rate-masked",
                    str(dstdir / "EPIC_rate_masked.fits"),
                ]
            run(cmd)
    run(
        [
            sys.executable,
            str(scripts / "image_postproc.py"),
            "qc",
            "--config",
            str(Path(cfg["CONFIG_FILE"])),
            "--stillsky-root",
            str(stillsky_root),
            "--comet-root",
            str(comet_root),
            "--detect-srclist",
            str(mask_srclist or ""),
            "--grid-json",
            str(workdir / "images" / "grid.json"),
            "--track",
            str(workdir / "track" / "comet_track.fits"),
            "--mask-radius-arcsec",
            cfg.get(
                "IMAGE_MASK_RADIUS_ARCSEC",
                cfg.get("FIELD_SOURCE_MASK_R_ARCSEC", "20.0"),
            ),
        ]
    )


def build_spotcheck(
    cfg: dict[str, str], scripts: Path, force: bool, band_label: str | None
) -> str:
    need("evselect", "eexpmap")
    workdir = Path(cfg["WORKDIR"])
    env = sas_env(load_env(cfg["CONFIG_FILE"], str(workdir / "sas_setup.env")), workdir)
    moved_atthk = workdir / "comet" / "moved_atthk.dat"
    seg_csv = workdir / "track" / "motion_segments.csv"
    if not moved_atthk.is_file():
        raise RuntimeError(f"Missing {moved_atthk}; run pipeline comet stage first")
    if not seg_csv.is_file():
        raise RuntimeError(f"Missing {seg_csv}; run pipeline track stage first")
    label, pimin, pimax = qc_band(cfg, band_label)
    grid = ensure_grid(cfg, workdir, scripts)
    outroot = workdir / "qc" / "spotchecks" / label
    outroot.mkdir(parents=True, exist_ok=True)
    only_insts = os.environ.get("ONLY_INSTS", "PN M1 M2").split()
    for segid, tstart, tstop, utc_start, utc_stop in segment_rows(
        seg_csv, cfg.get("SPOTCHECK_SEGMENTS")
    ):
        segdir = outroot / f"seg{segid}"
        segdir.mkdir(parents=True, exist_ok=True)
        meta = segdir / "segment_info.txt"
        if force or not meta.is_file():
            meta.write_text(
                f"segid={segid}\ntstart={tstart}\ntstop={tstop}\nutc_start={utc_start}\nutc_stop={utc_stop}\nband_label={label}\npimin={pimin}\npimax={pimax}\n",
                encoding="utf-8",
            )
        counts_list: list[str] = []
        expos_list: list[str] = []
        for inst in only_insts:
            evt = workdir / "comet" / f"{inst}_comet.fits"
            if not evt.is_file():
                continue
            img = segdir / f"{inst}_counts.fits"
            exp = segdir / f"{inst}_exp.fits"
            if not skip(img, force):
                run(
                    [
                        "evselect",
                        f"table={evt}:EVENTS",
                        "withimageset=yes",
                        f"imageset={img}",
                        "xcolumn=X",
                        "ycolumn=Y",
                        "imagebinning=binSize",
                        f'ximagebinsize={grid["bin_phys"]}',
                        f'yimagebinsize={grid["bin_phys"]}',
                        "withxranges=yes",
                        f'ximagemin={grid["x_min_phys"]}',
                        f'ximagemax={grid["x_max_phys"]}',
                        "withyranges=yes",
                        f'yimagemin={grid["y_min_phys"]}',
                        f'yimagemax={grid["y_max_phys"]}',
                        "writedss=yes",
                        f"expression=(TIME in [{tstart}:{tstop}])&&(PI in [{pimin}:{pimax}])",
                    ],
                    env,
                )
            if not skip(exp, force):
                run(
                    [
                        "eexpmap",
                        f"imageset={img}",
                        f"attitudeset={moved_atthk}",
                        f"eventset={evt}",
                        f"expimageset={exp}",
                        f"pimin={pimin}",
                        f"pimax={pimax}",
                        f'attrebin={cfg.get("QC_EEXPMAP_ATTREBIN", cfg.get("SCIENCE_EEXPMAP_ATTREBIN", "0.020626481"))}',
                    ],
                    env,
                )
            counts_list.append(str(img))
            expos_list.append(str(exp))
        if not counts_list:
            continue
        mask_srclist = first_existing(
            workdir / "detect" / "field_sources_all.fits",
            workdir / "detect" / "field_sources_curated.fits",
        )
        cmd = [
            sys.executable,
            str(scripts / "combine_epic_images.py"),
            "--counts",
            *counts_list,
            "--exposure",
            *expos_list,
            "--out-counts",
            str(segdir / "EPIC_counts.fits"),
            "--out-exposure",
            str(segdir / "EPIC_exp.fits"),
            "--out-rate",
            str(segdir / "EPIC_rate.fits"),
        ]
        if mask_srclist:
            cmd += [
                "--grid-json",
                str(workdir / "images" / "grid.json"),
                "--track",
                str(workdir / "track" / "comet_track.fits"),
                "--track-tmin-sec",
                str(tstart),
                "--track-tmax-sec",
                str(tstop),
                "--srclist",
                str(mask_srclist),
                "--mask-radius-arcsec",
                cfg.get(
                    "IMAGE_MASK_RADIUS_ARCSEC",
                    cfg.get("FIELD_SOURCE_MASK_R_ARCSEC", "20.0"),
                ),
                "--out-mask",
                str(segdir / "EPIC_trail_mask.fits"),
                "--out-rate-masked",
                str(segdir / "EPIC_rate_masked.fits"),
                "--out-json",
                str(segdir / "EPIC_image_summary.json"),
            ]
        if not skip(segdir / "EPIC_rate.fits", force):
            run(cmd)
    return label


def qc_run(config_file: str) -> None:
    scripts = Path(__file__).resolve().parent
    cfg = load_env(config_file)
    cfg["CONFIG_FILE"] = str(Path(config_file).resolve())
    workdir = Path(cfg["WORKDIR"])
    (workdir / "qc").mkdir(parents=True, exist_ok=True)
    band = qc_band(cfg, cfg.get("QC_SPOTCHECK_BAND") or None)[0]
    horizons = ["--horizons"] if truthy(cfg.get("QC_COMPARE_HORIZONS", "no")) else []
    force = truthy(os.environ.get("FORCE", "0"))
    first_label, first_pimin, first_pimax = parse_bands(
        cfg.get("IMAGE_BANDS", "broad:300:2000")
    )[0]
    run(
        [
            sys.executable,
            str(scripts / "xmm_qc_manifest.py"),
            "manifest",
            "--config",
            config_file,
            "--outdir",
            str(workdir / "qc"),
        ]
    )
    if list((workdir / "comet").glob("*_comet.fits")):
        try:
            run(
                [
                    sys.executable,
                    str(scripts / "xmm_comet_quick_diag.py"),
                    str(workdir),
                    "--pimin",
                    str(first_pimin),
                    "--pimax",
                    str(first_pimax),
                    "--nbins",
                    cfg.get("QC_QUICKDIAG_NBINS", "10"),
                    "--bin-arcsec",
                    cfg.get("QC_QUICKDIAG_BIN_ARCSEC", "8.0"),
                    "--radius-arcsec",
                    cfg.get("IMAGE_RADIUS_ARCSEC", "900"),
                    "--src-r-arcsec",
                    cfg.get("SRC_R_ARCSEC", "500"),
                ]
            )
        except subprocess.CalledProcessError:
            pass
    run(
        [
            sys.executable,
            str(scripts / "xmm_comet_quick_checks.py"),
            "--config",
            config_file,
            "--outdir",
            str(workdir / "qc"),
            *horizons,
            "clean",
            "exposures",
            "track",
            "detect",
            "contam",
            "image",
            "lcurve",
            "spectrum",
        ]
    )
    for inst in ("EPIC", "PN", "M1", "M2"):
        try:
            run(
                [
                    sys.executable,
                    str(scripts / "xmm_comet_quick_checks.py"),
                    "--config",
                    config_file,
                    "--outdir",
                    str(workdir / "qc"),
                    "--inst",
                    inst,
                    "region-support",
                ]
            )
        except subprocess.CalledProcessError:
            pass
    try:
        build_mosaics(cfg, scripts, force)
        run(
            [
                sys.executable,
                str(scripts / "xmm_comet_quick_checks.py"),
                "--config",
                config_file,
                "--outdir",
                str(workdir / "qc"),
                "mosaics",
            ]
        )
    except Exception:
        pass
    if (workdir / "track" / "motion_segments.csv").is_file() and (
        workdir / "comet" / "moved_atthk.dat"
    ).is_file():
        try:
            band = build_spotcheck(cfg, scripts, force, band)
            run(
                [
                    sys.executable,
                    str(scripts / "xmm_comet_quick_checks.py"),
                    "--config",
                    config_file,
                    "--outdir",
                    str(workdir / "qc"),
                    "--band",
                    band,
                    "spotcheck",
                ]
            )
        except Exception:
            pass
    run(
        [
            sys.executable,
            str(scripts / "xmm_qc_manifest.py"),
            "index",
            "--config",
            config_file,
            "--outdir",
            str(workdir / "qc"),
        ]
    )
    print(f'QC outputs written under {workdir / "qc"}')
    print("Start here:")
    for path in (
        "STAGE_QC.md",
        "00_manifest.png",
        "00_clean.png",
        "00b_exposures.png",
        "01_track.png",
        "02_detect.png",
        "03_contam.png",
        "04_image.png",
        "05_lcurve.png",
        "06_spectrum.png",
        f"07_spotchecks_{band}.png",
        "08_mosaics.png",
    ):
        print(f'  {workdir / "qc" / path}')
    print(f'  {workdir / "qc" / "quick_diag" / "mosaic.png"}')
    print(f'  {workdir / "qc" / "quick_diag" / "comparison.png"}')


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    args = ap.parse_args(argv)
    qc_run(args.config)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
