#!/usr/bin/env python3
from __future__ import annotations
import argparse
import json
import math
import shlex
import subprocess
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

def load_shell_env(path: str) -> dict[str, str]:
    quoted = shlex.quote(str(path))
    cmd = ['bash', '-lc', f'set -a; source {quoted} >/dev/null 2>&1; env -0']
    proc = subprocess.run(cmd, capture_output=True, check=True)
    env: dict[str, str] = {}
    for chunk in proc.stdout.split(b'\x00'):
        if not chunk or b'=' not in chunk:
            continue
        k, v = chunk.split(b'=', 1)
        env[k.decode()] = v.decode()
    return env

def read_list(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [line.strip() for line in path.read_text(encoding='utf-8').splitlines() if line.strip()]

def fits_rows(path: Path) -> int | None:
    if not path.exists():
        return None
    try:
        with fits.open(path) as hdul:
            if 'EVENTS' in hdul:
                data = hdul['EVENTS'].data
                return 0 if data is None else int(len(data))
            for hdu in hdul[1:]:
                if isinstance(hdu, fits.BinTableHDU) and hdu.data is not None:
                    return int(len(hdu.data))
    except Exception:
        return None
    return None

def nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0

def stage_status(ok: bool, warn: bool=False) -> str:
    if ok and (not warn):
        return 'ok'
    if ok and warn:
        return 'warn'
    return 'missing'

def cmd_manifest(args: argparse.Namespace) -> int:
    env = load_shell_env(args.config)
    workdir = Path(env['WORKDIR'])
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stages: dict[str, dict[str, object]] = {}
    init_ok = nonempty(workdir / 'sas_setup.env') and nonempty(workdir / 'init' / 'ccf.cif')
    stages['init'] = {'status': stage_status(init_ok), 'sas_setup_env': str(workdir / 'sas_setup.env'), 'ccf_cif': str(workdir / 'init' / 'ccf.cif')}
    repro_data: dict[str, dict[str, object]] = {}
    raw_counts = []
    for inst in ('PN', 'M1', 'M2'):
        files = read_list(workdir / 'repro' / 'manifest' / f'{inst}_raw.txt')
        raw_counts.append(len(files))
        repro_data[inst] = {'n_raw_eventlists': len(files), 'example': files[:3]}
    attitude_ok = nonempty(workdir / 'repro' / 'atthk.dat') and nonempty(workdir / 'attitude.env')
    stages['repro'] = {'status': stage_status(any((c > 0 for c in raw_counts)) and attitude_ok), 'attitude_ok': attitude_ok, 'per_instrument': repro_data}
    clean_data: dict[str, dict[str, object]] = {}
    clean_counts = []
    merged_rows = []
    for inst in ('PN', 'M1', 'M2'):
        files = read_list(workdir / 'clean' / f'{inst}_clean_files.txt')
        merged = workdir / 'clean' / f'{inst}_clean_merged.fits'
        rows = fits_rows(merged)
        clean_counts.append(len(files))
        merged_rows.append(rows if rows is not None else math.nan)
        clean_data[inst] = {'n_clean_eventlists': len(files), 'merged_eventlist': str(merged), 'merged_rows': rows}
    clean_warn = any((raw_counts[i] != clean_counts[i] and raw_counts[i] > 0 for i in range(3)))
    stages['clean'] = {'status': stage_status(any((c > 0 for c in clean_counts)), warn=clean_warn), 'per_instrument': clean_data}
    track_fits = workdir / 'track' / 'comet_track.fits'
    seg_csv = workdir / 'track' / 'motion_segments.csv'
    n_segments = 0
    if seg_csv.exists():
        n_segments = max(0, len(seg_csv.read_text(encoding='utf-8').splitlines()) - 1)
    stages['track'] = {'status': stage_status(nonempty(track_fits) and nonempty(seg_csv)), 'track_fits': str(track_fits), 'motion_segments_csv': str(seg_csv), 'n_motion_segments': n_segments}
    pseudoexps = read_list(workdir / 'detect' / 'stack_eventsets.txt')
    detected_rows = fits_rows(workdir / 'detect' / 'field_sources_all.fits')
    curated_rows = fits_rows(workdir / 'detect' / 'field_sources_curated.fits')
    detect_mosaics = [str(p) for p in sorted((workdir / 'detect').glob('EPIC_*_mosaic.fits'))]
    detect_ok = nonempty(workdir / 'detect' / 'stack_srclist.fits') and len(pseudoexps) > 0
    stages['detect'] = {'status': stage_status(detect_ok), 'n_pseudoexposures': len(pseudoexps), 'n_detected_sources': detected_rows, 'n_curated_sources': curated_rows, 'mosaic_images': detect_mosaics}
    comet_ok = nonempty(workdir / 'comet' / 'moved_sas_setup.env')
    comet_data: dict[str, object] = {'status': stage_status(comet_ok), 'moved_sas_setup_env': str(workdir / 'comet' / 'moved_sas_setup.env')}
    for inst in ('PN', 'M1', 'M2'):
        rows = fits_rows(workdir / 'comet' / f'{inst}_comet.fits')
        comet_data[inst] = {'rows': rows, 'path': str(workdir / 'comet' / f'{inst}_comet.fits')}
    stages['comet'] = comet_data
    contam_csv = workdir / 'contam' / 'contamination_report.csv'
    timeline_csv = workdir / 'contam' / 'contamination_timeline.csv'
    gti = workdir / 'contam' / 'science_gti.fits'
    contam_rows = 0
    if contam_csv.exists():
        contam_rows = max(0, len(contam_csv.read_text(encoding='utf-8').splitlines()) - 1)
    summary_json = workdir / 'contam' / 'contamination_summary.json'
    summary = load_json(summary_json)
    stages['contam'] = {'status': stage_status(nonempty(gti)), 'science_gti': str(gti), 'science_policy': summary.get('science_policy'), 'timeline_csv': str(timeline_csv), 'n_contam_sources': contam_rows}
    image_bands = []
    for entry in env.get('IMAGE_BANDS', 'broad:300:2000').split(';'):
        entry = entry.strip()
        if not entry:
            continue
        image_bands.append(entry.split(':')[0].strip())
    image_products = {}
    image_ok = False
    for band in image_bands:
        rate = workdir / 'images' / band / 'EPIC_rate.fits'
        exp = workdir / 'images' / band / 'EPIC_exp.fits'
        masked = workdir / 'images' / band / 'EPIC_rate_masked.fits'
        ok = nonempty(rate) and nonempty(exp)
        image_ok = image_ok or ok
        image_products[band] = {'rate': str(rate), 'exp': str(exp), 'rate_masked': str(masked), 'ok': ok}
    stages['image'] = {'status': stage_status(image_ok), 'bands': image_products}
    lc_summary = workdir / 'lcurve' / 'lc_summary.tsv'
    mode_file = workdir / 'lcurve' / 'EPIC_total_corr_mode.txt'
    combined_mode = None
    if mode_file.exists():
        for line in mode_file.read_text(encoding='utf-8').splitlines():
            if line.startswith('mode='):
                combined_mode = line.split('=', 1)[1].strip()
    req_abs = env.get('LC_APPLY_ABSOLUTE_CORRECTIONS', 'no').strip().lower() in {'1', 'y', 'yes', 'true'}
    lcurve_warn = combined_mode == 'raw' or (req_abs and combined_mode not in {'abs'})
    stages['lcurve'] = {'status': stage_status(nonempty(workdir / 'lcurve' / 'EPIC_total_corr_lc.fits'), warn=lcurve_warn), 'combined_mode': combined_mode, 'lc_summary_tsv': str(lc_summary)}
    spec_ok = nonempty(workdir / 'spectra' / 'EPIC_src_combined.fits') or nonempty(workdir / 'spectra' / 'EPIC_src_combined_grp.fits')
    stages['spectrum'] = {'status': stage_status(spec_ok), 'combined_src': str(workdir / 'spectra' / 'EPIC_src_combined.fits'), 'combined_grp': str(workdir / 'spectra' / 'EPIC_src_combined_grp.fits')}
    merge_ok = nonempty(workdir / 'final' / 'EPIC_comet_merged_sciencegti.fits')
    stages['merge'] = {'status': stage_status(merge_ok), 'science_eventlist': str(workdir / 'final' / 'EPIC_comet_merged_sciencegti.fits'), 'full_eventlist': str(workdir / 'final' / 'EPIC_comet_merged_full.fits')}
    payload = {'workdir': str(workdir), 'stages': stages}
    (outdir / '00_manifest.json').write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    md_lines = ['# Stage manifest summary', '', f'Workdir: `{workdir}`', '', '| Stage | Status | Key check |', '|---|---|---|']
    for name in ['init', 'repro', 'clean', 'track', 'detect', 'comet', 'contam', 'image', 'lcurve', 'spectrum', 'merge']:
        st = stages[name]
        status = str(st.get('status', 'missing'))
        if name == 'repro':
            key = ', '.join((f"{inst}:{st['per_instrument'][inst]['n_raw_eventlists']} raw" for inst in ('PN', 'M1', 'M2')))
        elif name == 'clean':
            key = ', '.join((f"{inst}:{st['per_instrument'][inst]['n_clean_eventlists']} clean" for inst in ('PN', 'M1', 'M2')))
        elif name == 'detect':
            key = f"{st.get('n_pseudoexposures', 0)} pseudoexps, {st.get('n_detected_sources', 'n/a')} detected, {st.get('n_curated_sources', 'n/a')} curated"
        elif name == 'lcurve':
            key = f"combined mode={st.get('combined_mode', 'n/a')}"
        elif name == 'track':
            key = f"{st.get('n_motion_segments', 0)} motion segments"
        else:
            key = 'outputs present' if status != 'missing' else 'missing'
        md_lines.append(f'| {name} | {status} | {key} |')
    (outdir / '00_manifest.md').write_text('\n'.join(md_lines) + '\n', encoding='utf-8')
    labels = ['PN raw', 'M1 raw', 'M2 raw', 'PN clean', 'M1 clean', 'M2 clean', 'pseudoexp', 'curated src']
    values = [raw_counts[0], raw_counts[1], raw_counts[2], clean_counts[0], clean_counts[1], clean_counts[2], len(pseudoexps), int(curated_rows or 0)]
    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    x = np.arange(len(labels))
    ax.bar(x, values)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha='right')
    ax.set_ylabel('count')
    ax.set_title('Pipeline manifest overview')
    for i, v in enumerate(values):
        ax.text(i, v + max(values + [1]) * 0.02, str(v), ha='center', va='bottom', fontsize=8)
    txt = [f"lcurve mode : {combined_mode or 'n/a'}", f'motion segs : {n_segments}', f"merged rows : PN={clean_data['PN']['merged_rows']} M1={clean_data['M1']['merged_rows']} M2={clean_data['M2']['merged_rows']}"]
    ax.text(0.99, 0.98, '\n'.join(txt), transform=ax.transAxes, ha='right', va='top', fontsize=8, family='monospace', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    fig.tight_layout()
    fig.savefig(outdir / '00_manifest.png', dpi=150)
    print(outdir / '00_manifest.png')
    return 0

def load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding='utf-8'))

def cmd_index(args: argparse.Namespace) -> int:
    env = load_shell_env(args.config)
    workdir = Path(env['WORKDIR'])
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = load_json(outdir / '00_manifest.json')
    stages = manifest.get('stages', {})
    report_md = outdir / 'qc_report.md'
    region_reports = {inst: outdir / f'region_support_{inst}.json' for inst in ('EPIC', 'PN', 'M1', 'M2')}
    lines = ['# Stage-by-stage QC index', '', f'Workdir: `{workdir}`', '', 'This file is intended to be the one place to inspect after a run. Each section points to the main image and/or text artifact for that stage.', '', '## Summary table', '', '| Stage | Status | Main artifact(s) |', '|---|---|---|']
    mapping = {'init': ['00_manifest.md'], 'repro': ['00_manifest.md', '00_manifest.png'], 'clean': ['00_clean.png', '00_manifest.md', '00_manifest.png'], 'track': ['01_track.png', 'qc_report.md'], 'detect': ['02_detect.png', '00_manifest.md'], 'comet': [f"07_spotchecks_{env.get('QC_SPOTCHECK_BAND', '') or env.get('IMAGE_BANDS', 'broad:300:2000').split(';')[0].split(':')[0]}.png", '00_manifest.md'], 'contam': ['04_contam.png', 'qc_report.md'], 'image': ['03_image.png', 'qc_report.md'], 'lcurve': ['05_lcurve.png', 'region_support_PN.json', 'qc_report.md'], 'spectrum': ['06_spectrum.png', 'qc_report.md'], 'merge': ['00_manifest.md']}
    for stage in ['init', 'repro', 'clean', 'track', 'detect', 'comet', 'contam', 'image', 'lcurve', 'spectrum', 'merge']:
        status = str(stages.get(stage, {}).get('status', 'missing'))
        arts = ', '.join((f'`{a}`' for a in mapping.get(stage, [])))
        lines.append(f'| {stage} | {status} | {arts} |')
    lines += ['', '## What to look for', '']
    lines += ['### init / repro / clean', 'Use `00_manifest.md` and `00_manifest.png` to confirm that all raw event lists were found, cleaned, and merged per instrument. For a multi-pointing observation, the raw counts should exceed one per instrument if the ODF really contains multiple pointings.', '', '### track', 'Inspect `01_track.png`. The track should be smooth, and the motion-segmentation curve should not show pathological jumps except possibly a small final-segment remainder effect.', '', '### detect', 'Inspect `02_detect.png`. The visible portion of the comet track should lie over the detect mosaics, and the soft/hard panels should show that source overlays track real sky structure rather than a clipped or shifted footprint. Compare the full detected-source list against the more conservative track-excluded list.', '', '### comet', 'Inspect the spot-check montage `07_spotchecks_*.png`. In a good comet frame the comet stays compact while stationary field sources trail. This is the single most important visual sanity check for the moving-target reprojection.', '', '### contam', 'Inspect `04_contam.png` and the contamination section in `qc_report.md`. The policy rows should make it obvious how much exposure would be lost by source-only, background-only, or strict contamination cuts, and the overlap-count panel should show whether masking is a better choice than discarding whole intervals.', '', '### image', 'Inspect `03_image.png`. The counts, exposure, rate, and masked-rate panels should all be finite in the comet-centered footprint; masked-rate images should suppress obvious trail contamination.', '', '### lcurve', 'Inspect `05_lcurve.png`, `lcurve/lc_summary.tsv`, and the `region_support_*.json` files. The strict default is that absolute-correction failures stop the pipeline unless you explicitly allow a rel-only fallback.', '', '### spectrum', 'Inspect `06_spectrum.png`. The combined source spectrum should sit above the scaled background over at least part of the configured band, and the approximate net panel should not be dominated by wild oscillations from a broken background selection.', '']
    lines.append('## Region-support quick notes')
    lines.append('')
    for inst, path in region_reports.items():
        if not path.exists():
            continue
        try:
            data = json.loads(path.read_text(encoding='utf-8'))
        except Exception:
            continue
        src_med = data.get('source', {}).get('median')
        bkg_med = data.get('background', {}).get('median')
        src_zero = data.get('source', {}).get('zero_fraction')
        bkg_zero = data.get('background', {}).get('zero_fraction')
        lines.append(f'- **{inst}**: source median exp={src_med}, background median exp={bkg_med}, source zero-fraction={src_zero}, background zero-fraction={bkg_zero}')
    lines.append('')
    if report_md.exists():
        lines.append('## Primary QC report')
        lines.append('')
        lines.append('The companion `qc_report.md` contains the detailed stage-level measurements produced by `xmm_comet_quick_checks.py`.')
        lines.append('')
    (outdir / 'STAGE_QC.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print(outdir / 'STAGE_QC.md')
    return 0

def find_candidates(repro: Path, patterns: list[str]) -> list[str]:
    out: list[str] = []
    for patt in patterns:
        out.extend((str(p) for p in sorted(repro.rglob(patt)) if p.is_file()))
    seen: set[str] = set()
    uniq: list[str] = []
    for p in out:
        if p not in seen:
            uniq.append(p)
            seen.add(p)
    return uniq

def cmd_repro(args: argparse.Namespace) -> int:
    env = load_shell_env(args.config)
    workdir = Path(env['WORKDIR'])
    repro = workdir / 'repro'
    data = {'PN': find_candidates(repro, ['*EPN*ImagingEvts.ds', '*EPN*ImagingEvts*.FIT*', '*PIEVLI*.FIT*', '*EPN*EVLI*.FIT*']), 'M1': find_candidates(repro, ['*EMOS1*ImagingEvts.ds', '*EMOS1*ImagingEvts*.FIT*', '*M1EVLI*.FIT*', '*EMOS1*EVLI*.FIT*']), 'M2': find_candidates(repro, ['*EMOS2*ImagingEvts.ds', '*EMOS2*ImagingEvts*.FIT*', '*M2EVLI*.FIT*', '*EMOS2*EVLI*.FIT*'])}
    print(json.dumps({k: {'count': len(v), 'files': v} for k, v in data.items()}, indent=2))
    return 0

def main() -> int:
    ap = argparse.ArgumentParser(description='Manifest, QC index, and repro diagnostics for the XMM comet package.')
    sub = ap.add_subparsers(dest='command')
    sub.required = True
    p_manifest = sub.add_parser('manifest', help='Write stage-level manifest JSON/MD/PNG')
    p_manifest.add_argument('--config', required=True)
    p_manifest.add_argument('--outdir', required=True)
    p_index = sub.add_parser('index', help='Build stage-by-stage QC navigation index')
    p_index.add_argument('--config', required=True)
    p_index.add_argument('--outdir', required=True)
    p_repro = sub.add_parser('repro', help='Report per-instrument calibrated event-list candidates')
    p_repro.add_argument('--config', required=True)
    args = ap.parse_args()
    if args.command == 'manifest':
        return cmd_manifest(args)
    elif args.command == 'index':
        return cmd_index(args)
    elif args.command == 'repro':
        return cmd_repro(args)
    return 1
if __name__ == '__main__':
    raise SystemExit(main())
