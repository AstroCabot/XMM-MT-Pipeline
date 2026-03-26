#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits

def event_count_fast(filepath):
    with fits.open(str(filepath)) as hdul:
        for ext in hdul:
            if ext.name == 'EVENTS':
                return ext.header.get('NAXIS2', 0)
    return 0

def main():
    workdir = Path('/mnt/c/Users/cabot/OneDrive/Documents/Research/Cambridge/X3I/reduction_v4/output')
    mt_dir = Path('/mnt/c/Users/cabot/OneDrive/Documents/Research/Data/xmm_3iatlas_updated/MT_PPS')
    outdir = Path('/mnt/c/Users/cabot/OneDrive/Documents/Research/Cambridge/X3I/reduction_v4/xmm_comet_pipeline_refactor_release/pps_analytics_output')
    outdir.mkdir(parents=True, exist_ok=True)
    lines = []
    lines.append('=' * 72)
    lines.append('PPS vs PIPELINE - Side-by-side comparison')
    lines.append('=' * 72)
    lines.append(f'MT_PPS:   {mt_dir}')
    lines.append(f'Pipeline: {workdir}')
    lines.append('')
    lines.append('=' * 72)
    lines.append('1. HOUSEKEEPING GTIs (STDGTI01 - per-CCD)')
    lines.append('=' * 72)
    hdr = f"{'Inst':<10} {'PPS #int':>10} {'PPS ks':>10} {'Pipe #int':>10} {'Pipe ks':>10} {'Match':>7}"
    lines.append(hdr)
    lines.append('-' * len(hdr))
    pps_evts = {'PN S': 'P0963720201PNS003PIEVLI0000.FTZ', 'MOS1 S': 'P0963720201M1S001MIEVLI0000.FTZ', 'MOS2 S': 'P0963720201M2S002MIEVLI0000.FTZ'}
    pipe_evts = {'PN S': 'PN/4759_0963720201_EPN_S003_ImagingEvts.clean.fits', 'MOS1 S': 'M1/4759_0963720201_EMOS1_S001_ImagingEvts.clean.fits', 'MOS2 S': 'M2/4759_0963720201_EMOS2_S002_ImagingEvts.clean.fits'}
    for label in ['PN S', 'MOS1 S', 'MOS2 S']:
        with fits.open(str(mt_dir / pps_evts[label])) as hdul:
            g = hdul['STDGTI01'].data
            pn, pk = (len(g), np.sum(g['STOP'] - g['START']) / 1000)
        with fits.open(str(workdir / 'clean' / pipe_evts[label])) as hdul:
            g = hdul['STDGTI01'].data
            pipn, pipk = (len(g), np.sum(g['STOP'] - g['START']) / 1000)
        ok = 'YES' if pn == pipn and abs(pk - pipk) < 0.1 else 'NO'
        lines.append(f'{label:<10} {pn:>10d} {pk:>10.1f} {pipn:>10d} {pipk:>10.1f} {ok:>7}')
    lines.append('')
    lines.append('=' * 72)
    lines.append('2. FLARE FILTERING')
    lines.append('   PPS:  SRC_GTIS from FBKTSR (bkg rate screening)')
    lines.append('   Pipe: espfilt (Gaussian fit to rate histogram)')
    lines.append('   espfilt params: PN rangescale=25, MOS=20, allowsigma=3.0')
    lines.append('=' * 72)
    hdr = f"{'Inst':<10} {'PPS ks':>10} {'Pipe ks':>10} {'Ratio':>8}"
    lines.append(hdr)
    lines.append('-' * len(hdr))
    pps_fbk = {'PN': 'P0963720201PNS003FBKTSR0000.FTZ', 'MOS1': 'P0963720201M1S001FBKTSR0000.FTZ', 'MOS2': 'P0963720201M2S002FBKTSR0000.FTZ'}
    pipe_espfilt = {'PN': 'PN/4759_0963720201_EPN_S003_ImagingEvts.clean.espfilt/pnS003-gti.fits', 'MOS1': 'M1/4759_0963720201_EMOS1_S001_ImagingEvts.clean.espfilt/mos1S001-gti.fits', 'MOS2': 'M2/4759_0963720201_EMOS2_S002_ImagingEvts.clean.espfilt/mos2S002-gti.fits'}
    for inst in ['PN', 'MOS1', 'MOS2']:
        with fits.open(str(mt_dir / pps_fbk[inst])) as hdul:
            g = hdul['SRC_GTIS'].data
            pk = np.sum(g['STOP'] - g['START']) / 1000
        esp = workdir / 'clean' / pipe_espfilt[inst]
        with fits.open(str(esp)) as hdul:
            g = hdul['STDGTI'].data
            pipk = np.sum(g['STOP'] - g['START']) / 1000
        ratio = pk / pipk if pipk > 0 else float('inf')
        lines.append(f'{inst:<10} {pk:>10.1f} {pipk:>10.1f} {ratio:>7.1f}x')
    lines.append('')
    lines.append('*** espfilt is 2-3x more aggressive than PPS flare screening. ***')
    lines.append('    Consider raising rangescale further or using the PPS FBKTSR GTIs.')
    lines.append('')
    lines.append('=' * 72)
    lines.append('3. EVENT COUNTS (scheduled)')
    lines.append('   PPS:  unfiltered calibrated events')
    lines.append('   Pipe: PI 200-12000 + STDGTI + espfilt GTI')
    lines.append('=' * 72)
    hdr = f"{'Inst':<10} {'PPS events':>15} {'Pipe events':>15} {'Cut %':>8}"
    lines.append(hdr)
    lines.append('-' * len(hdr))
    for label in ['PN S', 'MOS1 S', 'MOS2 S']:
        pn = event_count_fast(mt_dir / pps_evts[label])
        pipn = event_count_fast(workdir / 'clean' / pipe_evts[label])
        pct = (1 - pipn / pn) * 100 if pn > 0 else 0
        lines.append(f'{label:<10} {pn:>15,} {pipn:>15,} {pct:>7.1f}%')
    lines.append('')
    lines.append('=' * 72)
    lines.append('4. SOURCE DETECTION')
    lines.append('=' * 72)
    pps_src_files = sorted(mt_dir.glob('*OBSMLI*.FTZ'))
    for sf in pps_src_files:
        if 'EPX' in sf.name:
            with fits.open(str(sf)) as hdul:
                d = hdul['SRCLIST'].data
                ml = d['EP_DET_ML']
                lines.append(f'PPS:      {len(d)} EPIC sources (ML {ml.min():.1f} - {ml.max():.1f})')
    pipe_src = workdir / 'detect' / 'field_sources_curated.fits'
    if pipe_src.exists():
        with fits.open(str(pipe_src)) as hdul:
            d = hdul[1].data
            ml = d['DET_ML']
            lines.append(f'Pipeline: {len(d)} curated sources (ML {ml.min():.1f} - {ml.max():.1f}), band PI 1000-12000')
    lines.append('')
    lines.append('Pipeline finds more sources: wider mosaic from pseudoexposure')
    lines.append('stacking across the comet track covers larger sky area.')
    lines.append('')
    lines.append('=' * 72)
    lines.append('5. ATTITUDE')
    lines.append('=' * 72)
    att_files = sorted(mt_dir.glob('*ATTTSR*.FTZ'))
    if att_files:
        with fits.open(str(att_files[0])) as hdul:
            d = hdul['ATTHK'].data
            ra, dec = (d['AHFRA'], d['AHFDEC'])
            dt = (d['TIME'][-1] - d['TIME'][0]) / 1000
            lines.append(f'PPS attitude: {dt:.1f} ks, RA span {(ra.max() - ra.min()) * 3600:.1f}", DEC span {(dec.max() - dec.min()) * 3600:.1f}"')
    track_dir = workdir / 'track'
    if track_dir.exists():
        track_files = sorted(track_dir.iterdir())
        lines.append(f'Pipeline track: {len(track_files)} files in track/')
        for tf in track_files[:5]:
            lines.append(f'  {tf.name}')
    lines.append('')
    lines.append('=' * 72)
    lines.append('6. FLARE RATE STATISTICS (scheduled)')
    lines.append('=' * 72)
    for inst in ['PN', 'MOS1', 'MOS2']:
        with fits.open(str(mt_dir / pps_fbk[inst])) as hdul:
            r = hdul['RATE'].data['RATE']
            good = np.isfinite(r)
            r = r[good]
            lines.append(f'PPS {inst}: median={np.median(r):.2f}, mean={np.mean(r):.2f}, max={r.max():.1f} ct/s, P90={np.percentile(r, 90):.1f}, P95={np.percentile(r, 95):.1f}, P99={np.percentile(r, 99):.1f}')
    lines.append('')
    lines.append('=' * 72)
    lines.append('7. PIPELINE STAGE STATUS')
    lines.append('=' * 72)
    stages = ['init', 'repro', 'clean', 'track', 'detect', 'comet', 'contam', 'image', 'lcurve', 'spectrum', 'merge']
    for s in stages:
        sdir = workdir / s
        if sdir.exists():
            n = len(list(sdir.iterdir()))
            status = 'DONE' if n > 2 else 'partial'
        else:
            n = 0
            status = 'MISSING'
        lines.append(f'  {s:<12} {status:<10} ({n} files)')
    lines.append('')
    text = '\n'.join(lines)
    outpath = outdir / 'pps_vs_pipeline_comparison.txt'
    outpath.write_text(text)
    print(text)
    print(f'\nSaved to: {outpath}')
    return 0
if __name__ == '__main__':
    sys.exit(main())
