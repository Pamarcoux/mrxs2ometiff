#!/usr/bin/env python3
"""Benchmark our converter vs bioformats2raw+raw2ometiff."""

import subprocess
import sys
import time
import os
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIOFORMATS = REPO / 'bioformats2raw-0.12.0' / 'bin' / 'bioformats2raw'
RAW2OMETIFF = REPO / 'raw2ometiff-0.10.0' / 'bin' / 'raw2ometiff'
OUR_SCRIPT = REPO / 'mrxs2ometiff.py'

RESULTS = REPO / 'verification' / 'results'
RESULTS.mkdir(parents=True, exist_ok=True)
CSV_PATH = REPO / 'verification' / 'benchmark_results.csv'


def human(n):
    for unit in ['B', 'KiB', 'MiB', 'GiB']:
        if abs(n) < 1024:
            return f'{n:.1f} {unit}'
        n /= 1024
    return f'{n:.1f} TiB'


def run(cmd, label):
    """Run a command, return (real_sec, peak_mb)."""
    cmd_str = ' '.join(str(c) for c in cmd)
    print(f'\n  [{label}] Running: {cmd_str}')
    start = time.perf_counter()
    proc = subprocess.run(
        ['/usr/bin/time', '-v'] + [str(c) for c in cmd],
        capture_output=True, text=True, timeout=1800
    )
    elapsed = time.perf_counter() - start
    rc = proc.returncode

    # Parse /usr/bin/time -v output (on stderr)
    peak_mb = 0
    for line in (proc.stderr or '').split('\n'):
        if 'Maximum resident set size' in line:
            peak_kb = int(line.split(':')[1].strip())
            peak_mb = peak_kb / 1024
        if 'Elapsed (wall clock) time' in line:
            # Format: h:mm:ss or m:ss or mm:ss.ss
            part = line.split(':')[1].strip()
            pass  # We use perf_counter for wall time

    if rc != 0:
        print(f'  WARN: exit code {rc}')
        if proc.stderr:
            for line in proc.stderr.split('\n')[-5:]:
                print(f'  stderr: {line}')

    out_size = 0
    return elapsed, peak_mb


def main():
    slides = sys.argv[1:] if len(sys.argv) > 1 else [
        '295209_PM_PatientCevi040534-Panel2_Glofi2-1',   # small
        '295182_PM_PatientCevi040534-Panel1_UT1.1',      # medium
        '295199_PM_PatientCevi040534-Panel2_UT1-1',      # large
    ]

    existing = []
    if CSV_PATH.exists():
        with open(CSV_PATH) as f:
            existing = [l.split(',')[0] for l in f.readlines()[1:]]

    for slide in slides:
        mrxs = REPO / 'MRXS' / f'{slide}.mrxs'
        if not mrxs.exists():
            print(f'SKIP: {mrxs} not found')
            continue

        print(f'\n{"="*70}')
        print(f'  Slide: {slide}')
        print(f'{"="*70}')

        # 1. Our converter
        our_tiff = RESULTS / f'{slide}.test.ome.tif'
        t1, m1 = run(
            [sys.executable, str(OUR_SCRIPT), str(mrxs), '-o', str(our_tiff)],
            'mrxs2ometiff.py'
        )
        our_size = our_tiff.stat().st_size if our_tiff.exists() else 0

        # 2. bioformats2raw
        zarr_dir = Path(f'/tmp/zarr_{slide}')
        if zarr_dir.exists():
            import shutil
            shutil.rmtree(zarr_dir)
        t2a, m2a = run([str(BIOFORMATS), str(mrxs), str(zarr_dir), '--no-minmax'],
                       'bioformats2raw')

        # 3. raw2ometiff
        ref_tiff = RESULTS / f'{slide}.ref.ome.tif'
        t2b, m2b = run([str(RAW2OMETIFF), str(zarr_dir), str(ref_tiff)],
                       'raw2ometiff')
        ref_size = ref_tiff.stat().st_size if ref_tiff.exists() else 0
        t2 = t2a + t2b
        m2 = max(m2a, m2b)

        # Cleanup zarr
        if zarr_dir.exists():
            import shutil
            shutil.rmtree(zarr_dir)

        # Write CSV row
        with open(CSV_PATH, 'a') as f:
            f.write(f'{slide},mrxs2ometiff.py,{t1:.1f},{m1:.0f},{human(our_size)}\n')
            f.write(f'{slide},bioformats2raw+raw2ometiff,{t2:.1f},{m2:.0f},{human(ref_size)}\n')

    print(f'\n{"="*70}')
    print(f'  Results appended to {CSV_PATH}')
    if CSV_PATH.exists():
        with open(CSV_PATH) as f:
            print(f.read())


if __name__ == '__main__':
    main()
