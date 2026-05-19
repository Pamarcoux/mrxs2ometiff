#!/usr/bin/env python3
"""Pixel-level comparison of two OME-TIFFs (reference vs test)."""

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import tifffile

OME_NS = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}


def _read_ome_meta(tif):
    root = ET.fromstring(tif.ome_metadata)
    pixels = root.find('.//ome:Pixels', OME_NS)
    if pixels is None:
        return {}
    return dict(pixels.attrib)


def compare(ref_path, test_path, outdir):
    ref_path = Path(ref_path)
    test_path = Path(test_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    stem = test_path.stem
    print(f"\n  Reference: {ref_path.name}")
    print(f"  Test:      {test_path.name}")

    with tifffile.TiffFile(ref_path) as rf, tifffile.TiffFile(test_path) as tf:
        ref_series = rf.series[0]
        test_series = tf.series[0]
        ref_shape = ref_series.shape
        test_shape = test_series.shape
        ref_dtype = ref_series.dtype
        test_dtype = test_series.dtype

        print(f"\n  Reference shape: {ref_shape}  dtype: {ref_dtype}")
        print(f"  Test      shape: {test_shape}  dtype: {test_dtype}")

        # Metadata comparison
        ref_meta = _read_ome_meta(rf)
        test_meta = _read_ome_meta(tf)
        meta_keys = ['PhysicalSizeX', 'PhysicalSizeY', 'SizeC', 'SizeZ', 'SizeT']
        for k in meta_keys:
            ref_v = ref_meta.get(k, '?')
            test_v = test_meta.get(k, '?')
            status = 'OK' if ref_v == test_v else 'MISMATCH'
            print(f"  {k}: ref={ref_v} test={test_v} [{status}]")

        if ref_dtype != test_dtype:
            print(f"\n  *** DTYPE MISMATCH: {ref_dtype} vs {test_dtype} ***")
            return False

        ref = ref_series.asarray()
        test = test_series.asarray()

    # Handle dimension mismatch — compare cropped intersection
    n_channels = min(ref_shape[0], test_shape[0])
    common_h = min(ref_shape[1], test_shape[1])
    common_w = min(ref_shape[2], test_shape[2])
    same_dims = ref_shape == test_shape

    ref_cropped = ref[:n_channels, :common_h, :common_w]
    test_cropped = test[:n_channels, :common_h, :common_w]

    # Full-image comparison (only if same size)
    if same_dims:
        print(f"\n  --- Full-image comparison ---")
        _report_diff(ref[:n_channels], test[:n_channels], n_channels, outdir, stem, label='full')

    # Cropped comparison (overlap region)
    print(f"\n  --- Cropped common region ({common_h} x {common_w}) ---")
    all_identical = _report_diff(ref_cropped, test_cropped, n_channels, outdir, stem, label='crop')

    # Write summary
    _write_summary(outdir, stem, ref_path, test_path, ref_shape, test_shape,
                   ref_dtype, meta_keys, ref_meta, test_meta,
                   same_dims, common_h, common_w, all_identical)

    return all_identical


def _report_diff(ref_arr, test_arr, n_channels, outdir, stem, label):
    total_pixels = ref_arr.shape[1] * ref_arr.shape[2]

    print(f"  {'Ch':<5} {'Ref min':>8} {'Ref max':>8} {'Ref mean':>9} "
          f"{'Test min':>8} {'Test max':>8} {'Test mean':>9} "
          f"{'Diff max':>8} {'Diff mean':>8} {'Identical%':>10}")
    print(f"  {'-'*5} {'-'*8} {'-'*8} {'-'*9} "
          f"{'-'*8} {'-'*8} {'-'*9} "
          f"{'-'*8} {'-'*8} {'-'*10}")

    all_ok = True
    diff_stats = {}
    for c in range(n_channels):
        ref_ch = ref_arr[c].astype(np.float64)
        test_ch = test_arr[c].astype(np.float64)
        diff = np.abs(ref_ch - test_ch)

        stats = {
            'ref_min': int(ref_arr[c].min()), 'ref_max': int(ref_arr[c].max()),
            'ref_mean': float(ref_arr[c].mean()),
            'test_min': int(test_arr[c].min()), 'test_max': int(test_arr[c].max()),
            'test_mean': float(test_arr[c].mean()),
            'diff_min': float(diff.min()), 'diff_max': float(diff.max()),
            'diff_mean': float(diff.mean()), 'diff_std': float(diff.std()),
            'pct_identical': round(100.0 * np.count_nonzero(diff == 0) / total_pixels, 2),
            'pct_diff_gt_1': round(100.0 * np.count_nonzero(diff > 1) / total_pixels, 4),
        }
        if ref_arr[c].dtype == np.uint16:
            stats['pct_diff_gt_10'] = round(100.0 * np.count_nonzero(diff > 10) / total_pixels, 4)

        diff_stats[f'ch{c}'] = stats

        print(f"  {f'ch{c}':<5} {stats['ref_min']:>8} {stats['ref_max']:>8} {stats['ref_mean']:>9.1f} "
              f"{stats['test_min']:>8} {stats['test_max']:>8} {stats['test_mean']:>9.1f} "
              f"{stats['diff_max']:>8.0f} {stats['diff_mean']:>8.2f} "
              f"{stats['pct_identical']:>9.2f}%")
        if stats['pct_identical'] < 100:
            all_ok = False

    verdict = '>>> PIXEL-IDENTICAL (common region)' if all_ok else '>>> HAS DIFFERENCES'
    print(f"  {verdict}")

    # Write diff image if there are differences (small crops only)
    if not all_ok and ref_arr.shape[1] <= 4096 and ref_arr.shape[2] <= 4096:
        diff_stack = []
        for c in range(n_channels):
            ref_ch = ref_arr[c].astype(np.float64)
            test_ch = test_arr[c].astype(np.float64)
            d = np.abs(ref_ch - test_ch)
            d_norm = np.clip(d / (d.max() + 1e-8) * 65535, 0, 65535).astype(np.uint16)
            diff_stack.append(d_norm)
        diff_img = np.stack(diff_stack, axis=0)
        diff_path = outdir / f'diff_{stem}_{label}.ome.tif'
        tifffile.imwrite(diff_path, diff_img, bigtiff=True,
                         photometric='minisblack', metadata={'axes': 'CYX'})
        print(f"  Diff image: {diff_path}")

    return all_ok


def _write_summary(outdir, stem, ref_path, test_path, ref_shape, test_shape,
                   dtype, meta_keys, ref_meta, test_meta, same_dims,
                   common_h, common_w, all_identical):
    summary_path = outdir / f'summary_{stem}.txt'
    with open(summary_path, 'w') as f:
        f.write(f"Reference: {ref_path.name}\n")
        f.write(f"Test:      {test_path.name}\n")
        f.write(f"Shape:     ref={ref_shape} test={test_shape}\n")
        f.write(f"Dtype:     {dtype}\n")
        f.write(f"Same dims: {same_dims}\n")
        f.write(f"Common region: {common_h} x {common_w}\n\n")
        for k in meta_keys:
            f.write(f"{k}: ref={ref_meta.get(k, '?')} test={test_meta.get(k, '?')}\n")
        f.write(f"\nPixel-identical (common region): {all_identical}\n")
    print(f"  Summary:   {summary_path}")


def main():
    parser = argparse.ArgumentParser(description='Compare two OME-TIFFs pixel-level')
    parser.add_argument('reference', help='Reference TIFF (bioformats2raw+raw2ometiff)')
    parser.add_argument('test', help='Test TIFF (our converter)')
    parser.add_argument('--outdir', default='.', help='Output directory for results')
    args = parser.parse_args()

    ok = compare(args.reference, args.test, args.outdir)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
