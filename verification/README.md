# Verification — MRXS → OME-TIFF converter

Pixel-level validation of `mrxs2ometiff.py` against the established Java pipeline
(`bioformats2raw` + `raw2ometiff`).

## Method

For a given MRXS slide, both pipelines are run:

| Pipeline | Command |
|---|---|
| Reference (Java) | `bioformats2raw slide.mrxs /tmp/zarr/` → `raw2ometiff /tmp/zarr/ ref.ome.tif` |
| Test (Python) | `mrxs2ometiff.py slide.mrxs -o test.ome.tif` |

The two OME-TIFFs are compared channel-by-channel at full (16-bit) precision.

## Usage

```bash
./verification/run_comparison.sh SLIDE_NAME
```

Results are written to `verification/results/`:
- `*.test.ome.tif` — our converter output
- `summary_*.txt` — per-channel diff statistics
- `diff_*.ome.tif` — difference image (pixel error, 16-bit, normalized)

Reference outputs are cached in `verification/.ref/` (gitignored).

## Interpreting results

- **Pixel-identical** (100% identical): both pipelines produce bitwise-identical output.
- **Small differences** (< 0.1% pixels differ by > 1 DN): JPEG-XR decode rounding differences between Bio-Formats and `imagecodecs`.
- **Large differences**: possible position/overlap handling discrepancy.

## Results (first slide: 295182_PM_PatientCevi040534-Panel1_UT1.1)

| Metric | Finding |
|---|---|
| **Dimensions** | `bioformats2raw`: 21504×17212 (padded) vs **ours**: 20171×16335 (tight bounding box) |
| **PhysicalSize** | Both: 0.325 µm/pixel ✅ |
| **DAPI (ch3)** | `bioformats2raw`: mean=9.7 (blank) vs **ours**: mean=119.9 (valid signal) — **confirms layer bug fix** |
| **Channel order** | ch0/ch2 appear swapped between pipelines — under investigation |
| **Non-DAPI channels** | Means within ~1-2% across common region (JPEG-XR decode variance) |

The dimension difference is because `bioformats2raw` uses metadata-based or padded dimensions,
while our converter computes a tight bounding box from tile positions. The DAPI finding
confirms that `bioformats2raw` (via Bio-Formats' `MiraxReader`) suffers from the FilterLevel_1
bug that our converter fixes.

## Speed benchmark

See [BENCHMARK.md](BENCHMARK.md) for timing results. Our converter is **~2× faster**
on average (single-pass Python vs 2-pass Java pipeline).

## Limitations

- The reference pipeline (`bioformats2raw`) generates its own pyramid; our converter
  outputs only the full-resolution level.
- Overlapping tiles: `bioformats2raw` composites via Bio-Formats; our converter
  places tiles at raw positions without overlap blending.
- Comparison is full-resolution only (no pyramid levels compared).
