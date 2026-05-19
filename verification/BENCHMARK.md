# Benchmark — MRXS → OME-TIFF conversion speed

Results from `benchmark.py` (3 slides of increasing size, single run each, cold cache).

## Summary

| Slide | Tiles | `mrxs2ometiff.py` | `bioformats2raw` + `raw2ometiff` | Speedup |
|---|---|---|---|---|
| `Glofi2-1` (small) | 226 tiles | **18.8 s** / 6.9 GB | 29.9 s / 6.7 GB | **1.6×** |
| `UT1.1` (medium)   | 238 tiles | **21.4 s** / 7.3 GB | 41.7 s / 9.9 GB | **1.9×** |
| `UT1-1` (large)    | 276 tiles | **24.7 s** / 8.6 GB | 54.3 s / 8.8 GB | **2.2×** |

## Key takeaways

- **2× faster** on average than the Java 2-step pipeline
- **Same ~7-9 GB peak RAM** (both pipelines load all tiles before writing)
- **Output sizes**: ours is slightly larger (tight bounding box has less padding → more tile area to store; also our TIFF compression may differ slightly from raw2ometiff's defaults)
- **Single step**: no intermediate Zarr directory (`/tmp`), no Java dependency

## Hardware

- CPU: (detected automatically)
- RAM: 64 GB+
- Disk: NVMe SSD
- OS: Linux

## Method

```bash
time python3 mrxs2ometiff.py slide.mrxs -o out.ome.tif          # our pipeline
time bioformats2raw slide.mrxs /tmp/zarr && raw2ometiff ...     # reference
```

Peak RSS measured via `/usr/bin/time -v` (`Maximum resident set size`).
