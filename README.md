# mrxs2ometiff — MRXS to OME-TIFF converter

Convert 3DHistech MRXS (MIRAX) whole-slide images to 16-bit OME-TIFF with multi-resolution pyramids.

Python-native, single-step alternative to the `bioformats2raw + raw2ometiff` Java pipeline.
Reads **all filter levels** and preserves excitation/emission wavelengths.

## Features

- **Single-step conversion** — no Java, no intermediate Zarr files
- **Multi-resolution OME-TIFF** — 8-level pyramid sub-IFDs (like raw2ometiff)
- **All channels preserved** — reads both `FilterLevel_0` and `FilterLevel_1`
- **Full OME metadata** — channel names, excitation/emission wavelengths, physical pixel size
- **Streaming via memmap** — memory usage ~3.7 GB regardless of slide size
- **Multi-threaded decode** — JPEG-XR tile decoding runs in parallel (8 workers, ~4× speedup)
- **No padding** — output dimensions match the slide's tight bounding box (no waste)
- **16-bit zlib-compressed BigTIFF** — ready for downstream analysis

## Benchmark

| Metric | `mrxs2ometiff` | `bioformats2raw + raw2ometiff` |
|---|---|---|
| **Batch 20 slides** | **~20 s / slide** | ~40 s / slide |
| **Output size** | 0.8–1.8 GB (tight bounds) | ~4.0 GB (padded) |
| **Ex/Em wavelengths** | ✅ Preserved | Lost |

Tested on Linux, 64 GB RAM, NVMe SSD.

## Differences with the Java pipeline

The standard `bioformats2raw + raw2ometiff` pipeline does not handle all the features of the MRXS format. In particular, slides with more than 3 channels can store data on multiple filter levels — the Java pipeline only reads `FilterLevel_0`, so any channels stored on `FilterLevel_1` are lost or garbled. Our tool reads all filter levels and preserves every channel correctly.

Other differences:

- **Excitation/emission wavelengths** — Java outputs empty values while we keep the metadata from `Slidedat.ini`.
- **Channel ordering** — the Java pipeline reorders some channels when reading from multiple filter levels.
- **Output padding** — Java pads images to tile-grid boundaries (21504×17212), making files ~4× larger than needed.

## Installation

```bash
pip install git+https://github.com/paulmarcoux/mrxs2ometiff
```

Requires Python ≥ 3.9 with `numpy`, `tifffile`, and [imagecodecs](https://github.com/cgohlke/imagecodecs) (for JPEG-XR decode and zlib compression).

## Usage

```bash
# Convert a single slide
mrxs2ometiff slide.mrxs -o slide.ome.tif

# Batch convert all .mrxs in a directory (outputs to ../TIFF/)
mrxs2ometiff MRXS_DIR/

# Verify an existing OME-TIFF
mrxs2ometiff output.ome.tif --verify
mrxs2ometiff TIFF/ --verify

# Skip pyramids (faster, smaller, single resolution)
mrxs2ometiff slide.mrxs -o slide.ome.tif --no-pyramid
```

## How it works

### MRXS format

An `.mrxs` file is actually a directory with that extension, containing:

| File | Purpose |
|---|---|
| `Slidedat.ini` | Metadata: slide ID, channels, zoom levels, tile dimensions |
| `Index.dat` | Binary index: page-linked lists mapping image indices to Data file offsets |
| `Data*.dat` | JPEG-XR compressed tile data (typically 4–10 files) |
| `Thumbnail.dat` | Overview thumbnail (optional) |

#### Slidedat.ini structure

```
[GENERAL]
SLIDE_ID = 295182_PM_Panel1_UT1.1
IMAGENUMBER_X = 21, IMAGENUMBER_Y = 17  ← camera grid size
CameraImageDivisionsPerSide = 2         ← each camera position → 2×2 tiles

[HIERARCHICAL]
HIER_0_COUNT = 1      ← only level 0 has actual image data
HIER_1_COUNT = 4      ← 4 filter levels (channels)

[LAYER_0_LEVEL_0_SECTION]
DIGITIZER_WIDTH = 1024, DIGITIZER_HEIGHT = 1024  ← tile dimensions (pixels)
MICROMETER_PER_PIXEL_X = 0.325, _Y = 0.325       ← physical pixel size

[LAYER_1_LEVEL_0_SECTION]  ← channel 0
FILTER_NAME = AX SpGold
DATA_IN_THIS_FILTER_LEVEL = FilterLevel_0
STORING_CHANNEL_NUMBER = 0
EXCITATION_WAVELENGTH = ... / EMISSION_WAVELENGTH = ...

[LAYER_1_LEVEL_1_SECTION]  ← channel 1
FILTER_NAME = AX DA+FI+TR+Cy5-2
DATA_IN_THIS_FILTER_LEVEL = FilterLevel_0
STORING_CHANNEL_NUMBER = 1

[LAYER_1_LEVEL_2_SECTION]  ← channel 2
FILTER_NAME = AX DA+FI+TR+Cy5-4
DATA_IN_THIS_FILTER_LEVEL = FilterLevel_0
STORING_CHANNEL_NUMBER = 2

[LAYER_1_LEVEL_3_SECTION]  ← channel 3
FILTER_NAME = AX DA+FI+TR+Cy5-1
DATA_IN_THIS_FILTER_LEVEL = FilterLevel_1   ← different filter level!
STORING_CHANNEL_NUMBER = 0
```

**Channels can span different filter levels**, stored as separate sections in `Index.dat`.

### Pipeline

```
MRXS directory
  ├─ Slidedat.ini    → parse_ini()        → channel list, tile geometry
  ├─ Index.dat       → read_records()     → tile index (offset, size, file)
  ├─ Index.dat       → read_tile_positions() → camera position grid
  └─ Data*.dat       → JPEG-XR decode     → raw pixel tiles
       │
       ▼
  _build_tile_index()   → tile_map: pixel coordinate → channel list
       │
       ▼
  memmap (CHW format)   → decode_tile() writes each tile at its absolute position
       │
       ▼
  _downsample() × N     → generate pyramid levels via 2×2 block averaging
       │
       ▼
  TiffWriter.write()    → BigTIFF with sub-IFDs + OME-XML metadata
```

### Channel ordering

Channels are grouped by filter level (`FilterLevel_0` first, then `FilterLevel_1`), preserving the `filter_level` ordering from `Slidedat.ini`. Within each group, channels maintain their `storing_ch` index. This means the output channel order may differ from the Java pipeline.

### Pyramids

When `--no-pyramid` is not set, the converter generates additional resolution levels by 2×2 block averaging of the previous level, stopping when both dimensions are ≤ 256 pixels. Each level is written as a sub-IFD in the BigTIFF, compatible with OME-TIFF consumers (Napari, QuPath, OMERO, etc.).

Levels are generated sequentially: each level is downsampled from the level above, written to a temporary memmap file, and linked as a sub-IFD. Temp files are cleaned up after the full TIFF is written.

### Memory streaming

Rather than loading all tiles into a Python list and assembling in-memory (which peaks at ~7 GB for a typical slide), the converter:

1. Creates a temp file for the full-resolution array
2. Opens it as a writable `numpy.memmap` in CHW format (C×H×W, uint16)
3. Decodes tiles in parallel, each writing directly to its memmap region
4. Reads from memmap during pyramid generation

This keeps peak RSS at ~3.7 GB (mostly the JPEG-XR decoder buffers and Python overhead) regardless of slide size. The memmap file is automatically cleaned up on completion.

### Tile decoding

Tiles are JPEG-XR compressed and stored in `.dat` files. The `Index.dat` page-linked lists map each tile to `(data_file, offset, size)`. Camera position data (from `VIMSLIDE_POSITION_BUFFER` or `StitchingIntensityLayer`) gives the absolute pixel coordinate for each grid position.

Tile overlap is recorded in the metadata but not composited — tiles are placed at their recorded positions; overlapping regions will contain the last-written tile's data (no blending). This matches Java's behavior.

## Limitations

- MRXS format only (not SVS, NDPI, CZI, etc.)
- Tile overlap is not composited (matching Java's behavior)
- JPEG-XR dependency (`imagecodecs` via `glymur`) — requires `libjpegxr` on Linux

## License

LGPL-2.1-only
