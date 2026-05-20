# mrxs2ometiff — MRXS to OME-TIFF converter

Convert 3DHistech MRXS (MIRAX) whole-slide images to 16-bit OME-TIFF.

Python-native alternative to the bioformats2raw + raw2ometiff Java pipeline.
No Java or intermediate Zarr files required.

## Features

- Multi-resolution OME-TIFF with 8-level pyramid (like raw2ometiff)
- 4-channel 16-bit output with zlib compression
- Reads all filter levels (no loss of channels stored on FilterLevel_1)
- Full OME metadata (PhysicalSize, channel names, excitation/emission wavelengths)
- Streaming via memmap — memory usage ~3.7 GB regardless of slide size
- Multi-threaded JPEG-XR tile decode
- Single-step conversion — no Java or intermediate Zarr files
- 2–4× faster than the bioformats2raw + raw2ometiff pipeline (tested on 20 slides)

## Limitations

- Tile overlap is not composited — raw tile positions are used
- MRXS format only (not SVS, NDPI, etc.)

## Installation

```bash
pip install mrxs2ometiff
```

Requires Python ≥ 3.9 with numpy, tifffile, and imagecodecs.

## Usage

```bash
# Single slide — specify output path
mrxs2ometiff slide.mrxs -o slide.ome.tif

# Batch mode — convert all .mrxs in a directory
# Output TIFFs are written to ./TIFF/ with the same base names
mrxs2ometiff MRXS_DIR/

# Verify an existing OME-TIFF (dimensions, pixel ranges, OME metadata)
mrxs2ometiff output.ome.tif --verify
mrxs2ometiff TIFF/ --verify     # verify a directory of TIFFs

# Skip pyramid sub-IFDs (smaller output, no multi-resolution)
mrxs2ometiff slide.mrxs -o slide.ome.tif --no-pyramid
```

### Default output path

When `-o` is omitted with a single file, the output is written to `../TIFF/`
relative to the slide directory (one level above the `MRXS/` directory).

### What the converter does

1. Reads `Slidedat.ini` and `Index.dat` from the MRXS directory
2. Parses tile records and camera position maps
3. Decodes JPEG-XR tiles and places them at their absolute pixel coordinates
4. Writes a 16-bit multi-channel BigTIFF with OME-XML metadata

## License

LGPL-2.1-only
