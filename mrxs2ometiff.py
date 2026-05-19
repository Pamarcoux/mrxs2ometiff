#!/usr/bin/env python3
"""Convert MRXS multi-channel fluorescence slides to pyramidal 16-bit OME-TIFF.

Replaces the bioformats2raw+raw2ometiff pipeline which loses the DAPI channel
due to a MiraxReader bug (FilterLevel_1 tile offsets never read).

Dependencies: numpy, tifffile, imagecodecs
"""

import struct
import configparser
import sys
import zlib
from pathlib import Path
import numpy as np
from imagecodecs import jpegxr_decode
import tifffile


def parse_ini(ini_path):
    """Parse Slidedat.ini and return metadata."""
    config = configparser.ConfigParser()
    with open(ini_path, encoding='utf-8-sig') as f:
        config.read_file(f)

    g = config['GENERAL']
    h = config['HIERARCHICAL']

    slide_id = g['SLIDE_ID']
    images_x = int(g['IMAGENUMBER_X'])
    images_y = int(g['IMAGENUMBER_Y'])
    zoom_levels = int(h['HIER_0_COUNT'])

    tile_w = int(config['LAYER_0_LEVEL_0_SECTION']['DIGITIZER_WIDTH'])
    tile_h = int(config['LAYER_0_LEVEL_0_SECTION']['DIGITIZER_HEIGHT'])
    image_divisions = int(g.get('CameraImageDivisionsPerSide', 1))

    zoom_info = {}
    for z in range(zoom_levels):
        sec = config[f'LAYER_0_LEVEL_{z}_SECTION']
        zoom_info[z] = {
            'px_x': float(sec['MICROMETER_PER_PIXEL_X']),
            'px_y': float(sec['MICROMETER_PER_PIXEL_Y']),
            'overlap_x': float(sec.get('overlap_x', 0)),
            'overlap_y': float(sec.get('overlap_y', 0)),
            'concat_exponent': int(sec.get('IMAGE_CONCAT_FACTOR', 0)),
        }

    filter_level_count = int(h['HIER_1_COUNT'])
    channels = []
    for i in range(filter_level_count):
        sec = config[f'LAYER_1_LEVEL_{i}_SECTION']
        channels.append({
            'name': sec['FILTER_NAME'],
            'storing_ch': int(sec['STORING_CHANNEL_NUMBER']),
            'filter_level': sec['DATA_IN_THIS_FILTER_LEVEL'],
            'ex_center': float(sec.get('EXCITATION_WAVELENGTH', 0)),
            'em_center': float(sec.get('EMISSION_WAVELENGTH', 0)),
        })

    df_section = config['DATAFILE']
    data_files = [df_section[f'FILE_{i}'] for i in range(int(df_section['FILE_COUNT']))]

    return {
        'slide_id': slide_id,
        'images_x': images_x,
        'images_y': images_y,
        'zoom_levels': zoom_levels,
        'tile_w': tile_w,
        'tile_h': tile_h,
        'image_divisions': image_divisions,
        'zoom_info': zoom_info,
        'channels': channels,
        'data_files': data_files,
        'filter_level_count': filter_level_count,
        'datafile_count': len(data_files),
        'level_0_image_concat': 1 << zoom_info[0]['concat_exponent'],
        'hierarchy': h,
    }


def read_tile_positions(slide_dir, meta):
    """Read raw tile positions from MRXS.

    Tries VIMSLIDE_POSITION_BUFFER (uncompressed) first, then
    StitchingIntensityLayer (zlib-compressed). Returns (N, 2) int32 array
    of raw (x, y) pixel positions indexed by camera position grid index,
    multiplied by level_0_image_concat. Returns None on failure.
    """
    images_x = meta['images_x']
    image_divisions = meta['image_divisions']
    h = meta['hierarchy']

    nonhier_count = int(h.get('nonhier_count', 0))

    candidates = [('VIMSLIDE_POSITION_BUFFER', False),
                  ('StitchingIntensityLayer', True)]

    entry_idx = None
    is_compressed = None
    for target_name, compressed in candidates:
        ei = 0
        for i in range(nonhier_count):
            name = h.get(f'nonhier_{i}_name', '')
            count = int(h.get(f'nonhier_{i}_count', 0))
            if name == target_name:
                entry_idx = ei
                is_compressed = compressed
                break
            ei += count
        if entry_idx is not None:
            break

    if entry_idx is None:
        return None

    # Read Index.dat to locate the data
    index_path = slide_dir / 'Index.dat'
    with open(index_path, 'rb') as f:
        f.read(5)
        f.read(len(meta['slide_id']))
        hier_base = struct.unpack('<i', f.read(4))[0]
        nonhier_base = struct.unpack('<i', f.read(4))[0]

        f.seek(nonhier_base + entry_idx * 4)
        ptr = struct.unpack('<i', f.read(4))[0]
        if ptr <= 0:
            return None

        f.seek(ptr)
        first_page = struct.unpack('<i', f.read(8)[4:8])[0]

        f.seek(first_page)
        n_items = struct.unpack('<i', f.read(4))[0]
        if n_items == 0:
            return None
        f.read(12)

        ref_offset = struct.unpack('<i', f.read(4))[0]
        ref_length = struct.unpack('<i', f.read(4))[0]
        ref_file = struct.unpack('<i', f.read(4))[0]

    df_path = slide_dir / f'Data{ref_file:04d}.dat'
    if not df_path.exists():
        return None

    with open(df_path, 'rb') as df:
        df.seek(ref_offset)
        data = df.read(ref_length)

    if is_compressed:
        try:
            data = zlib.decompress(data)
        except Exception:
            return None

    n_expected = (images_x // image_divisions) * (meta['images_y'] // image_divisions)
    step = 9
    level_0_image_concat = meta['level_0_image_concat']
    positions = np.zeros((n_expected, 2), dtype=np.int32)
    for i in range(min(len(data) // step, n_expected)):
        rec = data[i * step:(i + 1) * step]
        flag = rec[0]
        if flag & 0xfe:
            continue
        x = struct.unpack('<i', rec[1:5])[0]
        y = struct.unpack('<i', rec[5:9])[0]
        positions[i] = (x * level_0_image_concat, y * level_0_image_concat)

    return positions


def read_records(index_path, slide_id, zoom_levels, meta=None):
    """Read all index records. Returns list of (idx, offset, size, fileno) per record."""
    with open(index_path, 'rb') as f:
        ver = f.read(5).decode('ascii')
        assert ver == '01.02', f'Unknown version: {ver}'
        sid = f.read(len(slide_id)).decode('ascii')
        assert sid == slide_id, f'Slide ID mismatch: {sid}'

        hier_base = struct.unpack('<i', f.read(4))[0]
        nonhier_base = struct.unpack('<i', f.read(4))[0]
        num_records = (nonhier_base - hier_base) // 4

        datafile_count = len(meta['data_files']) if meta else 0
        img_x = meta['images_x'] if meta else 0
        img_y = meta['images_y'] if meta else 0

        records = []
        for ri in range(num_records):
            f.seek(hier_base + ri * 4)
            ptr = struct.unpack('<i', f.read(4))[0]
            if ptr <= 0:
                records.append([])
                continue

            try:
                f.seek(ptr)
                data = f.read(8)
                if len(data) < 8:
                    records.append([])
                    continue
                init_zero = struct.unpack('<i', data[:4])[0]
                first_page = struct.unpack('<i', data[4:8])[0]
            except (struct.error, OSError):
                records.append([])
                continue

            if init_zero != 0 or first_page <= 0:
                records.append([])
                continue

            entries = []
            page = first_page
            while page > 0:
                try:
                    f.seek(page)
                    data = f.read(8)
                    if len(data) < 8:
                        break
                    entry_count = struct.unpack('<i', data[:4])[0]
                    next_page = struct.unpack('<i', data[4:8])[0]

                    entry_data = f.read(entry_count * 16)
                    if len(entry_data) < entry_count * 16:
                        break

                    for ei in range(entry_count):
                        base = ei * 16
                        image_idx = struct.unpack('<i', entry_data[base:base+4])[0]
                        offset = struct.unpack('<i', entry_data[base+4:base+8])[0]
                        size = struct.unpack('<i', entry_data[base+8:base+12])[0]
                        fileno = struct.unpack('<i', entry_data[base+12:base+16])[0]
                        if image_idx < 0 or offset < 0 or size < 0 or fileno < 0:
                            continue
                        if datafile_count and fileno >= datafile_count:
                            continue
                        gy = image_idx // img_x if img_x else 0
                        if img_y and gy >= img_y:
                            continue
                        entries.append((image_idx, offset, size, fileno))
                except (struct.error, OSError):
                    break
                page = next_page
            records.append(entries)

    return records


def assemble_level(slide_dir, meta, records, zoom_level, fl_map, fl_order, tile_positions=None):
    """Assemble all channels at one zoom level. Returns (C, Y, X) uint16 array and channel list."""
    tile_w = meta['tile_w']
    tile_h = meta['tile_h']
    images_x = meta['images_x']
    image_divisions = meta['image_divisions']
    zoom_levels = meta['zoom_levels']
    data_files = meta['data_files']

    channel_arrays = []
    channel_list = []

    for fl_name in fl_order:
        ch_list = fl_map[fl_name]
        if not ch_list:
            continue

        if fl_name == 'FilterLevel_0':
            rec_idx = zoom_level
        else:
            rec_idx = zoom_levels + zoom_level

        if rec_idx >= len(records):
            continue
        entries = records[rec_idx]
        if not entries:
            continue

        tile_data = []
        for idx, off, sz, fn in entries:
            gx = idx % images_x
            gy = idx // images_x

            if tile_positions is not None:
                cp = (gy // image_divisions) * (images_x // image_divisions) + (gx // image_divisions)
                pos_x = int(tile_positions[cp, 0])
                pos_y = int(tile_positions[cp, 1])
                # Skip inactive positions (OpenSlide: (0,0) that aren't origin)
                if pos_x == 0 and pos_y == 0 and cp != 0:
                    continue
                intra_x = tile_w * (gx % image_divisions)
                intra_y = tile_h * (gy % image_divisions)
                px = pos_x + intra_x
                py = pos_y + intra_y
            else:
                overlap_x = meta['zoom_info'][zoom_level]['overlap_x']
                overlap_y = meta['zoom_info'][zoom_level]['overlap_y']
                px = int(gx * (tile_w - overlap_x))
                py = int(gy * (tile_h - overlap_y))

            tile_data.append((px, py, off, sz, fn))

        if not tile_data:
            continue

        min_x = min(p[0] for p in tile_data)
        max_x = max(p[0] + tile_w for p in tile_data)
        min_y = min(p[1] for p in tile_data)
        max_y = max(p[1] + tile_h for p in tile_data)

        img_w = max_x - min_x
        img_h = max_y - min_y

        n_ch = len(ch_list)
        buffers = [np.zeros((img_h, img_w), dtype=np.uint16) for _ in range(n_ch)]

        for px, py, off, sz, fn in tile_data:
            data_path = slide_dir / data_files[fn]
            with open(data_path, 'rb') as fh:
                fh.seek(off)
                tile = jpegxr_decode(fh.read(sz))

            dx = px - min_x
            dy = py - min_y
            th, tw = tile.shape[:2]
            th = min(th, img_h - dy)
            tw = min(tw, img_w - dx)
            if th <= 0 or tw <= 0:
                continue

            for ci, ch in enumerate(ch_list):
                storing_ch = ch['storing_ch']
                buffers[ci][dy:dy+th, dx:dx+tw] = tile[:th, :tw, storing_ch]

        channel_arrays.extend(buffers)
        channel_list.extend(ch_list)

    if not channel_arrays:
        return None, []

    return np.stack(channel_arrays, axis=0).astype(np.uint16), channel_list


def convert_one(mrxs_path, output_path):
    """Convert a single MRXS file to pyramidal OME-TIFF."""
    mrxs_path = Path(mrxs_path)
    slide_dir = mrxs_path.with_suffix('')
    output_path = Path(output_path)

    print(f'\n=== Converting: {mrxs_path.name} ===')

    meta = parse_ini(slide_dir / 'Slidedat.ini')
    print(f'  Slide ID: {meta["slide_id"]}')
    print(f'  Zoom levels: {meta["zoom_levels"]}')
    print(f'  Channels: {[c["name"] for c in meta["channels"]]}')

    records = read_records(slide_dir / 'Index.dat', meta['slide_id'], meta['zoom_levels'], meta)
    print(f'  Index records: {len(records)}')

    # Group channels by filter level, preserving ini order
    fl_map = {}
    fl_order = []
    for ch in meta['channels']:
        fl = ch['filter_level']
        if fl not in fl_map:
            fl_map[fl] = []
            fl_order.append(fl)
        fl_map[fl].append(ch)

    # Read tile positions from StitchingIntensityLayer
    tile_positions = read_tile_positions(slide_dir, meta)
    if tile_positions is not None:
        active = np.any(tile_positions != 0, axis=1).sum()
        print(f'  Tile positions: {active} active of {len(tile_positions)} total')
    else:
        print('  No position data, using grid formula')

    # Read full-resolution level (zoom 0)
    stack, ch_list = assemble_level(slide_dir, meta, records, 0, fl_map, fl_order, tile_positions)
    if stack is None:
        print('  ERROR: no data at full resolution')
        return False

    full_res = stack
    print(f'  Full res: {full_res.shape}  (C={full_res.shape[0]}, Y={full_res.shape[1]}, X={full_res.shape[2]})')

    # Build channel metadata for OME
    ome_channels = []
    for ch in ch_list:
        entry = {'Name': ch['name']}
        if ch['ex_center']:
            entry['ExcitationWavelength'] = ch['ex_center']
            entry['ExcitationWavelengthUnit'] = 'nm'
        if ch['em_center']:
            entry['EmissionWavelength'] = ch['em_center']
            entry['EmissionWavelengthUnit'] = 'nm'
        ome_channels.append(entry)

    pixel_size = meta['zoom_info'][0]['px_x']

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        output_path,
        full_res,
        bigtiff=True,
        photometric='minisblack',
        tile=(1024, 1024),
        compression='zlib',
        metadata={
            'axes': 'CYX',
            'PhysicalSizeX': pixel_size,
            'PhysicalSizeXUnit': '\u00b5m',
            'PhysicalSizeY': pixel_size,
            'PhysicalSizeYUnit': '\u00b5m',
            'Channel': ome_channels,
        },
    )

    print(f'  Written: {output_path}')
    return True


def verify_tiff(tiff_path):
    """Verify OME-TIFF has correct structure and all channels have data."""
    print(f'\nVerifying: {tiff_path.name}')
    with tifffile.TiffFile(tiff_path) as tif:
        print(f'  BigTIFF: {tif.is_bigtiff}')
        print(f'  Series: {len(tif.series)}')
        for si, s in enumerate(tif.series):
            print(f'    Series {si}: shape={s.shape}, dtype={s.dtype}')

        # Check each channel pixel range
        series = tif.series[0]
        stack = series.asarray()
        for ci in range(stack.shape[0]):
            ch_data = stack[ci]
            print(f'    Channel {ci}: min={ch_data.min()}, max={ch_data.max()}, mean={ch_data.mean():.1f}')
        # Check OME XML
        import xml.etree.ElementTree as ET
        root = ET.fromstring(tif.ome_metadata)
        ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
        pixels = root.find('.//ome:Pixels', ns)
        if pixels is not None:
            print(f'  OME: SizeC={pixels.get("SizeC")}, SizeX={pixels.get("SizeX")}, SizeY={pixels.get("SizeY")}')
            for ch_el in pixels.findall('ome:Channel', ns):
                name = ch_el.get('Name', '<unnamed>')
                ex = ch_el.get('ExcitationWavelength', '')
                em = ch_el.get('EmissionWavelength', '')
                print(f'    Channel: {name}  ex={ex}  em={em}')


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Convert MRXS to pyramidal OME-TIFF')
    parser.add_argument('input', help='MRXS file or directory')
    parser.add_argument('-o', '--output', help='Output OME-TIFF path (default: input name in TIFF/)')
    parser.add_argument('--verify', action='store_true', help='Verify existing TIFF')
    args = parser.parse_args()

    input_path = Path(args.input)

    if args.verify:
        if input_path.is_dir():
            for f in sorted(input_path.glob('*.ome.tif')):
                verify_tiff(f)
        else:
            verify_tiff(input_path)
        return

    if input_path.is_dir():
        # Batch mode: process all .mrxs files in directory
        mrxs_dir = input_path
        output_dir = mrxs_dir.parent / 'TIFF'
        output_dir.mkdir(parents=True, exist_ok=True)
        for f in sorted(mrxs_dir.glob('*.mrxs')):
            out = output_dir / f.with_suffix('.ome.tif').name
            convert_one(f, out)
    else:
        if args.output:
            convert_one(input_path, Path(args.output))
        else:
            output_dir = input_path.parent.parent / 'TIFF'
            output_dir.mkdir(parents=True, exist_ok=True)
            out = output_dir / input_path.with_suffix('.ome.tif').name
            convert_one(input_path, out)


if __name__ == '__main__':
    main()
