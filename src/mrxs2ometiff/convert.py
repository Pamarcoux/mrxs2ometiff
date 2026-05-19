from pathlib import Path

import numpy as np
from imagecodecs import jpegxr_decode
import tifffile

from .reader import parse_ini, read_records, read_tile_positions


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
    """Convert a single MRXS file to OME-TIFF."""
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

    fl_map = {}
    fl_order = []
    for ch in meta['channels']:
        fl = ch['filter_level']
        if fl not in fl_map:
            fl_map[fl] = []
            fl_order.append(fl)
        fl_map[fl].append(ch)

    tile_positions = read_tile_positions(slide_dir, meta)
    if tile_positions is not None:
        active = np.any(tile_positions != 0, axis=1).sum()
        print(f'  Tile positions: {active} active of {len(tile_positions)} total')
    else:
        print('  No position data, using grid formula')

    stack, ch_list = assemble_level(slide_dir, meta, records, 0, fl_map, fl_order, tile_positions)
    if stack is None:
        print('  ERROR: no data at full resolution')
        return False

    print(f'  Full res: {stack.shape}  (C={stack.shape[0]}, Y={stack.shape[1]}, X={stack.shape[2]})')

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
        stack,
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
