import configparser
import struct
import zlib
from pathlib import Path

import numpy as np


def parse_ini(ini_path):
    """Parse Slidedat.ini and return metadata dictionary."""
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
    camera_rotation = float(g.get('CAMERA_ROTATION', 0))

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
        'camera_rotation': camera_rotation,
        'hierarchy': h,
    }


def read_stitching_layer(slide_dir, meta):
    """Read stride values from StitchingLayer.

    Returns (stride_x, stride_y) or None.
    """
    h = meta['hierarchy']
    nonhier_count = int(h.get('nonhier_count', 0))

    ei = 0
    found = False
    for i in range(nonhier_count):
        name = h.get(f'nonhier_{i}_name', '')
        count = int(h.get(f'nonhier_{i}_count', 0))
        if name == 'StitchingLayer':
            found = True
            break
        ei += count

    if not found:
        return None

    index_path = slide_dir / 'Index.dat'
    with open(index_path, 'rb') as f:
        f.read(5)
        f.read(len(meta['slide_id']))
        hier_base = struct.unpack('<i', f.read(4))[0]
        nonhier_base = struct.unpack('<i', f.read(4))[0]
        f.seek(nonhier_base + ei * 4)
        ptr = struct.unpack('<i', f.read(4))[0]
        if ptr <= 0:
            return None
        f.seek(ptr)
        first_page = struct.unpack('<i', f.read(8)[4:8])[0]
        f.seek(first_page)
        n_items = struct.unpack('<i', f.read(4))[0]
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

    if len(data) < 28:
        return None

    stride_x = struct.unpack_from('<H', data, 10)[0]
    stride_y = struct.unpack_from('<H', data, 26)[0]
    return stride_x, stride_y


def read_tile_positions(slide_dir, meta):
    """Read tile positions from MRXS.

    Tries VIMSLIDE_POSITION_BUFFER (uncompressed) first, then
    StitchingIntensityLayer (zlib-compressed). Returns (N, 2) int32 array
    indexed by camera position grid index, or None on failure.
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
