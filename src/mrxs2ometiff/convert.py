import os
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from imagecodecs import jpegxr_decode
import tifffile

from .reader import parse_ini, read_records, read_tile_positions, read_stitching_layer


def _compute_stride(tile_positions, images_x, images_y, tile_w, tile_h):
    """Compute median stride from non-zero tile positions."""
    poss = tile_positions[np.any(tile_positions != 0, axis=1)]
    if len(poss) < 2:
        return None, None

    gx_est = np.round(poss[:, 0] / (tile_w * 0.5)).astype(int)
    gy_est = np.round(poss[:, 1] / (tile_h * 0.5)).astype(int)

    dxs, dys = [], []
    for gy in np.unique(gy_est):
        mask = gy_est == gy
        xs = np.sort(poss[mask, 0])
        if len(xs) > 1:
            dxs.extend(np.diff(xs).tolist())
    for gx in np.unique(gx_est):
        mask = gx_est == gx
        ys = np.sort(poss[mask, 1])
        if len(ys) > 1:
            dys.extend(np.diff(ys).tolist())

    stride_x = int(round(np.median(dxs))) if dxs else None
    stride_y = int(round(np.median(dys))) if dys else None
    return stride_x, stride_y


def _register_overlap(new_band, existing_band, max_shift=4):
    """Find vertical shift that maximizes NCC between two overlap bands.

    new_band, existing_band: (H, W) arrays — typically a horizontal
    overlap strip (~112 px wide, full tile height).
    Returns best integer shift in [-max_shift, max_shift].
    """
    best_corr = -1.0
    best_shift = 0
    th = new_band.shape[0]
    for shift in range(-max_shift, max_shift + 1):
        if shift < 0:
            n = new_band[:shift, :].ravel()
            e = existing_band[-shift:, :].ravel()
        elif shift > 0:
            n = new_band[shift:, :].ravel()
            e = existing_band[:-shift, :].ravel()
        else:
            n = new_band.ravel()
            e = existing_band.ravel()

        nm = n.mean()
        em = e.mean()
        nc = n.astype(np.float32) - nm
        ec = e.astype(np.float32) - em
        denom = np.linalg.norm(nc) * np.linalg.norm(ec)
        corr = np.dot(nc, ec) / denom if denom > 0 else 0.0
        if corr > best_corr:
            best_corr = corr
            best_shift = shift
    return best_shift


def _blend_horiz(existing, new):
    """Linear blend in horizontal overlap — left (existing) → right (new)."""
    h, w = existing.shape
    alpha = np.linspace(0.0, 1.0, w, dtype=np.float32)
    return (existing * (1 - alpha) + new * alpha).astype(np.uint16)


def _blend_vert(existing, new):
    """Linear blend in vertical overlap — top (existing) → bottom (new)."""
    h, w = existing.shape
    alpha = np.linspace(0.0, 1.0, h, dtype=np.float32)[:, None]
    return (existing * (1 - alpha) + new * alpha).astype(np.uint16)


def _build_tile_index(meta, records, fl_map, fl_order, tile_positions,
                      stitch_stride=None):
    tile_w = meta['tile_w']
    tile_h = meta['tile_h']
    images_x = meta['images_x']
    image_divisions = meta['image_divisions']
    zoom_levels = meta['zoom_levels']

    stride_x, stride_y = None, None
    if stitch_stride is not None:
        stride_x, stride_y = stitch_stride
    elif tile_positions is not None:
        stride_x, stride_y = _compute_stride(
            tile_positions, images_x, meta['images_y'], tile_w, tile_h,
        )

    ch_list = []
    tile_map = defaultdict(list)
    min_x = min_y = float('inf')
    max_x = max_y = float('-inf')

    for fl_name in fl_order:
        ch_group = fl_map[fl_name]
        if not ch_group:
            continue

        if fl_name == 'FilterLevel_0':
            rec_idx = 0
        else:
            rec_idx = zoom_levels

        if rec_idx >= len(records):
            ch_list.extend(ch_group)
            continue
        record_entries = records[rec_idx]
        if not record_entries:
            ch_list.extend(ch_group)
            continue

        base_idx = len(ch_list)
        ch_list.extend(ch_group)

        for idx, off, sz, fn in record_entries:
            gx = idx % images_x
            gy = idx // images_x

            if stride_x is not None and stride_y is not None:
                px = gx * stride_x
                py = gy * stride_y
            elif tile_positions is not None:
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
                overlap_x = meta['zoom_info'][0]['overlap_x']
                overlap_y = meta['zoom_info'][0]['overlap_y']
                px = int(gx * (tile_w - overlap_x))
                py = int(gy * (tile_h - overlap_y))

            min_x = min(min_x, px)
            min_y = min(min_y, py)
            max_x = max(max_x, px + tile_w)
            max_y = max(max_y, py + tile_h)

            key = (off, sz, fn, py, px)
            for ci, ch in enumerate(ch_group):
                tile_map[key].append((base_idx + ci, ch['storing_ch']))

    img_w = int(max_x - min_x)
    img_h = int(max_y - min_y)
    return tile_map, len(ch_list), img_w, img_h, min_x, min_y, ch_list, stride_x, stride_y


def _pyramid_levels(H, W):
    levels = 1
    while max(H, W) > 256:
        levels += 1
        H //= 2
        W //= 2
    return levels


def _downsample(src_mm, C):
    C, H, W = src_mm.shape
    sh, sw = H // 2, W // 2
    pyr = np.empty((C, sh, sw), dtype=np.uint16)
    for y in range(0, sh, 1024):
        bh = min(1024, sh - y)
        src_block = np.array(
            src_mm[:, y*2:y*2+bh*2, :sw*2], dtype=np.uint16
        )
        pyr[:, y:y+bh, :] = (
            src_block.reshape(C, bh, 2, sw, 2).mean(axis=(2, 4)).astype(np.uint16)
        )
    return pyr


def _ome_metadata(ch_list, pixel_size):
    channels = []
    for ch in ch_list:
        entry = {'Name': ch['name']}
        if ch['ex_center']:
            entry['ExcitationWavelength'] = ch['ex_center']
            entry['ExcitationWavelengthUnit'] = 'nm'
        if ch['em_center']:
            entry['EmissionWavelength'] = ch['em_center']
            entry['EmissionWavelengthUnit'] = 'nm'
        channels.append(entry)
    return {
        'PhysicalSizeX': pixel_size,
        'PhysicalSizeXUnit': '\u00b5m',
        'PhysicalSizeY': pixel_size,
        'PhysicalSizeYUnit': '\u00b5m',
        'Channel': channels,
    }


def convert_one(mrxs_path, output_path, pyramid=True):
    mrxs_path = Path(mrxs_path)
    slide_dir = mrxs_path.with_suffix('')
    output_path = Path(output_path)

    print(f'\n=== Converting: {mrxs_path.name} ===')

    meta = parse_ini(slide_dir / 'Slidedat.ini')
    print(f'  Slide ID: {meta["slide_id"]}')
    print(f'  Zoom levels: {meta["zoom_levels"]}')
    print(f'  Channels: {[c["name"] for c in meta["channels"]]}')

    records = read_records(
        slide_dir / 'Index.dat', meta['slide_id'], meta['zoom_levels'], meta,
    )
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

    stitch_stride = read_stitching_layer(slide_dir, meta)
    if stitch_stride is not None:
        sx, sy = stitch_stride
        print(f'  StitchingLayer stride: X={sx}, Y={sy}')
    else:
        print('  No StitchingLayer data')

    tile_map, C, img_w, img_h, min_x, min_y, ch_list, stride_x, stride_y = (
        _build_tile_index(meta, records, fl_map, fl_order, tile_positions,
                          stitch_stride=stitch_stride)
    )

    if stitch_stride is not None:
        print(f'  Using grid stride: X={stride_x}, Y={stride_y}')
    elif stride_x is not None:
        print(f'  Computed stride: X={stride_x}, Y={stride_y}')
    else:
        print('  Using overlap-based formula')

    if C == 0:
        print('  ERROR: no channels')
        return False

    n_levels = _pyramid_levels(img_h, img_w) if pyramid else 1
    gb = C * img_h * img_w * 2 / 1e9
    print(
        f'  Full res: C={C}, Y={img_h}, X={img_w} '
        f'({gb:.1f} GB raw, streamed via temp file)'
    )
    if n_levels > 1:
        print(f'  Pyramids: {n_levels - 1} additional levels')

    tmpdir = output_path.parent
    tmpdir.mkdir(parents=True, exist_ok=True)

    data_file = tempfile.NamedTemporaryFile(
        dir=tmpdir, prefix=f'.{output_path.stem}_', suffix='.raw', delete=False,
    )
    tmp_path = data_file.name
    data_file.close()

    try:
        mm = np.memmap(
            tmp_path, dtype='uint16', mode='write', shape=(C, img_h, img_w),
        )

        data_files = meta['data_files']
        tile_w = meta['tile_w']
        tile_h = meta['tile_h']

        def decode_tile(item):
            key, channels = item
            off, sz, fn, py, px = key
            data_path = slide_dir / data_files[fn]
            with open(data_path, 'rb') as fh:
                fh.seek(off)
                tile = jpegxr_decode(fh.read(sz))
            return key, channels, tile

        # Sort tiles row-major so neighbors are already placed
        sorted_items = sorted(
            tile_map.items(), key=lambda kv: (kv[0][3], kv[0][4])
        )

        overlap_x = tile_w - stride_x if stride_x else 0
        overlap_y = tile_h - stride_y if stride_y else 0

        print(f'  Overlap: X={overlap_x} px, Y={overlap_y} px')
        do_register = (stride_x is not None and stride_y is not None
                       and overlap_x >= 8 and overlap_y >= 8)

        if do_register:
            print('  Registration + blending enabled')
        else:
            print('  Registration disabled (no valid stride or small overlap)')

        with ThreadPoolExecutor(max_workers=8) as pool:
            for key, channels, tile in pool.map(decode_tile, sorted_items):
                _, _, _, py, px = key
                dx = px - min_x
                dy_orig = py - min_y
                th, tw = tile.shape[:2]
                th = min(th, img_h - dy_orig)
                tw = min(tw, img_w - dx)
                if th <= 0 or tw <= 0:
                    continue

                shift_v = 0
                shift_h = 0

                if do_register and overlap_x > 0 and overlap_y > 0:
                    ref_ch = channels[0][1]
                    # --- Horizontal overlap (left neighbor) ---
                    # Left neighbor's right edge overlaps our left edge at
                    # mosaic columns [dx, dx+overlap_x]; its data is already
                    # placed (row-major order).
                    if dx > 0:
                        overlap_w = min(overlap_x, tw, dx)
                        if overlap_w >= 8:
                            new_band = tile[:th, :overlap_w, ref_ch]
                            exist_band = mm[ref_ch, dy_orig:dy_orig+th,
                                            dx:dx+overlap_w]
                            if (new_band.shape == exist_band.shape
                                    and new_band.size > 0):
                                shift_v = _register_overlap(
                                    new_band, exist_band, max_shift=4
                                )

                    # --- Vertical overlap (top neighbor) ---
                    # Top neighbor's bottom edge overlaps our top edge at
                    # mosaic rows [dy_orig, dy_orig+overlap_y].
                    if dy_orig > 0:
                        overlap_h = min(overlap_y, th, dy_orig)
                        if overlap_h >= 8:
                            new_band = tile[:overlap_h, :tw, ref_ch]
                            exist_band = mm[ref_ch,
                                            dy_orig:dy_orig+overlap_h,
                                            dx:dx+tw]
                            if (new_band.shape == exist_band.shape
                                    and new_band.size > 0):
                                shift_h = _register_overlap(
                                    new_band.T, exist_band.T, max_shift=4
                                )

                # Apply shift to placement
                dy = dy_orig + shift_v
                if dy < 0:
                    tile = tile[-dy:, :, :]
                    th -= -dy
                    dy = 0
                if dy + th > img_h:
                    th = img_h - dy
                th = min(th, tile.shape[0])

                # --- Write with blending ---
                for ch_idx, storing_ch in channels:
                    tile_ch = tile[:th, :tw, storing_ch]

                    # Save existing overlap regions before we overwrite
                    saved = {}
                    if do_register:
                        if dx > 0 and overlap_x >= 8:
                            ow = min(overlap_x, tw, dx)
                            if ow >= 4:
                                saved['h'] = (
                                    ow,
                                    mm[ch_idx, dy:dy+th, dx:dx+ow].copy(),
                                )
                        if dy > 0 and overlap_y >= 8:
                            oh = min(overlap_y, th, dy)
                            if oh >= 4:
                                saved['v'] = (
                                    oh,
                                    mm[ch_idx, dy:dy+oh, dx:dx+tw].copy(),
                                )

                    # Write full tile (overlap regions will be re-blended)
                    mm[ch_idx, dy:dy+th, dx:dx+tw] = tile_ch

                    # Blend horizontal overlap
                    if 'h' in saved:
                        ow, exist = saved['h']
                        blended = _blend_horiz(exist, tile_ch[:, :ow])
                        mm[ch_idx, dy:dy+th, dx:dx+ow] = blended

                    # Blend vertical overlap
                    if 'v' in saved:
                        oh, exist = saved['v']
                        blended = _blend_vert(exist, tile_ch[:oh, :])
                        mm[ch_idx, dy:dy+oh, dx:dx+tw] = blended

        mm.flush()

        pyr_arrays = []
        if n_levels > 1:
            src = mm
            for level in range(1, n_levels):
                pyr_h, pyr_w = src.shape[1] // 2, src.shape[2] // 2
                print(f'  Pyramid level {level}: Y={pyr_h}, X={pyr_w}')
                pyr_arr = _downsample(src, C)
                pyr_arrays.append(pyr_arr)
                src = pyr_arr

        pixel_size = meta['zoom_info'][0]['px_x']
        ome_meta = _ome_metadata(ch_list, pixel_size)

        with tifffile.TiffWriter(output_path, bigtiff=True) as tif:
            tif.write(
                data=mm,
                tile=(1024, 1024),
                subifds=n_levels - 1,
                photometric='minisblack',
                compression='zlib',
                compressionargs={'level': 6},
                metadata=ome_meta,
            )
            for pyr_arr in pyr_arrays:
                tif.write(
                    data=pyr_arr,
                    tile=(1024, 1024),
                    subfiletype=1,
                    photometric='minisblack',
                    compression='zlib',
                    compressionargs={'level': 1},
                )

        size_gb = output_path.stat().st_size / 1e9
        print(f'  Written: {output_path} ({size_gb:.1f} GB)')
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)

    return True
