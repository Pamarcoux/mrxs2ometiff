import xml.etree.ElementTree as ET

import tifffile

OME_NS = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}


def verify_tiff(tiff_path):
    """Verify OME-TIFF structure and per-channel pixel ranges."""
    print(f'\nVerifying: {tiff_path.name}')
    with tifffile.TiffFile(tiff_path) as tif:
        print(f'  BigTIFF: {tif.is_bigtiff}')
        print(f'  Series: {len(tif.series)}')
        for si, s in enumerate(tif.series):
            print(f'    Series {si}: shape={s.shape}, dtype={s.dtype}')

        series = tif.series[0]
        stack = series.asarray()
        for ci in range(stack.shape[0]):
            ch_data = stack[ci]
            print(f'    Channel {ci}: min={ch_data.min()}, max={ch_data.max()}, mean={ch_data.mean():.1f}')

        root = ET.fromstring(tif.ome_metadata)
        pixels = root.find('.//ome:Pixels', OME_NS)
        if pixels is not None:
            print(f'  OME: SizeC={pixels.get("SizeC")}, SizeX={pixels.get("SizeX")}, SizeY={pixels.get("SizeY")}')
            for ch_el in pixels.findall('ome:Channel', OME_NS):
                name = ch_el.get('Name', '<unnamed>')
                ex = ch_el.get('ExcitationWavelength', '')
                em = ch_el.get('EmissionWavelength', '')
                print(f'    Channel: {name}  ex={ex}  em={em}')
