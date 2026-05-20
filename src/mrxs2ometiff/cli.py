import argparse
from pathlib import Path

from .reader import parse_ini, read_records, read_tile_positions
from .convert import convert_one
from .verify import verify_tiff


def main():
    parser = argparse.ArgumentParser(description='Convert MRXS whole-slide images to 16-bit OME-TIFF')
    parser.add_argument('input', help='MRXS file or directory')
    parser.add_argument('-o', '--output', help='Output OME-TIFF path (default: input name in TIFF/)')
    parser.add_argument('--verify', action='store_true', help='Verify existing TIFF')
    parser.add_argument('--no-pyramid', action='store_true', help='Skip pyramid sub-IFDs')
    args = parser.parse_args()

    pyramid = not args.no_pyramid
    input_path = Path(args.input)

    if args.verify:
        if input_path.is_dir():
            for f in sorted(input_path.glob('*.ome.tif')):
                verify_tiff(f)
        else:
            verify_tiff(input_path)
        return

    if input_path.is_dir():
        mrxs_dir = input_path
        output_dir = mrxs_dir.parent / 'TIFF'
        output_dir.mkdir(parents=True, exist_ok=True)
        for f in sorted(mrxs_dir.glob('*.mrxs')):
            out = output_dir / f.with_suffix('.ome.tif').name
            convert_one(f, out, pyramid=pyramid)
    else:
        if args.output:
            convert_one(input_path, Path(args.output), pyramid=pyramid)
        else:
            output_dir = input_path.parent.parent / 'TIFF'
            output_dir.mkdir(parents=True, exist_ok=True)
            out = output_dir / input_path.with_suffix('.ome.tif').name
            convert_one(input_path, out, pyramid=pyramid)
