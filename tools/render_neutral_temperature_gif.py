#!/usr/bin/env python3
"""Render saved centered neutral-velocity temperatures without evolving them."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.meta_path[:]=[finder for finder in sys.meta_path
                 if type(finder).__module__!='_editable_skbc_mhdlab']
from mhdlab.neutral_temperature_render import render_temperature


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',required=True,type=Path,help='saved kinetic-moment HDF5')
    parser.add_argument('--metadata',required=True,type=Path)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--fps',type=float,default=12.)
    parser.add_argument('--dpi',type=int,default=110)
    parser.add_argument('--temperature-max-k',type=float,help='explicit shared color maximum; input values stay unchanged')
    parser.add_argument('--view-bounds-mm',type=float,nargs=4,metavar=('XMIN','XMAX','YMIN','YMAX'))
    args=parser.parse_args()
    try:
        report=render_temperature(args.input,args.metadata,args.output_dir,fps=args.fps,dpi=args.dpi,
            temperature_max_k=args.temperature_max_k,
            view_bounds_m=None if args.view_bounds_mm is None else [v/1000 for v in args.view_bounds_mm])
    except (OSError,ValueError,TypeError,KeyError) as error:
        parser.exit(2,f'neutral temperature GIF: {error}\n')
    print(json.dumps({'gif':str(args.output_dir/'neutral_temperature.gif'),
                     'frames':report['frame_count'],'status':report['status']}))


if __name__=='__main__': main()
