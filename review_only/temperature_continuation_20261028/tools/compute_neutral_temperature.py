"""Compute streamed local neutral kinetic temperatures from saved packets."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.meta_path[:] = [f for f in sys.meta_path if type(f).__module__!='_editable_skbc_mhdlab']
from mhdlab.neutral_kinetic_maps import compute_neutral_temperature


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packets',type=Path,required=True)
    parser.add_argument('--maps',type=Path,required=True)
    parser.add_argument('--metadata',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--coarsen-factor',type=int,default=1)
    parser.add_argument('--minimum-effective-sample-size',type=float,default=50.)
    parser.add_argument('--minimum-expected-cell-particles',type=float,default=1.)
    args=parser.parse_args()
    try:
        report=compute_neutral_temperature(args.packets,args.maps,args.metadata,args.output_dir,
            coarsen_factor=args.coarsen_factor,minimum_effective_sample_size=args.minimum_effective_sample_size,
            minimum_expected_cell_particles=args.minimum_expected_cell_particles)
    except (OSError,ValueError,KeyError,TypeError) as exc:
        parser.exit(2,f'neutral kinetic temperature: {exc}\n')
    print(json.dumps({'hdf':str(args.output_dir/report['hdf_file']),'shape':report['shape'],'status':report['status']}))


if __name__=='__main__': main()
