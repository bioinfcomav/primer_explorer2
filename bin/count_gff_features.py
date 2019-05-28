#!/usr/bin/env python
import argparse
from collections import OrderedDict
from os import listdir
from os.path import join
from primer_explorer.gff import get_gff_intersect_results
from primer_explorer.report import write_gff_report


def get_beds_fpaths(bed_dir):
    bed_fpaths = OrderedDict()
    for bed_fname in listdir(path=bed_dir):
        if bed_fname.endswith(".bed"):
            primers_pair = str(bed_fname).split(".")[0]
            bed_fpaths[primers_pair] = join(bed_dir, bed_fname)
    return bed_fpaths


def parse_arguments():
    parser = argparse.ArgumentParser(description="Count gff features from bed_files")
    parser.add_argument('-g', '--gff_fpath', help='gff file', required=True)
    parser.add_argument('-b', '--bed_fdir', help='directory to input beds',
                        required=True)
    parser.add_argument('-o', '--output_fpath', help='report output',
                        type=argparse.FileType('wb'),
                        required=True)

    return parser


def get_args():
    parser = parse_arguments()
    args = parser.parse_args()
    gff_fpath = args.gff_fpath
    bed_fdir = args.bed_fdir
    out_fpath = args.output_fpath
    return {'gff_fpath': gff_fpath, 'bed_fdir': bed_fdir,
            'output_fpath': out_fpath}


def main():
    args = get_args()
    bed_fdir = args['bed_fdir']
    gff_fpath = args['gff_fpath']
    out_fpath = args['output_fpath']
    beds_fpaths = get_beds_fpaths(bed_fdir)
    gff_results = get_gff_intersect_results(beds_fpaths, gff_fpath)
    write_gff_report(gff_results, out_fpath)


if __name__ == "__main__":
    main()
