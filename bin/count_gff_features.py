#!/usr/bin/env python
import argparse
from collections import OrderedDict
from os import listdir
from os.path import join

from primer_explorer.gff import count_regions_in_gff_features
from primer_explorer.report import write_gff_report
from primer_explorer import config


def parse_arguments():
    desc = "Count gff features from bed_files"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-g', '--gff_fpath', help='gff file', required=True)
    parser.add_argument('-b', '--bed_fdir', help='directory to input beds',
                        required=True)
    parser.add_argument('-o', '--output_fpath', help='report output',
                        type=argparse.FileType('wb'),
                        required=True)
    parser.add_argument('-f', '--features', help='Features you want to count',
                        default=config.GFF_FEATURES, nargs='*')
    return parser


def get_args():
    parser = parse_arguments()
    args = parser.parse_args()
    gff_fpath = args.gff_fpath
    bed_fdir = args.bed_fdir
    out_fpath = args.output_fpath
    features = args.features
    return {'gff_fpath': gff_fpath, 'bed_fdir': bed_fdir,
            'output_fpath': out_fpath, 'features': features}


def get_beds_fpaths(bed_dir):
    bed_fpaths = OrderedDict()
    for bed_fname in listdir(path=bed_dir):
        if bed_fname.endswith(".bed"):
            primers_pair = str(bed_fname).split(".")[0]
            bed_fpaths[primers_pair] = join(bed_dir, bed_fname)
    return bed_fpaths


def _count_lines(fpath):
    f = open(fpath, 'rb')
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.raw.read

    buf = read_f(buf_size)
    while buf:
        lines += buf.count(b'\n')
        buf = read_f(buf_size)

    return lines


def main():
    args = get_args()
    bed_fdir = args['bed_fdir']
    gff_fpath = args['gff_fpath']
    out_fpath = args['output_fpath']
    features_to_count = args['features']

    fpaths_by_pair = get_beds_fpaths(bed_fdir)
    counts_by_pair = {}
    for pair, fpath in fpaths_by_pair.items():
        counts_by_pair[pair] = {}
        counts_by_pair[pair]["num_pcr_products"] = _count_lines(fpath)
        feature_counts = count_regions_in_gff_features(fpath, gff_fpath,
                                                       features_to_count=features_to_count)
        counts_by_pair[pair].update(feature_counts)

    write_gff_report(counts_by_pair, out_fpath)


if __name__ == "__main__":
    main()
