#!/usr/bin/env python
import argparse
import pickle

from pathlib import Path
from primer_explorer.regions import write_primer_regions_in_bed_format


def get_pcr_products_sets(pcr_products_fhand):
    for pcr_products_set in pcr_products_fhand:
        primers, pcr_products = pcr_products_set
        yield primers, pcr_products


def parse_arguments():
    parser = argparse.ArgumentParser(description="Create bed files for pcr_products")
    parser.add_argument('-p', '--pcr_products', help='path to pcr products',
                        type=argparse.FileType('rb'), required=True)
    parser.add_argument('-o', '--out_dir', help='directory to output beds',
                        required=True)
    return parser


def get_args():
    parser = parse_arguments()
    args = parser.parse_args()
    pcr_products_fpath = args.pcr_products
    out_dir = args.out_dir
    return {'pcr_products_fpath': pcr_products_fpath, 'out_dir': out_dir}


def main():
    args = get_args()
    pcr_products_sets = pickle.load(args['pcr_products_fpath'])
    out_dir = Path(args['out_dir'])
    if not out_dir.exists():
        out_dir.mkdir(exist_ok=True)

    write_primer_regions_in_bed_format(pcr_products_sets, out_dir)


if __name__ == "__main__":
    main()
