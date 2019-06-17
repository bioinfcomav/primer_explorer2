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
    msg = 'Size of the read to generate regions. Depends on sequence technology you want to predict'
    parser.add_argument('-r', '--read_size', help=msg, required=True,
                        type=int)
    return parser


def get_args():
    parser = parse_arguments()
    args = parser.parse_args()
    pcr_products_fhand = args.pcr_products
    out_dir = args.out_dir
    read_length = args.read_size
    return {'pcr_products_fhand': pcr_products_fhand, 'out_dir': out_dir,
            'read_length': read_length}


def main():
    args = get_args()
    pcr_products_sets = pickle.load(args['pcr_products_fhand'])
    out_dir = Path(args['out_dir'])
    read_length = args['read_length']
    if not out_dir.exists():
        out_dir.mkdir(exist_ok=True)

    write_primer_regions_in_bed_format(pcr_products_sets, out_dir,
                                       min_product_length=100,
                                       max_product_length=700,
                                       read_length=read_length)


if __name__ == "__main__":
    main()
