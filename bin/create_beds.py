#!/usr/bin/env python
import argparse
import pickle

from primer_explorer.bed import write_bed
from primer_explorer.regions import GenomeRegions


def get_pcr_products_sets(pcr_products_fhand):
    for pcr_products_set in pcr_products_fhand:
        primers, pcr_products = pcr_products_set
        yield primers, pcr_products


def parse_arguments():
    parser = argparse.ArgumentParser(description="Create bed files for pcr_products")
    parser.add_argument('-p', '--pcr_products', help='path to pcr products',
                        type=argparse.FileType('rb'), required=True)
    parser.add_argument('-o', '--out_bed', help='directory to output beds',
                        required=True)
    return parser


def get_args():
    parser = parse_arguments()
    args = parser.parse_args()
    pcr_products_fpath = args.pcr_products
    out_fpath = args.out_bed
    return {'pcr_products_fpath': pcr_products_fpath,
            'out_fpath': out_fpath}


def main():
    args = get_args()
    pcr_products_sets = pickle.load(args['pcr_products_fpath'])
    out_fdir = args['out_fpath']
    for pcr_products_set in pcr_products_sets:
        for pair, pcr_products in pcr_products_set['products'].items():
            pcr_products_regions = GenomeRegions(pcr_products=pcr_products)
            write_bed(pair, pcr_products_regions, out_fdir)


if __name__ == "__main__":
    main()
