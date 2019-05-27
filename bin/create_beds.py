#!/usr/bin/env python
import argparse
import pickle


from primer_explorer.bed import write_bed
from primer_explorer.pcr import PrimerAndLocation
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
    args = parser.get_args()
    pcr_products_fpath = args.pcr_products
    out_fpath = args.out_bed
    return {'pcr_products_fpath': pcr_products_fpath,
            'out_fpath': out_fpath}


def main():
    parser = get_args()
    args = parser.parse_args()
    pcr_products_fhand = pickle.load(args['pcr_products_fpath'])
    out_fdir = args['out_fpath']
    for _, primer_set in get_pcr_products_sets(pcr_products_fhand):
        for primer_pair, pcr_products in primer_set.items():
            products_regions = GenomeRegions(pcr_products)
            write_bed(primer_pair, products_regions, out_fdir)


if __name__ == "__main__":
    main()
