#!/usr/bin/env python
import argparse
import pickle

from primer_explorer.stats import get_product_lengths_by_pair
from primer_explorer.plot import plot_length_distribs


def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot pcr lengths distributions")
    parser.add_argument('-p', '--pcr_products', help='path to pcr products',
                        type=argparse.FileType('rb'))
    parser.add_argument('-n', '--num_sets', help='number of sets to look',
                        type=int, default=1)
    parser.add_argument('-o', '--output', required=True,
                        help='Path to write the output',
                        type=argparse.FileType('wb'))
    parser.add_argument('-r', '--primers',
                        help='use only the primers of this file',
                        required=False, type=argparse.FileType('rb'))
    return parser


def get_args():
    parser = parse_arguments()
    args = parser.parse_args()
    pcr_products_fhand = args.pcr_products
    num_sets = args.num_sets
    output_fhand = args.output

    primers_fhand = args.primers
    return {'pcr_products': pcr_products_fhand,
            'num_sets': num_sets, 'output': output_fhand,
            'primers_fhand': primers_fhand}


def main():
    arguments = get_args()
    pcr_products_fhand = arguments["pcr_products"]
    out_fhand = arguments['output']
    primers_fhand = arguments['primers_fhand']
    if primers_fhand:
        selected_primers = [line.strip() for line in primers_fhand]
    else:
        selected_primers = None

    num_sets_to_represent = arguments['num_sets']

    pcr_products = pickle.load(pcr_products_fhand)

    stats = get_product_lengths_by_pair(pcr_products, selected_primers,
                                        num_sets_to_represent)

    plot_length_distribs(stats, out_fhand)


if __name__ == '__main__':
    main()
