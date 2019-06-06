#!/usr/bin/env python
import argparse
import pickle
import itertools

from primer_explorer.plot import draw_histograms
from primer_explorer.stats import IntCounter

DEFAULT_NUMBER_OF_BINS = 20


def _calc_product_length(pair):
    if pair[0].chrom_location[0] != pair[1].chrom_location[0]:
        raise RuntimeError("the products doesn't have the same chrom")
    start = pair[0].chrom_location[1]
    end = pair[1].chrom_location[1]
    return abs(start - end)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot pcr lengths distributions")
    parser.add_argument('-p', '--pcr_products', help='path to pcr products',
                        type=argparse.FileType('rb'))
    parser.add_argument('-n', '--num_sets', help='number of sets to look',
                        type=int, default=1)
    parser.add_argument('-b', '--num_bin', help='number of bins to draw',
                        type=int, default=DEFAULT_NUMBER_OF_BINS)
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
    num_bin = args.num_bin
    output_fhand = args.output

    primers_fhand = args.primers
    return {'pcr_products': pcr_products_fhand, 'num_sets': num_sets,
            'num_bin': num_bin, 'output': output_fhand,
            'primers_fhand': primers_fhand}


def main():
    arguments = get_args()
    pcr_products_fhand = arguments["pcr_products"]
    output_fhand = arguments['output']
    primers_fhand = arguments['primers_fhand']
    if primers_fhand:
        selected_primers = [line.strip() for line in primers_fhand]
    else:
        selected_primers = None

    num_sets_to_represent = arguments['num_sets']
    num_bins = arguments['num_bin']

    pcr_products = pickle.load(pcr_products_fhand)

    stats = get_stats(pcr_products, selected_primers, num_sets_to_represent)
    titles = []
    counters = []

    for pair, lengths in stats.items():
        title = "PCR length products for pair {}, {} (all)"
        titles.append(title.format(pair[0].decode(), pair[1].decode()))
        counters.append(IntCounter(lengths))

        title = "PCR length products for pair {}, {} (max 1000 lenth)"
        titles.append(title.format(pair[0].decode(), pair[1].decode()))
        lengths = [len_ for len_ in lengths if len_ <= 1000]
        counters.append(IntCounter(lengths))
    # xlabel="Product_size"
    xlabel = None
    draw_histograms(counters, output_fhand, xlabel=xlabel,
                    titles=titles, ylabel="Number of products",
                    plots_per_chart=1, kind='bar', num_bins=num_bins,
                    xtickslabel_rotation=60, vlines=100)


def get_stats(pcr_products, selected_primers, num_sets_to_represent):
    stats = {}

    for set_index in range(num_sets_to_represent):
        primer_set = pcr_products[set_index]['products']
        if selected_primers is not None:
            combinations = itertools.combinations(selected_primers, 2)
        else:
            combinations = primer_set.keys()
        combinations = sorted(combinations)
        for combination in combinations:
            if combination in stats:
                continue
            products = primer_set.get(combination, None)
            products2 = primer_set.get((combination[1], combination[0], None))
            products_to_calculate = []

            if products is not None:
                products_to_calculate += products
            if products2 is not None:
                products_to_calculate += products2

            if products_to_calculate:
                lengths = [_calc_product_length(pair) for pair in primer_set[combination]]
                stats[combination] = lengths

    return stats


if __name__ == '__main__':
    main()
