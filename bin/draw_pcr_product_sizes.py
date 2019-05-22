#!/usr/bin/env python
import argparse
import pickle
import itertools

from primer_explorer.plot import draw_histograms, IntCounter

BINS = 10


def _calc_product_length(pair):
    if pair[0].chrom_location[0] != pair[1].chrom_location[0]:
        raise RuntimeError("the products doesn't have the same chrom")
    start = pair[0].chrom_location[1]
    end = pair[1].chrom_location[1]
    return abs(start - end)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot pcr lengths distributions")
    parser.add_argument('-l', '--locations', help='path to pcr products',
                        type=argparse.FileType('rb'))
    parser.add_argument('-n', '--num_sets', help='number of sets to look',
                        type=int, default=1)
    parser.add_argument('-b', '--num_bin', help='number of bins to draw',
                        type=int, default=BINS)
    parser.add_argument('-o', '--output', required=True,
                        help='Path to write the output', type=argparse.FileType('wb'))
    parser.add_argument('-p', '--primers',
                        help='use only the primers of this file',
                        required=False, type=argparse.FileType('rb'))
    return parser


def get_args():
    parser = parse_arguments()
    args = parser.parse_args()
    location_fhand = args.locations
    num_sets = args.num_sets
    num_bin = args.num_bin
    output_fhand = args.output

    primers_fhand = args.primers
    return {'locations': location_fhand, 'num_sets': num_sets,
            'num_bin': num_bin, 'output': output_fhand,
            'primers_fhand': primers_fhand}


def main():
    arguments = get_args()
    products_fhand = arguments["locations"]
    output_fhand = arguments['output']
    primers_fhand = arguments['primers_fhand']
    if primers_fhand:
        selected_primers = [line.strip() for line in primers_fhand]
    else:
        selected_primers = None

    num_sets_to_represent = arguments['num_sets']
    num_bins = arguments['num_bin']

    pcr_products = pickle.load(products_fhand)

    stats = get_stats(pcr_products, selected_primers, num_sets_to_represent)

    titles = []
    counters = []

    for pair, lengths in stats.items():
        titles.append("PCR length products for pair {}, {}".format(pair[0], pair[1]))
        counters.append(IntCounter(lengths))

    draw_histograms(counters, output_fhand, xlabel="Product_size",
                    ylabel="Number of products", titles=titles,
                    plots_per_chart=1, kind='bar', num_bins=num_bins)


def get_stats(pcr_products, selected_primers, num_sets_to_represent):
    stats = {}
    for set_index in range(num_sets_to_represent):
        primer_set = pcr_products[set_index]

        if selected_primers is not None:
            combinations = itertools.combinations(selected_primers, 2)
        else:
            combinations = primer_set.keys()

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
                lengths = [_calc_product_length(pair) for pair in products_to_calculate]
                stats[combination] = lengths

    return stats


if __name__ == '__main__':
    main()
