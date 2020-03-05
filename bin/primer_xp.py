#!/usr/bin/env python
import pickle
import argparse

from pathlib import Path

from primer_explorer.kmer import get_kmers
from primer_explorer.pcr import (select_primers_combinations,
                                 get_pcr_products_in_sets)
from primer_explorer.stats import get_stats_by_pair_in_sets
from primer_explorer.report import write_stats_in_excel
from primer_explorer.primer3.primer3 import reverse_complement


def parse_arguments():
    desc = "Search primers to use with K-Seq"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-f', '--genome_path',
                        help='Path to the genome in fasta format',
                        type=argparse.FileType('rb'), required=True)
    msg = 'path to a bed files with regions yu want to avoid for your primers'
    parser.add_argument('-r', '--heterochromatin_regions', help=msg,
                        type=argparse.FileType('rb'))

    parser.add_argument('-k', '--kmer_size', type=int, default=8,
                        help='Size of the kmers to look for')
    parser.add_argument('-c', '--cache_dir', help='cache dir',
                        default='./cache')
    parser.add_argument('-d', '--pcr_products', required=True,
                        help='Path to write pcr_products database in pickle format',
                        type=argparse.FileType('wb'))
    parser.add_argument('-o', '--report', required=True,
                        help='Path to write report in excel format',
                        type=argparse.FileType('wb'))
    parser.add_argument('-t', '--top_kmers', help='Max umber of best kmers to use',
                        type=int, default=1000)

    parser.add_argument('-m', '--num_sets', help="number of sets of primers to calculate",
                        type=int, default=3)

    msg = "Force to look only for this kmers"
    parser.add_argument('-l', '--forced_kmers', nargs="*", help=msg)

    return parser


def get_args():
    parser = parse_arguments()
    args = parser.parse_args()
    genome_fhand = args.genome_path
    regions_fhand = args.heterochromatin_regions
    kmer_size = args.kmer_size
    cache_dir = args.cache_dir
    products_fhand = args.pcr_products
    top_kmers = args.top_kmers
    num_sets = args.num_sets
    report_fhand = args.report
    forced_kmers = args.forced_kmers
    return {'genome_fhand': genome_fhand, 'regions_fhand': regions_fhand,
            'kmer_size': kmer_size, 'cache_dir': cache_dir,
            'products_fhand': products_fhand, 'top_kmers': top_kmers,
            'num_sets': num_sets, 'report_fhand': report_fhand,
            'forced_kmers': forced_kmers}


def main():
    args = get_args()
    genome_fhand = args['genome_fhand']
    heterochromatic_regions_fhand = args['regions_fhand']
    kmer_len = args['kmer_size']
    cache_dir = Path(args['cache_dir'])
    top_kmers = args['top_kmers']
    num_sets = args['num_sets']
    report_fhand = args['report_fhand']
    forced_kmers = args['forced_kmers']
    if not cache_dir.exists():
        cache_dir.mkdir(exist_ok=True)
    pcr_products_fhand = args['products_fhand']

    kmers_to_keep = None
    if forced_kmers:
        forced_kmers = [kmer.encode() for kmer in forced_kmers]
        kmers_to_keep = forced_kmers + [reverse_complement(kmer) for kmer in forced_kmers]

    kmers, kmers_locations = get_kmers(genome_fhand.name,
                                       heterochromatic_regions_fhand.name,
                                       kmer_len, cache_dir, num_kmers_to_keep=top_kmers,
                                       kmers_to_keep=kmers_to_keep)
    if forced_kmers:
        kmers = forced_kmers

    primer_combinations = select_primers_combinations(kmers, num_compatible_groups=num_sets)

    product_results = get_pcr_products_in_sets(primer_combinations, kmers,
                                               kmers_locations,
                                               max_pcr_product_length=10000)

    pickle.dump(product_results, pcr_products_fhand, pickle.HIGHEST_PROTOCOL)

    stats = get_stats_by_pair_in_sets(product_results,
                                      kmers_locations=kmers_locations)

    write_stats_in_excel(report_fhand.name, stats)


if __name__ == '__main__':
    main()
