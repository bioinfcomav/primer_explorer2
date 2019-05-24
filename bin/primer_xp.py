#!/usr/bin/env python
import pickle
import argparse

from pathlib import Path

from primer_explorer.utils import get_fhand
from primer_explorer.kmer import get_kmers, count_kmers
from primer_explorer.pcr import (select_primers_combinations, get_pcr_products,
                                 annotate_products)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Look for working primers")
    parser.add_argument('-f', '--genome_path',
                        help='Path to the genome in fasta format',
                        type=argparse.FileType('rb'), required=True)
    parser.add_argument('-r', '--heterochromatin_regions',
                        help='Path to bed file with heterochromatin',
                        type=argparse.FileType('rb'))
    parser.add_argument('-s', '--kmer_size', type=int, default=8,
                        help='Size of the kmers to look for')
    parser.add_argument('-c', '--cache_dir', help='cache dir',
                        default='./cache')
    parser.add_argument('-o', '--pcr_products', required=True,
                        help='Path to write pcr_products result',
                        type=argparse.FileType('wb'))
    return parser


def get_args():
    parser = parse_arguments()
    args = parser.parse_args()
    genome_fhand = args.genome_path
    regions_fhand = args.heterochromatin_regions
    kmer_size = args.kmer_size
    cache_dir = args.cache_dir
    products_fhand = args.pcr_products

    return {'genome_fhand': genome_fhand, 'regions_fhand': regions_fhand,
            'kmer_size': kmer_size, 'cache_dir': cache_dir,
            'products_fhand': products_fhand}


def main():
    args = get_args()
    genome_fhand = args['genome_fhand']
    heterochromatic_regions_fhand = args['regions_fhand']
    kmer_len = args['kmer_size']
    cache_dir = Path(args['cache_dir'])
    if not cache_dir.exists():
        cache_dir.mkdir(exist_ok=True)

    pcr_products_fhand = args['products_fhand']

    genome_fhand = get_fhand(genome_fhand.name)
    regions_fhand = get_fhand(heterochromatic_regions_fhand.name)

    counters = count_kmers(genome_fhand, kmer_len,
                           regions_fhand=regions_fhand)

    most_freq_kmers = [k[0] for k in counters[False].most_common(1000)]

    kmers, kmers_locations = get_kmers(genome_fhand.name,
                                       heterochromatic_regions_fhand.name,
                                       kmer_len, cache_dir,
                                       kmers_to_keep=most_freq_kmers)

    primer_combinations = select_primers_combinations(kmers)
    product_results = []

    for primer_combination in primer_combinations:
        pcr_products = get_pcr_products(kmers_locations, primer_combination)
        product_results.append(pcr_products)

    pickle.dump(product_results, pcr_products_fhand, pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    main()
