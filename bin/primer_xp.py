#!/usr/bin/env python
import pickle
import argparse

from primer_explorer.kmer import get_kmers
from primer_explorer.pcr import (select_primers_combinations, get_pcr_products,
                                 annotate_products)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Look for working primers")
    parser.add_argument('-f', '--genome_path', help='Path to the genome in fasta format',
                        type=argparse.FileType('rb'), required=True)
    parser.add_argument('-r', '--heterochromatin_regions',
                        help='Path to bed file with heterochromatin',
                        type=argparse.FileType('rb'))
    parser.add_argument('-s', '--kmer_size', help='Size of the kmers to look for',
                        type=int, default=8)
    parser.add_argument('-c', '--cache_dir', help='cache dir',
                        default='./cache')
    parser.add_argument('-o', '--pcr_products', required=True,
                        help='Path to write pcr_products result',
                        type=argparse.FileType('wb'))
    parser.add_argument('-a', '--pcr_annotations', required=True,
                        help='Path to write pcr annotations result',
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
    annotations_fhand = args.pcr_annotations

    return {'genome_fhand': genome_fhand, 'regions_fhand': regions_fhand,
            'kmer_size': kmer_size, 'cache_dir': cache_dir,
            'products_fhand': products_fhand,
            'annotations_fhand': annotations_fhand}


def main():
    args = get_args()
    genome_fpath = args['genome_fhand'].name
    heterochromatic_regions_fpath = args['regions_fhand'].name
    kmer_len = args['kmer_size']
    cache_dir = args['cache_dir']
    pcr_products_fhand = args['products_fhand']
    pcr_annotation_fhand = args['annotations_fhand']

    kmers, kmers_locations = get_kmers(genome_fpath, heterochromatic_regions_fpath, kmer_len, cache_dir)

    primer_combinations = select_primers_combinations(kmers)
    product_results = []
    annotation_results = []

    for primer_combination in primer_combinations:
        pcr_products = get_pcr_products(kmers_locations, primer_combination)
        product_results.append(pcr_products)
        annotation = annotate_products(pcr_products)
        annotation_results.append(annotation)

    pickle.dump(product_results, pcr_products_fhand, pickle.HIGHEST_PROTOCOL)
    pickle.dump(annotation_results, pcr_annotation_fhand, pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    main()
