import pickle
import sys

from primer_explorer import config
from primer_explorer.kmer import get_kmers
from primer_explorer.pcr import (select_primers_combinations, get_pcr_products,
                                 annotate_products)

TOMATO_CHROM_FASTA_GZ = '/home/jope/devel3/primer_explorer_old_version/genome/S_lycopersicum_chromosomes.3.00.chrom1.fasta.gz'
HETEROCHROMATIN_BED = '/home/jope/devel3/primer_experience/test_heterochromatin.bed'

if __name__ == '__main__':
    genome_fpath = "/home/jope/genomes/SL3.0/S_lycopersicum_chromosomes.3.00.fa"
    heterochromatic_regions_fpath = "test_heterochromatin.bed"
    kmer_len = 8
    cache_dir = config.CACHE_DIR
    kmers, kmers_locations = get_kmers(genome_fpath, heterochromatic_regions_fpath, kmer_len, cache_dir)

    # # heterochromatic_regions = GenomeRegions(open("/home/jope/genomes/SL3.0/ITAG3.2_REPET_repeats_agressive.gff", 'rb'))
    # heterochromatic_regions = GenomeRegions(open("test_heterochromatin.bed", 'rb'))
    # # genome_fasta_fhand = io.BytesIO(GENOME_FASTA)
    # # genome_fasta_fhand = io.StringIO(GENOME_FASTA.decode())
    # # genome_fasta_fhand = gzip.open(TOMATO_CHROM_FASTA_GZ, 'rb')
    # # genome_fasta_fhand = open('chrom1.100k.fasta', 'rb')
    # genome_fasta_fhand = open("/home/jope/genomes/SL3.0/S_lycopersicum_chromosomes.3.00.fa", 'rb')

    # kmers, kmers_locations = generate_kmer_locations(genome_fasta_fhand, kmer_len=8,
    #                                                  heterochromatic_regions=heterochromatic_regions)

    sys.exit()
    for kmer in kmers:
        print(kmer.decode())
    primer_combinations = select_primers_combinations(kmers)
    product_results = []
    annotation_results = []
    for primer_combination in primer_combinations:
        pcr_products = get_pcr_products(kmers_locations, primer_combination)
        product_results.append(pcr_products)
        annotation = annotate_products(pcr_products)
        annotation_results.append(annotation)
    pcr_products_fhand = open(config.PCR_PRODUCTS_PICKLE, "wb")
    pcr_annotation_fhand = open(config.PCR_ANNOTATION, "wb")
    pickle.dump(product_results, pcr_products_fhand, pickle.HIGHEST_PROTOCOL)
    pickle.dump(annotation_results, pcr_annotation_fhand, pickle.HIGHEST_PROTOCOL)
