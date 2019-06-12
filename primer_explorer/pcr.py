from collections import namedtuple
from itertools import combinations

from primer_explorer.primer3.primer3 import reverse_complement
from primer_explorer.combinations import get_compatible_groups_of_primers
from primer_explorer.stats import filter_by_length

PrimerAndLocation = namedtuple('PrimerAndLocation', ('strand', 'seq', 'chrom_location', 'is_heterochromatic'))


def select_primers_combinations(kmers, num_compatible_groups=10):
    compatibility_primers_groups = get_compatible_groups_of_primers(kmers, num_compatible_groups=num_compatible_groups)
    return compatibility_primers_groups


def _get_pcr_products(fwd_primers_locations, rev_primers_locations, max_pcr_product_length):
    fwd_locs = [PrimerAndLocation(0, loc.seq, loc.chrom_location, loc.is_heterochromatic) for loc in fwd_primers_locations]
    locs = fwd_locs
    rev_locs = [PrimerAndLocation(1, loc.seq, loc.chrom_location, loc.is_heterochromatic) for loc in rev_primers_locations]
    locs.extend(rev_locs)
    locs = sorted(locs, key=lambda PrimerAndLocation: (PrimerAndLocation.chrom_location))
    pcr_products = [(loc1, loc2) for loc1, loc2 in zip(locs[:-1], locs[1:]) if (loc1.strand == 0) and (loc1.strand != loc2.strand) and (loc1.chrom_location[0] == loc2.chrom_location[0]) and abs(loc1.chrom_location[1] - loc2.chrom_location[1]) < max_pcr_product_length]
    return pcr_products


def get_pcr_products(kmer_locations, primers, max_pcr_product_length=10000):
    pcr_products = {}
    for combination in combinations(primers, 2):
        fwd_primers = combination
        fwd_primers_locations = []
        rev_primers = [reverse_complement(primer) for primer in combination]
        rev_primers_locations = []

        for primer in fwd_primers:
            fwd_primers_locations.extend(kmer_locations.get(primer, ""))

        for primer in rev_primers:
            rev_primers_locations.extend(kmer_locations.get(primer, ""))
        pcr_products[combination] = _get_pcr_products(fwd_primers_locations,
                                                      rev_primers_locations,
                                                      max_pcr_product_length=max_pcr_product_length)
    return pcr_products


def count_total_products(pcr_products):
    return len(pcr_products)


def get_pcr_products_in_sets(primer_combinations, kmers, kmers_locations,
                             max_pcr_product_length):
    product_results = []

    for primer_combination in primer_combinations:
        pcr_products = get_pcr_products(kmers_locations, primer_combination,
                                        max_pcr_product_length=max_pcr_product_length)
        primer_group = {'primers': primer_combination,
                        'products': pcr_products}
        product_results.append(primer_group)
    return product_results


def pcr_products_to_regions(pcr_products, min_length, max_length, read_length):
    pcr_products = filter_by_length(pcr_products, min_length=min_length, max_length=max_length)
    for pcr_product in pcr_products:
        chrom_fwd = pcr_product[0].chrom_location
        chrom_rev = pcr_product[1].chrom_location
        yield "{}\t{}\t{}".format(chrom_fwd[0].decode(), chrom_fwd[1], chrom_fwd[1] + read_length)
        yield "{}\t{}\t{}".format(chrom_rev[0].decode(), chrom_rev[1] - read_length, chrom_rev[1])
