from collections import namedtuple, Counter, defaultdict
from itertools import combinations

from primer_explorer.primer3.primer3 import reverse_complement
from primer_explorer import config
from primer_explorer.combinations import get_compatible_groups_of_primers

PrimerAndLocation = namedtuple('PrimerAndLocation', ('strand', 'seq', 'chrom_location', 'is_heterochromatic'))


def select_primers_combinations(kmers, num_possible_combinations=10):
    compatibility_primers_groups = get_compatible_groups_of_primers(kmers, num_compatible_groups=num_possible_combinations)
    return compatibility_primers_groups


def _get_pcr_products(fwd_primers_locations, rev_primers_locations, max_pcr_product_length=1000):
    fwd_locs = [PrimerAndLocation(0, loc.seq, loc.chrom_location, loc.is_heterochromatic) for loc in fwd_primers_locations]
    locs = fwd_locs
    rev_locs = [PrimerAndLocation(1, loc.seq, loc.chrom_location, loc.is_heterochromatic) for loc in rev_primers_locations]
    locs.extend(rev_locs)
    locs = sorted(locs, key=lambda PrimerAndLocation: (PrimerAndLocation.chrom_location))
    pcr_products = [(loc1, loc2) for loc1, loc2 in zip(locs[:-1], locs[1:]) if (loc1.strand != loc2.strand) and (loc1.chrom_location[0] == loc2.chrom_location[0]) and abs(loc1.chrom_location[1] - loc2.chrom_location[1]) < max_pcr_product_length]
    return pcr_products


def get_pcr_products(kmer_locations, primers):
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
        pcr_products[combination] = _get_pcr_products(fwd_primers_locations, rev_primers_locations)
    return pcr_products


def annotate_reverse_complement_dimers(pcr_products):
    rev_comp_dimers_counts = Counter()
    for product in pcr_products:
        fwd_primer = product[0].seq
        rev_primer = product[1].seq
        if fwd_primer == reverse_complement(rev_primer):
            rev_comp_dimers_counts[True] += 1
        else:
            rev_comp_dimers_counts[False] += 1
    return rev_comp_dimers_counts


def annotate_length_viable_products(pcr_products,
                                    min_viable_length=config.MIN_PCR_VIABLE_LENGTH):
    length_viable_counts = Counter()
    for product in pcr_products:
        fwd_primer_position = product[0].chrom_location[1]
        rev_primer_position = product[1].chrom_location[1]
        if abs(fwd_primer_position - rev_primer_position) < min_viable_length:
            length_viable_counts["inviable"] += 1
        else:
            length_viable_counts["viable"] += 1
    return length_viable_counts


def annotate_products_by_euchromatin_region(pcr_products):
    euchromatin_locations = Counter()
    for product in pcr_products:
        if all(primer.is_heterochromatic for primer in product):
            euchromatin_locations["heterochromatic"] += 1
        elif all(not primer.is_heterochromatic for primer in product):
            euchromatin_locations["euchromatic"] += 1
        else:
            euchromatin_locations["mixed"] += 1
    return euchromatin_locations


def count_total_products(pcr_products):
    return len(pcr_products)


def get_pcr_products_in_sets(primer_combinations, kmers, kmers_locations):
    product_results = []

    for primer_combination in primer_combinations:

        pcr_products = get_pcr_products(kmers_locations, primer_combination)
        primer_group = {'primers': primer_combination,
                        'products': pcr_products}
        product_results.append(primer_group)
    return product_results


def annotate_products(pcr_products):
    annotations = defaultdict(dict)
    for primers, products in pcr_products.items():
        annotations[primers]["total_products"] = count_total_products(products)
        annotations[primers]["pcr_revcomp_dimers"] = annotate_reverse_complement_dimers(products)
        annotations[primers]["length_viable_products"] = annotate_length_viable_products(products)
        annotations[primers]["euchromatin_products"] = annotate_products_by_euchromatin_region(products)
    return annotations
