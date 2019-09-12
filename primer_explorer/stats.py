import itertools
from collections import Counter

from primer_explorer.primer3.primer3 import reverse_complement

LABELS = {'title': 'histogram', 'xlabel': 'values',
          'ylabel': 'count', 'minimum': 'minimum',
          'maximum': 'maximum', 'average': 'average',
          'variance': 'variance', 'sum': 'sum',
          'items': 'items', 'quartiles': 'quartiles'}

PERCENTAGE_OF_SEQUENCIABLE_NUCLEOTIDES = "percentage_sequenciable_nucleotides"
ADJUSTED_PERCENTAGE_OF_SEQUENCIABLE_NUCLEOTIDES = "adjusted_percentage_sequenciable_nucleotides"
NUM_OF_POSSIBLE_PRODUCTS_10000 = "num_possible_products_10000"
NUM_OF_POSSIBLE_PRODUCTS_700 = "num_possible_products_700"
NUM_SEQUENCIABLE_PRODUCTS = "num_sequenciable_products"
NUM_SHORT_PRODUCTS = "num_short_products"
NUM_UNION_SITES = "num_union_sites"
NUM_REPETITIVE_REPETITIVE_PRODUCTS = "num_repetitive-repetitive_products"
NUM_UNIQUE_UNIQUE_PRODUCTS = "num_unique-unique_products"
NUM_UNIQUE_REPETITIVE_PRODUCTS = "num_unique-repetitive_products"
NUM_UNION_SITES_P1 = "num_union_sites_primer_1"
NUM_UNION_SITES_P2 = "num_union_sites_primer_2"


def get_total_nondimer_pcr_products(pcr_products):
    non_dimer_pcr_products = []
    for product in pcr_products:
        fwd_primer = product[0].seq
        rev_primer = product[1].seq
        if fwd_primer == reverse_complement(rev_primer):
            continue
        else:
            non_dimer_pcr_products.append(product)
    return non_dimer_pcr_products


def filter_by_length(pcr_products, min_length=None, max_length=None):
    viable_products = []
    for product in pcr_products:
        fwd_primer_position = product[0].chrom_location[1]
        rev_primer_position = product[1].chrom_location[1]
        length = abs(fwd_primer_position - rev_primer_position)
        min_ok = False
        if (min_length is None or length >= min_length):
            min_ok = True
        max_ok = False
        if (max_length is None or length <= max_length):
            max_ok = True
        if min_ok and max_ok:
            viable_products.append(product)

    return viable_products


def get_union_sites_for_primers(pcr_products):
    union_sites = Counter()
    for product in pcr_products:
        fwd_primer = product[0].seq
        rev_primer = product[1].seq
        union_sites[fwd_primer] += 1
        union_sites[rev_primer] += 1
    return {'primer1': union_sites[fwd_primer],
            'primer2': union_sites[rev_primer]}


def get_total_euchromatic_products(pcr_products):
    euchromatic_products = []
    for product in pcr_products:
        if all(not primer.is_heterochromatic for primer in product):
            euchromatic_products.append(product)
    return euchromatic_products


def get_total_heterochromatic_products(pcr_products):
    heterochromatic_products = []
    for product in pcr_products:
        if all(primer.is_heterochromatic for primer in product):
            heterochromatic_products.append(product)
    return heterochromatic_products


def get_total_mixed_products(pcr_products):
    mixed_products = []
    for product in pcr_products:
        if (not all(not primer.is_heterochromatic for primer in product) and
                not all(primer.is_heterochromatic for primer in product)):
            mixed_products.append(product)
    return mixed_products


def get_pcr_products_counts(pcr_products, min_length, max_length, genome_length):
    pcr_products_counts = Counter()
    total_products = get_total_nondimer_pcr_products(pcr_products)
    pcr_products_counts[NUM_OF_POSSIBLE_PRODUCTS_10000] = len(total_products)

    sequenciable_products = filter_by_length(total_products, min_length, max_length)
    pcr_products_counts[NUM_SEQUENCIABLE_PRODUCTS] = {'min': min_length,
                                                      'max': max_length,
                                                      'count': len(sequenciable_products)}
    short_product_count = filter_by_length(total_products, max_length=min_length)
    pcr_products_counts[NUM_SHORT_PRODUCTS] = {'max': min_length,
                                               'count': len(short_product_count)}
    pcr_products_counts[NUM_OF_POSSIBLE_PRODUCTS_700] = len(sequenciable_products) + len(short_product_count)
    euchromatin_products = get_total_euchromatic_products(sequenciable_products)
    pcr_products_counts[NUM_UNIQUE_UNIQUE_PRODUCTS] = len(euchromatin_products)
    heterochromatin_products = get_total_heterochromatic_products(sequenciable_products)
    pcr_products_counts[NUM_REPETITIVE_REPETITIVE_PRODUCTS] = len(heterochromatin_products)
    mixed_products = get_total_mixed_products(sequenciable_products)
    pcr_products_counts[NUM_UNIQUE_REPETITIVE_PRODUCTS] = len(mixed_products)
    # TODO JUNTO BREADTH
    if genome_length:
        perc_secuenciable_nucls = get_pcr_nucleotide_count(sequenciable_products, genome_length)
        adjusted_perct_of_sequenciable_nucleotides = get_pcr_nucleotide_count(sequenciable_products,
                                                                              genome_length,
                                                                              scale_to_illumina_sequences=True)
    else:
        perc_secuenciable_nucls = 0
        adjusted_perct_of_sequenciable_nucleotides = 0

    pcr_products_counts[PERCENTAGE_OF_SEQUENCIABLE_NUCLEOTIDES] = perc_secuenciable_nucls
    pcr_products_counts[ADJUSTED_PERCENTAGE_OF_SEQUENCIABLE_NUCLEOTIDES] = adjusted_perct_of_sequenciable_nucleotides
    return pcr_products_counts


def _get_union_sites_for_kmer(kmer, kmer_locations):
    forward_unions = len(kmer_locations[kmer])
    rev_kmer = reverse_complement(kmer)
    rev_unions = len(kmer_locations[rev_kmer])

    return forward_unions + rev_unions


def get_stats_by_pair_in_sets(products_sets, min_length=100, max_length=1000,
                              genome_length=None, kmers_locations=None):
    stats = {}
    for set_index, set_info in enumerate(products_sets):
        stats[set_index] = {'primers': [p.decode() for p in set_info['primers']],
                            'stats': {}}
        for pair in sorted(set_info['products'].keys()):
            pcr_products = set_info['products'][pair]
            counts = get_pcr_products_counts(pcr_products,
                                             min_length=min_length,
                                             max_length=max_length,
                                             genome_length=genome_length)
            decoded_pair = tuple([p.decode() for p in pair])
            counts[NUM_UNION_SITES_P1] = _get_union_sites_for_kmer(pair[0], kmers_locations)
            counts[NUM_UNION_SITES_P2] = _get_union_sites_for_kmer(pair[1], kmers_locations)

            stats[set_index]['stats'][decoded_pair] = counts
    return stats


def get_pcr_nucleotide_count(pcr_products, genome_length,
                             scale_to_illumina_sequences=False,
                             illumina_min_pair_length=300):
    nucleotides = 0
    for product in pcr_products:
        start = product[0].chrom_location[1]
        end = product[1].chrom_location[1]
        distance = abs(start - end)
        if scale_to_illumina_sequences and distance > illumina_min_pair_length:
            distance = illumina_min_pair_length
        nucleotides += distance
    return float(nucleotides / genome_length)


def _calc_product_length(pair):
    if pair[0].chrom_location[0] != pair[1].chrom_location[0]:
        raise RuntimeError("the products doesn't have the same chrom")
    start = pair[0].chrom_location[1]
    end = pair[1].chrom_location[1]
    return abs(start - end)


def get_product_lengths_by_pair(pcr_products, selected_primers, num_sets_to_represent):
    stats = {}

    for set_index in range(num_sets_to_represent):
        primer_set = pcr_products[set_index]['products']
        if selected_primers is not None:
            combinations = itertools.combinations(selected_primers, 2)
        else:
            combinations = primer_set.keys()
        combinations = sorted(combinations)
        for combination in combinations:
            combination = tuple(combination)
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
