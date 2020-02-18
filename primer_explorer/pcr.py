from collections import namedtuple
from itertools import combinations
import pyranges as pr
import pandas as pd

from primer_explorer.primer3.primer3 import reverse_complement
from primer_explorer.combinations import get_compatible_groups_of_primers
from primer_explorer.stats import filter_by_length

PrimerAndLocation = namedtuple('PrimerAndLocation', ('strand', 'seq', 'chrom_location', 'is_heterochromatic'))


def select_primers_combinations(kmers, num_compatible_groups=10):
    compatibility_primers_groups = get_compatible_groups_of_primers(kmers, num_compatible_groups=num_compatible_groups)
    return compatibility_primers_groups


def _get_pcr_products_old(fwd_primers_locations, rev_primers_locations, max_pcr_product_length):

    fwd_locs = [PrimerAndLocation(0, loc.seq, loc.chrom_location, loc.is_heterochromatic) for loc in fwd_primers_locations]
    locs = fwd_locs
    rev_locs = [PrimerAndLocation(1, loc.seq, loc.chrom_location, loc.is_heterochromatic) for loc in rev_primers_locations]
    locs.extend(rev_locs)
    locs = sorted(locs, key=lambda PrimerAndLocation: (PrimerAndLocation.chrom_location))
    # pcr_products = [(loc1, loc2) for loc1, loc2 in zip(locs[:-1], locs[1:]) if (loc1.strand == 0) and (loc1.strand != loc2.strand) and (loc1.chrom_location[0] == loc2.chrom_location[0]) and abs(loc1.chrom_location[1] - loc2.chrom_location[1]) < max_pcr_product_length]
    pcr_products = []

    for loc1, loc2 in zip(locs[:-1], locs[1:]):
        if ((loc1.strand == 0) and (loc1.strand != loc2.strand) and
            (loc1.chrom_location[0] == loc2.chrom_location[0]) and
            (loc1.seq != reverse_complement(loc2.seq)) and
                (abs(loc1.chrom_location[1] - loc2.chrom_location[1]) < max_pcr_product_length)):

            pcr_products.append((loc1, loc2))
    return pcr_products


def kmer_locations_to_ranges(kmer_locations, reverse=False, extend_range=1):
    chroms = []
    starts = []
    ends = []
    strand = '-' if reverse else '+'
    primer_strand = 1 if reverse else 0
    seqs = []
    is_hetero = []
    indexed_primers = {}
    for loc in kmer_locations:
        chroms.append(loc.chrom_location[0].decode())

        if reverse:
            ends.append(loc.chrom_location[1])
            start = loc.chrom_location[1] - extend_range
            if start < 0:
                start = 0
            starts.append(start)
        else:
            starts.append(loc.chrom_location[1])
            ends.append(loc.chrom_location[1] + extend_range)

        seqs.append(reverse_complement(loc.seq) if reverse else loc.seq)
        is_hetero = loc.is_heterochromatic
        indexed_primers[(loc.chrom_location)] = PrimerAndLocation(primer_strand, loc.seq, loc.chrom_location, loc.is_heterochromatic)

    df = pd.DataFrame({'Chromosome': chroms, 'Start': starts, 'End': ends,
                       'seq': seqs, 'Strand': strand, 'is_het': is_hetero})

    return pr.PyRanges(df), indexed_primers


def _get_pcr_products(fwd_primers_locations, rev_primers_locations, max_pcr_product_length):
    fwd_ranges, fwr_primers = kmer_locations_to_ranges(fwd_primers_locations, reverse=False,
                                                       extend_range=max_pcr_product_length)
    rev_ranges, rev_primers = kmer_locations_to_ranges(rev_primers_locations, reverse=True,
                                                       extend_range=max_pcr_product_length)
#     print(fwd_ranges)
    pcr_products = []
    intersection = fwd_ranges.intersect(rev_ranges, strandedness="opposite")
#     print(intersection)
    for chrom in intersection.chromosomes:
        df = intersection[chrom].df

        for index in range(df.shape[0]):
            line = df.iloc[index]
            try:
                frw_primer = fwr_primers[(line.Chromosome.encode(), line.Start)]
                rev_primer = rev_primers[(line.Chromosome.encode(), line.End)]
            except KeyError:
                continue
            if frw_primer.strand != rev_primer.strand and frw_primer.seq != reverse_complement(rev_primer.seq):
                pcr_products.append((frw_primer, rev_primer))
    return pcr_products


def get_pcr_products(kmer_locations, primers, max_pcr_product_length=5000):
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


def pcr_products_to_regions(pcr_products, min_length, max_length, read_length, kmer_length):
    pcr_products = filter_by_length(pcr_products, min_length=min_length, max_length=max_length)
    for pcr_product in pcr_products:
        chrom_fwd = pcr_product[0].chrom_location
        chrom_rev = pcr_product[1].chrom_location
        if chrom_rev[1] - chrom_fwd[1] < read_length:
            yield "{}\t{}\t{}".format(chrom_fwd[0].decode(), chrom_fwd[1], chrom_rev[1] + kmer_length)
        else:
            yield "{}\t{}\t{}".format(chrom_fwd[0].decode(), chrom_fwd[1] + kmer_length, chrom_fwd[1] + read_length)
            yield "{}\t{}\t{}".format(chrom_rev[0].decode(), chrom_rev[1] - read_length, chrom_rev[1])
