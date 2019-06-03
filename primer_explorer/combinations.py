from itertools import zip_longest
from primer_explorer.primer3.primer3 import (combinations_validated_by_primer3,
                                             reverse_complement)


def get_compatible_groups_of_primers(primers, seed_index=0, num_compatible_groups=10):
    compatibility_primers_groups = []
    seed_primer = _get_first_cadidate_primer_from_index_seed(
        seed_index, primers)
    compatible_primers_group, incompatible_primers = make_compatibility_group(
        seed_primer, primers)
    compatibility_primers_groups.append(compatible_primers_group)
    if num_compatible_groups == 1 or len(incompatible_primers) == 0:
        return compatibility_primers_groups
    else:
        for primer in incompatible_primers:
            seed_primer = [primer]
            compatible_primers_group, _ = make_compatibility_group(
                seed_primer, primers)
            compatibility_primers_groups.append(compatible_primers_group)
            if len(compatibility_primers_groups) >= num_compatible_groups:
                break
    return compatibility_primers_groups


def make_compatibility_group(seed_primer, primers, group_len=10):
    compatible_primers = seed_primer[:]
    incompatible_primers = []
    for primer in primers:
        if primer in compatible_primers:
            continue
        else:
            if combinations_validated_by_primer3(compatible_primers, primer.decode()):
                compatible_primers.append(primer)
                if len(compatible_primers) == group_len:
                    break
            else:
                incompatible_primers.append(primer)
    return compatible_primers, incompatible_primers


def divide_kmers_into_groups(kmers):
    group1 = [kmers[0]]
    group2 = []
    for kmers in grouper(4, kmers[1:]):
        group2 = _unpack_kmers(group2, kmers[:2])
        group1 = _unpack_kmers(group1, kmers[2:])
    return group1, group2


def _unpack_kmers(group, kmers):
    kmers = [kmer for kmer in kmers if kmer is not None]
    group += kmers
    return group


def _get_first_cadidate_primer_from_index_seed(idx, sorted_kmers):
    return [sorted_kmers[idx]]


def _primer_is_compatibility_group_seed(primer, compatibility_group_primers):
    return primer == compatibility_group_primers[0]


def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def filter_out_kmers_by_revcomp(combination, sorted_primer_locations):
    fwd_primer = combination[0]
    rev_primer = combination[1]
    if fwd_primer == reverse_complement(rev_primer):
        num_products_fwd = len(sorted_primer_locations[fwd_primer])
        num_products_rev = len(sorted_primer_locations[rev_primer])
        if num_products_fwd < num_products_rev:
            return fwd_primer
        else:
            return rev_primer
