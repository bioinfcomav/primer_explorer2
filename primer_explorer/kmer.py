import hashlib
import pickle
from itertools import groupby, combinations
from collections import Counter, namedtuple, OrderedDict, defaultdict
from operator import itemgetter

from primer_explorer.dust_score import dust_score_is_ok
from primer_explorer.utils import get_fhand
from primer_explorer.seq_filters import gc_is_ok, blacklisted_seqs_in_seq
from primer_explorer.regions import GenomeRegion, GenomeRegions
from primer_explorer.primer3.primer3 import (reverse_complement,
                                             kmer_validated_by_primer3)

GT_BIN = b'>'
GT_CHAR = ord(GT_BIN)

KmerAndLocation = namedtuple('KmerAndLocation', ('seq', 'chrom_location',
                                                 'is_heterochromatic'))


def get_first_that_complies(iterator, condition):
    for item in iterator:
        if condition(item):
            return item
    raise ValueError('No item matches the given condition')


class KmerLocationGenerator:

    def __init__(self, seqs, kmer_len, heterochromatic_regions=None,
                 min_gc=0.35, max_gc=0.75, dust_windowsize=64,
                 dust_windowstep=32, dust_threshold=7, kmers_to_keep=None):
        self.kmer_len = kmer_len
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.dust_windowstep = dust_windowstep
        self.dust_windowsize = dust_windowsize
        self.dust_threshold = dust_threshold

        # as we only need it to make checks with in, the speed improvement is
        # using a dict instead of a list is huge
        if kmers_to_keep is not None:
            kmers_to_keep = {kmer: '' for kmer in kmers_to_keep}
        self.kmer_to_keep = kmers_to_keep
        self.kmer_counters = {True: Counter(), False: Counter()}

        if heterochromatic_regions is None:
            heterochromatic_regions = []
        self._heterochromatic_regions = iter(heterochromatic_regions)
        self._current_heterochromatic_region = None
        self._seqs = seqs

    def _region_is_heterochromatic(self, kmer_region):
        current_heterochromatic_region = self._current_heterochromatic_region
        if current_heterochromatic_region is None:

            def heterochromatic_region_is_next(region):
                return kmer_region < region or kmer_region.overlaps(region)

            try:
                current_heterochromatic_region = get_first_that_complies(self._heterochromatic_regions,
                                                                         heterochromatic_region_is_next)
                self._current_heterochromatic_region = current_heterochromatic_region
            except ValueError:
                current_heterochromatic_region = None
        if current_heterochromatic_region is None:
            return False

        if current_heterochromatic_region.overlaps(kmer_region):
            return True
        else:
            if current_heterochromatic_region < kmer_region:
                self._current_heterochromatic_region = None
            return False

    def generate_kmer_locations(self):
        kmer_len = self.kmer_len
        min_gc = self.min_gc
        max_gc = self.max_gc
        dust_windowstep = self.dust_windowstep
        dust_windowsize = self.dust_windowsize
        dust_threshold = self.dust_threshold
        kmer_counters = self.kmer_counters
        kmers_to_keep = self.kmer_to_keep
        seqs = self._seqs

        for seq in seqs:
            for location in range(len(seq['seq']) - kmer_len + 1):
                kmer = seq['seq'][location: location + kmer_len]
                if kmers_to_keep and kmer not in kmers_to_keep:
                    continue
                if not dust_score_is_ok(kmer, windowsize=dust_windowsize,
                                        windowstep=dust_windowstep,
                                        threshold=dust_threshold):
                    continue

                if not gc_is_ok(kmer, min_gc=min_gc, max_gc=max_gc):
                    continue

                if blacklisted_seqs_in_seq(kmer):
                    continue

                chrom = seq['seq_id']
                # location_start = GenomeLocation(chrom, location)
                location_start = (chrom, location)
                kmer_region = GenomeRegion(chrom,
                                           location,
                                           location + kmer_len)
                is_heterochromatic = self._region_is_heterochromatic(
                    kmer_region)

                kmer_counters[is_heterochromatic][kmer] += 1
                yield KmerAndLocation(kmer, location_start, is_heterochromatic)


def parse_fasta(fhand, ignore_softmask=True):
    fasta_chunks = (x[1] for x in groupby(fhand,
                                          lambda line: line[0] == GT_CHAR))
    while True:
        try:
            header = list(next(fasta_chunks))[0][1:]
            seq_id = header.split()[0]
        except StopIteration:
            break
        try:
            seq_lines = next(fasta_chunks)
        except StopIteration:
            raise RuntimeError('Sequence has no sequence: ' + header.decode())

        seq = b''.join(line.strip() for line in seq_lines)

        if ignore_softmask:
            yield {'seq_id': seq_id, 'seq': seq.upper()}
        else:
            yield {'seq_id': seq_id, 'seq': seq}


def get_revcomp_locs(kmer, kmer_generator, get_heterochromatic_locs=False):
    return kmer_generator.kmer_counters[False][reverse_complement(kmer)]


def get_top_kmers_by_euchromatin_abundance(kmer_generator, max_num_kmers):
    kmer_locs = Counter()
    for kmer, num_locs in kmer_generator.kmer_counters[False].items():
        kmer_locs[kmer] = num_locs + get_revcomp_locs(kmer, kmer_generator)
    sorted_by_abundance_kmers = [kmer[0] for kmer in kmer_locs.most_common()]
    return sorted_by_abundance_kmers[:max_num_kmers]


def get_top_kmers_by_minimum_abundance(kmer_generator, min_abundance=1000,
                                       max_num_kmers=1000):
    grouped_counters = defaultdict(int)
    for _bool in [True, False]:
        for kmer, counts in kmer_generator.kmer_counters[_bool].items():
            grouped_counters[kmer] += counts

    sorted_kmers = [kmer[0] for kmer in sorted(grouped_counters.items(),
                                               key=itemgetter(1),
                                               reverse=True)
                    if kmer[1] >= min_abundance][:max_num_kmers]
    return sorted_kmers


def calculate_euchromatin_ratios(kmer_generator):
    euchromatin_ratio = {}
    euchromatin_abundance = {kmer: abundance for kmer, abundance in kmer_generator.kmer_counters[False].items()}
    for kmer in euchromatin_abundance:
        euchromatin_abundance[kmer] = euchromatin_abundance[kmer] + euchromatin_abundance.get(reverse_complement(kmer), 0)

    heterochromatin_abundance = {kmer: abundance for kmer, abundance in kmer_generator.kmer_counters[True].items()}
    for kmer in heterochromatin_abundance:
        heterochromatin_abundance[kmer] = heterochromatin_abundance[kmer] + heterochromatin_abundance.get(reverse_complement(kmer), 0)
    for kmer in euchromatin_abundance:
        if kmer not in heterochromatin_abundance:
            euchromatin_ratio[kmer] = 1
        else:
            euchromatin_ratio[kmer] = float(euchromatin_abundance[kmer] / (heterochromatin_abundance[kmer] + euchromatin_abundance[kmer]))
    return euchromatin_ratio


def get_top_kmers_by_euchromatin_ratio(kmer_generator, max_num_kmers=1000, min_abundance=1000):
    ratios = calculate_euchromatin_ratios(kmer_generator)
    sorted_kmers = [kmer[0] for kmer in sorted(ratios.items(), key=itemgetter(1), reverse=True)][:max_num_kmers]
    return sorted_kmers


def filter_kmers_by_heterochromatin_stats(kmer_generator, max_num_kmers=1000,
                                          criterion="euchromatin abundance",
                                          min_abundance=1000):
    if criterion == "euchromatin abundance":
        return get_top_kmers_by_euchromatin_abundance(kmer_generator,
                                                      max_num_kmers)
    if criterion == "euchromatin ratio":
        return get_top_kmers_by_euchromatin_ratio(kmer_generator,
                                                  max_num_kmers)
    if criterion == "minimun total abundance":
        return get_top_kmers_by_minimum_abundance(kmer_generator,
                                                  max_num_kmers=max_num_kmers,
                                                  min_abundance=min_abundance)
    if criterion == "total abundance":
        return get_top_kmers_by_minimum_abundance(kmer_generator,
                                                  max_num_kmers=max_num_kmers,
                                                  min_abundance=0)
    raise RuntimeError(
        "No available criterion for filtering: {}".format(criterion))


def filter_kmers_by_primer3(kmers, num_kmers_to_keep):
    filtered_kmers = []
    for kmer in kmers:
        if not kmer_validated_by_primer3(kmer.decode()):
            continue
        else:
            filtered_kmers.append(kmer)
            if len(filtered_kmers) == num_kmers_to_keep:
                break
    return filtered_kmers


def filter_kmers_by_revcomp(kmers, kmer_counters, count_heterochromatic=False):
    filt_out_kmers = []
    for combination in combinations(kmers, 2):
        fwd_primer = combination[0]
        rev_primer = combination[1]
        if fwd_primer.decode() == reverse_complement(rev_primer).decode():
            num_products_fwd = kmer_counters[False][fwd_primer]
            num_products_rev = kmer_counters[False][rev_primer]
            if count_heterochromatic:
                num_products_fwd += kmer_counters[True][fwd_primer]
                num_products_rev += kmer_counters[True][rev_primer]
            if num_products_fwd < num_products_rev:
                filt_out_kmers.append(fwd_primer)
            else:
                filt_out_kmers.append(rev_primer)
    kmers = [kmer for kmer in kmers if kmer not in filt_out_kmers]
    return kmers


def pack_kmer_locations(kmer_locations):
    packed_kmers = OrderedDict()
    for kmer in kmer_locations:
        if kmer[0] in packed_kmers:
            packed_kmers[kmer[0]].append(kmer)
        else:
            packed_kmers[kmer[0]] = [kmer]
    return packed_kmers


def generate_kmer_locations(genome_fhand, kmer_len, heterochromatic_regions,
                            num_kmers_to_keep=1000, kmers_to_keep=None):
    genome = parse_fasta(genome_fhand)
    kmer_generator = KmerLocationGenerator(genome, kmer_len, heterochromatic_regions,
                                           kmers_to_keep=kmers_to_keep)
    kmer_locations = kmer_generator.generate_kmer_locations()
    kmer_locations = list(kmer_locations)
    filt_kmers_by_het_stats = filter_kmers_by_heterochromatin_stats(kmer_generator,
                                                                    criterion="euchromatin abundance",
                                                                    max_num_kmers=num_kmers_to_keep)
    filt_kmers_by_primer3 = filter_kmers_by_primer3(filt_kmers_by_het_stats, num_kmers_to_keep)
    filtered_kmers = filter_kmers_by_revcomp(filt_kmers_by_primer3, kmer_generator.kmer_counters)
    packed_kmers = pack_kmer_locations(kmer_locations)
    return filtered_kmers, packed_kmers


def get_kmers(genome_fpath, heterochromatic_regions_fpath, kmer_len, cache_dir,
              num_kmers_to_keep=1000, kmers_to_keep=None):
    key = str(genome_fpath)
    key += str(heterochromatic_regions_fpath)
    key += str(kmer_len)
    if kmers_to_keep:
        key += str(kmers_to_keep)
    key = hashlib.md5(key.encode()).hexdigest()
    cache_fpath = cache_dir / ('kmers_' + key)

    if cache_fpath.exists():
        kmers, kmers_locations = pickle.load(cache_fpath.open('rb'))
        return kmers, kmers_locations
    else:
        if kmers_to_keep:
            most_freq_kmers = kmers_to_keep
        else:
            genome_fhand = get_fhand(genome_fpath)
            regions_fhand = get_fhand(heterochromatic_regions_fpath)
            counters = count_kmers(genome_fhand, kmer_len,
                                   regions_fhand=regions_fhand)
            genome_fhand.close()
            regions_fhand.close()
            most_freq_kmers = [k[0] for k in counters[False].most_common(num_kmers_to_keep)]

        regions_fhand = get_fhand(heterochromatic_regions_fpath)
        genome_fhand = get_fhand(genome_fpath)

        heterochromatic_regions = GenomeRegions(bed_fhand=regions_fhand)
        kmers, kmers_locations = generate_kmer_locations(genome_fhand,
                                                         kmer_len=kmer_len,
                                                         heterochromatic_regions=heterochromatic_regions,
                                                         kmers_to_keep=most_freq_kmers)
        cache_fhand = open(str(cache_fpath), "wb")
        pickle.dump((kmers, kmers_locations), cache_fhand, pickle.HIGHEST_PROTOCOL)
        cache_fhand.close()
        regions_fhand.close()
        genome_fhand.close()

        return kmers, kmers_locations


def count_kmers(fhand, kmer_len, regions_fhand=None):

    seqs = parse_fasta(fhand)
    if regions_fhand:
        heterochromatic_regions = GenomeRegions(bed_fhand=regions_fhand)
    else:
        heterochromatic_regions = None

    kmer_generator = KmerLocationGenerator(seqs, kmer_len, heterochromatic_regions=heterochromatic_regions)

    for _ in kmer_generator.generate_kmer_locations():
        _
    return kmer_generator.kmer_counters
