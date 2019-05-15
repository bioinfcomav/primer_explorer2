import hashlib
import operator
import pickle
import sys
import config
from itertools import groupby, combinations
from collections import Counter, namedtuple, OrderedDict, defaultdict

from combinations import get_compatible_groups_of_primers
from dust_score import dust_score_is_ok
from primer3.primer3 import kmer_validated_by_primer3, reverse_complement
from seq_filters import gc_is_ok, blacklisted_seqs_in_seq

GT = '>'
GT_BIN = b'>'
GT_CHAR = ord(GT_BIN)
GenomeLocation = namedtuple('GenomeLocation', ('chrom', 'location'))
KmerAndLocation = namedtuple('KmerAndLocation', ('seq', 'chrom_location', 'is_heterochromatic'))
PrimerAndLocation = namedtuple('PrimerAndLocation', ('strand', 'seq', 'chrom_location', 'is_heterochromatic'))


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


def location_is_in_region(location, region):
    if location.chrom != region.chrom:
        return False
    if region.start <= location.location < region.stop:
        return True
    else:
        return False


def regions_overlap(region1, region2):
    if region1.chrom != region2.chrom:
        return False


class KmerLocationGenerator:

    def __init__(self, seqs, kmer_len, heterochromatic_regions=None,
                 min_gc=0.35, max_gc=0.75,
                 dust_windowsize=64, dust_windowstep=32, dust_threshold=7):
        self.kmer_len = kmer_len
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.dust_windowstep = dust_windowstep
        self.dust_windowsize = dust_windowsize
        self.dust_threshold = dust_threshold
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
                return kmer_region.start < region.start or kmer_region.overlaps(region)

            try:
                current_heterochromatic_region = get_first_that_complies(self._heterochromatic_regions, heterochromatic_region_is_next)
                self._current_heterochromatic_region = current_heterochromatic_region
            except ValueError:
                current_heterochromatic_region = None
        if current_heterochromatic_region is None:
            return False

        if current_heterochromatic_region.overlaps(kmer_region):
            return True
        else:
            if current_heterochromatic_region.start < kmer_region.start or current_heterochromatic_region.chrom != kmer_region.chrom:
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
        seqs = self._seqs

        for seq_idx, seq in enumerate(seqs):
            for location in range(len(seq['seq']) - kmer_len + 1):
                kmer = seq['seq'][location: location + kmer_len]

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


class GenomeRegions:

    def __init__(self, bed_fhand):
        self._bed_fhand = bed_fhand

    def __iter__(self):
        return self

    def __next__(self):
        line_items = next(self._bed_fhand).split()
        chrom = line_items[0]
        start = int(line_items[1])
        stop = int(line_items[2])
        return GenomeRegion(chrom, start, stop)


class GenomeRegion:

    def __init__(self, chrom, start, stop):
        self.chrom = chrom
        self.start = start
        self.stop = stop

    def __lt__(self, region2):
        if self.chrom != region2.chrom:
            return True
        return self.start < region2.start

    def overlaps(self, region2):
        if self.chrom != region2.chrom:
            return False
        if self.stop <= region2.start:
            return False
        if self.start >= region2.stop:
            return False
        return True


def test_genome_region():
    region1 = GenomeRegion('chrom1', 1, 2)
    region2 = GenomeRegion('chrom2', 1, 2)
    assert region1 < region2

    assert region1 < GenomeRegion('chrom1', 2, 3)
    assert region1 < GenomeRegion('chrom1', 2, 3)

    assert not GenomeRegion('c', 0, 2).overlaps(GenomeRegion('c', 2, 4))
    assert not GenomeRegion('c', 2, 4).overlaps(GenomeRegion('c', 0, 2))
    assert not GenomeRegion('c', 0, 2).overlaps(GenomeRegion('c', 3, 4))
    assert not GenomeRegion('c', 3, 4).overlaps(GenomeRegion('c', 0, 2))
    assert GenomeRegion('c', 0, 2).overlaps(GenomeRegion('c', 0, 1))
    assert GenomeRegion('c', 10, 12).overlaps(GenomeRegion('c', 9, 11))
    assert GenomeRegion('c', 10, 12).overlaps(GenomeRegion('c', 11, 13))
    assert not GenomeRegion('c', 10, 12).overlaps(GenomeRegion('c', 12, 13))
    assert GenomeRegion('c', 10, 15).overlaps(GenomeRegion('c', 11, 13))


def get_first_that_complies(iterator, condition):
    for item in iterator:
        if condition(item):
            return item
    raise ValueError('No item matches the given condition')


def test_get_first_that_complies():

    def is_odd(x):
        return x % 2

    assert get_first_that_complies(iter([1, 3, 5]), is_odd) == 1
    assert get_first_that_complies(iter([2, 3, 4]), is_odd) == 3
    try:
        assert get_first_that_complies(iter([2, 4]), is_odd)
    except ValueError:
        pass


def pack_kmer_locations(kmer_locations):
    packed_kmers = OrderedDict()
    for kmer in kmer_locations:
        if kmer[0] in packed_kmers:
            packed_kmers[kmer[0]].append(kmer)
        else:
            packed_kmers[kmer[0]] = [kmer]
    return packed_kmers


def get_selected_kmer_locations(kmer_locations, selected_kmers):
    selected_kmers_locations = []
    for selected_kmer in selected_kmers:
        for kmer_location in kmer_locations:
            if kmer_location[0] == selected_kmer:
                selected_kmers_locations.append(kmer_location)
    return selected_kmers_locations


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


def get_revcomp_locs(kmer, kmer_generator, get_heterochromatic_locs=False):
    return kmer_generator.kmer_counters[False][reverse_complement(kmer)]


def get_top_kmers_by_euchromatin_abundance(kmer_generator, max_num_kmers):
    kmer_locs = Counter()
    for kmer, num_locs in kmer_generator.kmer_counters[False].items():
        kmer_locs[kmer] = num_locs + get_revcomp_locs(kmer, kmer_generator)
    sorted_by_abundance_kmers = [kmer[0] for kmer in kmer_locs.most_common()][:max_num_kmers]
    return sorted_by_abundance_kmers


def get_top_kmers_by_minimum_abundance(kmer_generator, min_abundance=1000, max_num_kmers=1000):
    grouped_counters = defaultdict(int)
    for _bool in [True, False]:
        for kmer, counts in kmer_generator.kmer_counters[_bool].items():
            grouped_counters[kmer] += counts

    sorted_kmers = [kmer[0] for kmer in sorted(grouped_counters.items(), key=operator.itemgetter(1), reverse=True) if kmer[1] >= min_abundance][:max_num_kmers]
    return sorted_kmers


def get_top_kmers_by_euchromatin_ratio(kmer_generator, max_num_kmers=1000, min_abundance=1000):
    ratios = calculate_euchromatin_ratios(kmer_generator)
    sorted_kmers = [kmer[0] for kmer in sorted(ratios.items(), key=operator.itemgetter(1), reverse=True)][:max_num_kmers]
    return sorted_kmers


def filter_kmers_by_heterochromatin_stats(kmer_generator, criterion="euchromatin abundance", max_num_kmers=1000, min_abundance=1000):
    if criterion == "euchromatin abundance":
        return get_top_kmers_by_euchromatin_abundance(kmer_generator, max_num_kmers)
    if criterion == "euchromatin ratio":
        return get_top_kmers_by_euchromatin_ratio(kmer_generator, max_num_kmers)
    if criterion == "minimun total abundance":
        return get_top_kmers_by_minimum_abundance(kmer_generator, max_num_kmers=max_num_kmers, min_abundance=min_abundance)
    if criterion == "total abundance":
        return get_top_kmers_by_minimum_abundance(kmer_locations, kmer_generator, max_num_kmers=max_num_kmers,
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


def generate_kmer_locations(genome_fhand, kmer_len, heterochromatic_regions, num_kmers_to_keep=1000):
    genome = parse_fasta(genome_fhand)
    kmer_generator = KmerLocationGenerator(
        genome, kmer_len, heterochromatic_regions)
    kmer_locations = kmer_generator.generate_kmer_locations()
    kmer_locations = list(kmer_locations)
    filt_kmers_by_het_stats = filter_kmers_by_heterochromatin_stats(kmer_generator,
                                                                    criterion="euchromatin abundance",
                                                                    max_num_kmers=num_kmers_to_keep)
    filt_kmers_by_primer3 = filter_kmers_by_primer3(filt_kmers_by_het_stats, num_kmers_to_keep)
    filtered_kmers = filter_kmers_by_revcomp(filt_kmers_by_primer3, kmer_generator.kmer_counters)
    packed_kmers = pack_kmer_locations(kmer_locations)
    return filtered_kmers, packed_kmers


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


def annotate_length_viable_products(pcr_products, min_viable_length=config.MIN_PCR_VIABLE_LENGTH):
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
        if all(primer.is_heterochromatic == True for primer in product):
            euchromatin_locations["heterochromatic"] += 1
        elif all(primer.is_heterochromatic == False for primer in product):
            euchromatin_locations["euchromatic"] += 1
        else:
            euchromatin_locations["mixed"] += 1
    return euchromatin_locations


def count_total_products(pcr_products):
    return len(pcr_products)


def annotate_products(pcr_products):
    annotations = defaultdict(dict)
    for primers, products in pcr_products.items():
        annotations[primers]["total_products"] = count_total_products(products)
        annotations[primers]["pcr_revcomp_dimers"] = annotate_reverse_complement_dimers(products)
        annotations[primers]["length_viable_products"] = annotate_length_viable_products(products)
        annotations[primers]["euchromatin_products"] = annotate_products_by_euchromatin_region(products)
    return annotations


def get_kmers(genome_fpath, heterochromatic_regions_fpath, kmer_len, cache_dir):
    key = str(genome_fpath)
    key += str(heterochromatic_regions_fpath)
    key += str(kmer_len)
    key = hashlib.md5(key.encode()).hexdigest()
    cache_fpath = cache_dir / ('kmers_' + key)
    if cache_fpath.exists():
        kmers, kmers_locations = pickle.load(cache_fpath.open('rb'))
        return kmers, kmers_locations
    else:
        cache_fhand = open(str(cache_fpath), "wb")
        heterochromatic_regions = GenomeRegions(open(heterochromatic_regions_fpath, 'rb'))
        genome_fasta_fhand = open(genome_fpath, 'rb')
        kmers, kmers_locations = generate_kmer_locations(genome_fasta_fhand, kmer_len=kmer_len, heterochromatic_regions=heterochromatic_regions)
        pickle.dump((kmers, kmers_locations), cache_fhand, pickle.HIGHEST_PROTOCOL)
        return kmers, kmers_locations


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
