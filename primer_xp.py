import io
from itertools import groupby
from collections import Counter, namedtuple

import gzip
from dust_score import dust_score_is_ok
from seq_filters import gc_is_ok, blacklisted_seqs_in_seq

GT = '>'
GT_BIN = b'>'
GT_CHAR = ord(GT_BIN)


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


GenomeLocation = namedtuple('GenomeLocation', ('chrom', 'location'))
KmerAndLocation = namedtuple(
    'KmerAndLocation', ('seq', 'chrom_location', 'is_heterochromatic'))


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
        # if current_heterochromatic_region is not None:
        #   print(current_heterochromatic_region.start)
        if current_heterochromatic_region is None:
            def heterochromatic_region_is_next(region):
                return kmer_region < region or kmer_region.overlaps(region)

            try:
                current_heterochromatic_region = get_first_that_complies(
                    self._heterochromatic_regions, heterochromatic_region_is_next)
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
                #location_start = GenomeLocation(chrom, location)
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
        if self.chrom < region2.chrom:
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


def generate_kmer_locations(genome_fhand, kmer_len, heterochromatic_regions):

    genome = parse_fasta(genome_fhand)
    kmer_generator = KmerLocationGenerator(
        genome, kmer_len, heterochromatic_regions)
    kmer_locations = kmer_generator.generate_kmer_locations()
    print(len(list(kmer_locations)))
    print(kmer_generator.kmer_counters[True].most_common(10))
    print(kmer_generator.kmer_counters[False].most_common(10))
    print('hola')


GENOME_FASTA = b'''>chrom1
ATGGTAGGGATAGATAGAATAGATAGATAGATTAGGGCTTCATCATACT
GATGGACGCTAGACTGATGCTAGCTGATGTCGCTAGCATATACGCAGAT
>chrom2
TATGGAGACTACGACCTAGCTACGGCATGATACGGCATACGCATACGCA
'''


TOMATO_CHROM_FASTA_GZ = '/home/jope/devel3/primer_explorer_old_version/genome/S_lycopersicum_chromosomes.3.00.chrom1.fasta.gz'
HETEROCHROMATIN_BED = '/home/jope/devel3/primer_experience/test_heterochromatin.bed'


if __name__ == '__main__':
    test_genome_region()
    test_get_first_that_complies()
    heterochromatic_regions = GenomeRegions(open(HETEROCHROMATIN_BED, 'rb'))

    genome_fasta_fhand = io.BytesIO(GENOME_FASTA)
    # genome_fasta_fhand = io.StringIO(GENOME_FASTA.decode())
    genome_fasta_fhand = gzip.open(TOMATO_CHROM_FASTA_GZ, 'rb')
    #genome_fasta_fhand = open('chrom1.500k.fasta', 'rb')
    #genome_fasta_fhand = open('chrom1.100k.fasta', 'rb')
    generate_kmer_locations(genome_fasta_fhand, kmer_len=8,
                            heterochromatic_regions=heterochromatic_regions)
