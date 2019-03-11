import io
import gzip
from itertools import (groupby, zip_longest, combinations, islice)
import itertools
import collections
from collections import Counter

#from parse_fasta import parse_fasta
#from generate_kmer_locations import generate_kmer_locations

GT = '>'
GT_BIN = b'>'
GT_CHAR = ord(GT_BIN)
G_CHAR = ord(b'G')
C_CHAR = ord(b'C')


def parse_fasta(fhand, ignore_softmask=True):
    fasta_chunks = (x[1] for x in groupby(fhand, lambda line: line[0] == GT_CHAR))
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


def n_in_seq(seq):
    return b'N' in seq


def calculate_gc(seq):
    count = collections.Counter(seq)
    g_count = count[G_CHAR]
    c_count = count[C_CHAR]
    gc = (g_count + c_count) / len(seq)
    return gc


_GC_RESULT_CACHE = {}


def gc_is_ok(kmer, min_gc, max_gc):
    try:
        return _GC_RESULT_CACHE[kmer]
    except KeyError:
        if min_gc <= calculate_gc(kmer) <= max_gc:
            result = True
        else:
            result = False
        _GC_RESULT_CACHE[kmer] = result
        #print(kmer, calculate_gc(kmer), result)
        return result


def _rolling_window_serie(serie, window, length_, step):
    '''It yields lists of items with a window number of elements'''
    return (serie[i:i + window] for i in range(0, length_ - window + 1, step))


def _rolling_window_iter(iterator, window, step):
    '''It yields lists of items with a window number of elements giving
     an iterator'''
    items = []
    for item in iterator:
        if len(items) >= window:
            yield items
            items = items[step:]
        items.append(item)
    else:
        if len(items) >= window:
            yield items


def rolling_window(iterator, window, step=1):
    'It yields lists of items with a window number of elements'
    try:
        length_ = len(iterator)
    except TypeError:
        length_ = None
    if length_ is None:
        return _rolling_window_iter(iterator, window, step)
    else:
        return _rolling_window_serie(iterator, window, length_, step)


def _calculate_rawscore(string):
    'It returns a non-normalized dustscore'
    triplet_counts = Counter()
    for triplet in rolling_window(string, 3):
        # It should do something with non ATCG, but we sacrifice purity for
        # speed. Maybe we should reconsider this
        triplet_counts[triplet.upper()] += 1

    return sum(tc * (tc - 1) * 0.5 for tc in triplet_counts.values())


def calculate_dust_score(seq, windowsize, windowstep):
    '''It returns the dust score.
    From: "A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA
    Sequences"
    doi:10.1089/cmb.2006.13.1028
    and re-implemented from PRINSEQ
    '''
    length = len(seq)
    if length == 3:
        return 0
    if length <= 5:
        return None

    dustscores = []
    if length > windowsize:
        windows = 0
        for seq_in_win in rolling_window(seq, windowsize, windowstep):
            score = _calculate_rawscore(seq_in_win)
            dustscores.append(score / (windowsize - 2))
            windows += 1
        remaining_seq = seq[windows * windowstep:]
    else:
        remaining_seq = seq

    if len(remaining_seq) > 5:
        length = len(remaining_seq)
        score = _calculate_rawscore(remaining_seq)
        dustscore = score / (length - 3) * (windowsize - 2) / (length - 2)
        dustscores.append(dustscore)

    # max score should be 100 not 31
    dustscore = sum(dustscores) / len(dustscores) * 100 / 31
    return dustscore


_DUST_CACHE = {}


def dust_score_is_ok(seq, windowsize, windowstep, threshold):
    try:
        return _DUST_CACHE[seq]
    except KeyError:
        dust = calculate_dust_score(seq, windowsize, windowstep)
        result = dust < threshold
        _DUST_CACHE[seq] = result
        return result


class KmerLocationGenerator():
    def __init__(self, kmer_len, min_gc=0.35, max_gc=0.75,
                 dust_windowsize=64, dust_windowstep=32, dust_threshold=7):
        self.kmer_len = kmer_len
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.dust_windowstep = dust_windowstep
        self.dust_windowsize = dust_windowsize
        self.dust_threshold = dust_threshold
        self.kmer_counter = Counter()

    def generate_kmer_locations(self, seqs):
        kmer_len = self.kmer_len
        min_gc = self.min_gc
        max_gc = self.max_gc
        dust_windowstep = self.dust_windowstep
        dust_windowsize = self.dust_windowsize
        dust_threshold = self.dust_threshold
        kmer_counter = self.kmer_counter

        for seq_idx, seq in enumerate(seqs):
            for location in range(len(seq['seq']) - kmer_len + 1):
                kmer = seq['seq'][location: location + kmer_len]

                if not dust_score_is_ok(kmer, windowsize=dust_windowsize, windowstep=dust_windowstep, threshold=dust_threshold):
                    continue

                if not gc_is_ok(kmer, min_gc=min_gc, max_gc=max_gc):
                    continue

                if n_in_seq(kmer):
                    continue

                location = (seq["seq_id"], location)
                kmer_counter[kmer] += 1
                yield kmer, location


def filter_kmers_by_seq_characteristics(kmers_to_check, filtering_criteria):
    results = {"total": len(kmers_to_check)}

    if 'gc' in filtering_criteria:
        filter_funct = partial(kmer_filter_by_gc,
                               min_gc=filtering_criteria['gc']['min_gc'],
                               max_gc=filtering_criteria['gc']['max_gc'])
        filtered_kmers = list(filter(filter_funct, kmers_to_check))
        results["filtered_by_gc"] = len(kmers_to_check) - len(filtered_kmers)
        total_num_kmers = len(filtered_kmers)

    if 'seqs_blacklist' in filtering_criteria:
        filter_funct = partial(kmer_filter_by_seq,
                               seqs=filtering_criteria["seqs_blacklist"])
        filtered_kmers = list(filter(filter_funct, filtered_kmers))
        results["filtered_by_blacklisted_seq"] = total_num_kmers - len(filtered_kmers)
        total_num_kmers = len(filtered_kmers)

    if "default_dust_threshold" in filtering_criteria:
        filter_funct = partial(kmer_filter_by_dust_score,
                               dustscore=filtering_criteria["default_dust_threshold"])
        filtered_kmers = list(filter(filter_funct, filtered_kmers))
        results["filtered_by_dust_score"] = total_num_kmers - len(filtered_kmers)
        total_num_kmers = len(filtered_kmers)

    return filtered_kmers, results


def no_name(genome_fhand, kmer_len, filtering_criteria):

    genome = parse_fasta(genome_fhand)
    kmer_generator = KmerLocationGenerator(kmer_len)
    kmer_location_tuples = kmer_generator.generate_kmer_locations(genome)
    print(len(list(kmer_location_tuples)))
    print(kmer_generator.kmer_counter.most_common(10))
    print('hola')
    print(calculate_dust_score(b'CCCAAAAT', windowsize=64, windowstep=32))

GENOME_FASTA = b'''>chrom1
ATGGTAGGGATAGATAGAATAGATAGATAGATTAGGGCTTCATCATACT
GATGGACGCTAGACTGATGCTAGCTGATGTCGCTAGCATATACGCAGAT
>chrom2
TATGGAGACTACGACCTAGCTACGGCATGATACGGCATACGCATACGCA
'''


TOMATO_CHROM_FASTA_GZ = '/home/jope/devel3/primer_explorer_old_version/genome/S_lycopersicum_chromosomes.3.00.chrom1.fasta.gz'

if __name__ == '__main__':
    genome_fasta_fhand = io.BytesIO(GENOME_FASTA)
    #genome_fasta_fhand = io.StringIO(GENOME_FASTA.decode())
    #genome_fasta_fhand = gzip.open(TOMATO_CHROM_FASTA_GZ, 'rb')
    genome_fasta_fhand = open('chrom1.500k.fasta', 'rb')
    genome_fasta_fhand = open('chrom1.100k.fasta', 'rb')
    no_name(genome_fasta_fhand, kmer_len=8, filtering_criteria='')
