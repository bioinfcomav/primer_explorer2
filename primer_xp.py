import io
from itertools import groupby
from collections import Counter
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

                if not dust_score_is_ok(kmer, windowsize=dust_windowsize,
                                        windowstep=dust_windowstep,
                                        threshold=dust_threshold):
                    continue

                if not gc_is_ok(kmer, min_gc=min_gc, max_gc=max_gc):
                    continue

                if blacklisted_seqs_in_seq(kmer):
                    continue

                location = (seq["seq_id"], location)
                kmer_counter[kmer] += 1
                yield kmer, location


def no_name(genome_fhand, kmer_len, filtering_criteria):

    genome = parse_fasta(genome_fhand)
    kmer_generator = KmerLocationGenerator(kmer_len)
    kmer_location_tuples = kmer_generator.generate_kmer_locations(genome)
    print(len(list(kmer_location_tuples)))
    print(kmer_generator.kmer_counter.most_common(10))
    print('hola')


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
    genome_fasta_fhand = open('chrom1.500k.fasta', 'rb')
    no_name(genome_fasta_fhand, kmer_len=8, filtering_criteria='')
