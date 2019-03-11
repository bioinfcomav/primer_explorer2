from primer_xp import (parse_fasta, KmerLocationGenerator)
import config
import unittest


class TestKmerGenerator(unittest.TestCase):

    def test_generate_kmers(self):
        genome_fasta_fhand = open(str(config.KMERS_GENERATOR_FASTA), 'rb')
        genome = parse_fasta(genome_fasta_fhand)
        kmer_locations = KmerLocationGenerator(config.KMER_LENGTH)
        kmer_location_tuples = list(kmer_locations.generate_kmer_locations(genome))

        assert kmer_location_tuples[0] == (b'ATGGTAGG', (b'test_kmers', 0))
        assert kmer_location_tuples[51] == (b'TGGACGCT', (b'test_kmers', 51))


if __name__ == "__main__":
    unittest.main()
