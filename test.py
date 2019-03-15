from primer_xp import (parse_fasta, KmerLocationGenerator)
import config
import unittest


class TestKmerGenerator(unittest.TestCase):

    def test_generate_kmers(self):
        genome_fasta_fhand = open(str(config.KMERS_GENERATOR_FASTA), 'rb')
        genome = parse_fasta(genome_fasta_fhand)
        kmer_locations = KmerLocationGenerator(config.KMER_LENGTH)
        kmers_locations_tuples = list(
            kmer_locations.generate_kmer_locations(genome))
        print(kmers_locations_tuples[51])

        assert kmers_locations_tuples[0] == (b'ATGGTAGG', (b'test_kmers', 0))
        assert kmers_locations_tuples[51] == (b'ATGTCGCT', (b'test_kmers', 75))


if __name__ == "__main__":
    unittest.main()
