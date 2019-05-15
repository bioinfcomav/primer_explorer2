import unittest
from unittest.mock import Mock
from pathlib import Path

from primer_explorer.regions import GenomeRegions
from primer_explorer.kmer import (KmerLocationGenerator, parse_fasta,
                                  filter_kmers_by_heterochromatin_stats,
                                  get_top_kmers_by_minimum_abundance,
                                  KmerAndLocation, filter_kmers_by_revcomp,
                                  get_first_that_complies)

from primer_explorer.primer3.primer3 import kmer_validated_by_primer3

TEST_DATA_PATHDIR = Path(__file__).parent.parent.joinpath('tests').joinpath('data')


class TestKmerGenerator(unittest.TestCase):

    def test_generate_kmers(self):
        with open(TEST_DATA_PATHDIR / 'heterochromatin_prueba.bed', 'rb') as regions_fhand, \
                open(TEST_DATA_PATHDIR / 'kmers.fa', 'rb') as genome_fhand:
            heterochromatic_regions = GenomeRegions(regions_fhand)
            genome = parse_fasta(genome_fhand)
            kmer_generator = KmerLocationGenerator(genome, 8, heterochromatic_regions)
            kmer_locations = kmer_generator.generate_kmer_locations()
            kmer_locations = list(kmer_locations)
            assert len(kmer_generator.kmer_counters[True].most_common(60)) == 18

    def test_get_first_that_complies(self):

        def is_odd(x):
            return x % 2

        assert get_first_that_complies(iter([1, 3, 5]), is_odd) == 1
        assert get_first_that_complies(iter([2, 3, 4]), is_odd) == 3
        try:
            assert get_first_that_complies(iter([2, 4]), is_odd)
        except ValueError:
            pass


class TestPrimer3Filter(unittest.TestCase):

    def test_primer3_kmers(self):
        kmers = [KmerAndLocation(b'GCAGACTTGTTCCTCCCCC', (b'chrom1', 1), True),
                 KmerAndLocation(b'GGGGGAGGAACAAGTCTGC', (b'chrom1', 2), False),
                 KmerAndLocation(b'TCCCCTTTAAAAAGGCTCAA', (b'chrom1', 3), True),
                 KmerAndLocation(b'TTGAGCCTTTTTAAAGGGGA', (b'chrom1', 4), False)]
        primer3_results = [kmer_validated_by_primer3(kmer[0].decode()) for kmer in kmers]
        assert primer3_results == [True, True, False, False]


class TestFilterKmersByReverseComplement(unittest.TestCase):

    def test_filter_by_revcomp(self):
        kmer_generator = Mock()
        kmers = [b'ATGC', b'GCAT', b'CCCC', b'GGGG']
        attrs = {'kmer_counters.return_value': {True: {b'ATGC': 0, b'GCAT': 0, b'CCCC': 10, b'GGGG': 0},
                                                False: {b'ATGC': 10, b'GCAT': 12, b'CCCC': 100, b'GGGG': 101}}}
        kmer_generator.configure_mock(**attrs)
        filtered_kmers = filter_kmers_by_revcomp(
            kmers, kmer_generator.kmer_counters(), count_heterochromatic=True)
        assert filtered_kmers == [b'GCAT', b'CCCC']


class TestSortByEuchromatinStats(unittest.TestCase):

    def test_sort_by_euchromatin_abundance(self):
        with open(TEST_DATA_PATHDIR / 'heterochromatin_abundance.bed', 'rb') as regions_fhand, \
                open(TEST_DATA_PATHDIR / "kmers_abundance.fa", 'rb') as genome_fhand:
            heterochromatic_regions = GenomeRegions(regions_fhand)
            genome = parse_fasta(genome_fhand)
            kmer_generator = KmerLocationGenerator(genome, 8, heterochromatic_regions)
            filtered_kmer_locations = [b'GGTAGGAT', b'TAGGATGG', b'GTAGGATG']
            kmer_locations = kmer_generator.generate_kmer_locations()
            kmer_locations = list(kmer_locations)
            filtered_kmers = filter_kmers_by_heterochromatin_stats(kmer_generator, criterion="euchromatin abundance", max_num_kmers=3)
            for filtered_kmer in filtered_kmers:
                assert filtered_kmer in filtered_kmer_locations

    def test_sort_by_euchromatin_ratio(self):
        filtered_kmers_to_test = [b'TGGTAGGG', b'GGCATATG', b'GTAGGGCA', b'AGGGCATA', b'TAGGGCAT', b'GGTAGGGC', b'GGGCATAT']
        with open(TEST_DATA_PATHDIR / 'heterochromatin_abundance.bed', 'rb') as regions_fhand, \
                open(TEST_DATA_PATHDIR / "kmers_abundance.fa", 'rb') as genome_fhand:
            heterochromatic_regions = GenomeRegions(regions_fhand)
            genome = parse_fasta(genome_fhand)
            kmer_generator = KmerLocationGenerator(genome, 8, heterochromatic_regions)
            kmer_locations = kmer_generator.generate_kmer_locations()
            kmer_locations = list(kmer_locations)
            filtered_kmers = filter_kmers_by_heterochromatin_stats(kmer_generator, criterion="euchromatin ratio", max_num_kmers=7)
            for kmer in filtered_kmers:
                assert kmer in filtered_kmers_to_test

    def test_sort_by_minimum_total_abundance(self):
        with open(TEST_DATA_PATHDIR / 'heterochromatin_abundance.bed', 'rb') as regions_fhand, \
                open(TEST_DATA_PATHDIR / "kmers_abundance.fa", 'rb') as genome_fhand:
            genome = parse_fasta(genome_fhand)
            heterochromatic_regions = GenomeRegions(regions_fhand)
            kmer_generator = KmerLocationGenerator(genome, 8, heterochromatic_regions)
            kmer_locations = kmer_generator.generate_kmer_locations()
            kmer_locations = list(kmer_locations)
            filtered_kmers = get_top_kmers_by_minimum_abundance(kmer_generator, max_num_kmers=2, min_abundance=0)
            filtered_kmer_locations = [b'ATGGTAGG', b'GCATATGT']
            for kmer in filtered_kmer_locations:
                assert kmer in filtered_kmers

            filtered_kmers = get_top_kmers_by_minimum_abundance(kmer_generator, max_num_kmers=1000, min_abundance=4)
            assert len(filtered_kmers) == 16


if __name__ == "__main__":
    unittest.main()
