import unittest
import config
from pathlib import Path
from collections import defaultdict
from unittest.mock import Mock

from combinations import get_compatible_groups_of_primers
from primer_xp import (parse_fasta, GenomeRegion, get_pcr_products,
                       get_first_that_complies, GenomeRegions, KmerAndLocation,
                       filter_kmers_by_heterochromatin_stats,
                       get_top_kmers_by_minimum_abundance,
                       filter_kmers_by_revcomp, KmerLocationGenerator,
                       PrimerAndLocation, annotate_reverse_complement_dimers,
                       annotate_length_viable_products,
                       annotate_products_by_euchromatin_region,
                       count_total_products)
from primer3.primer3 import kmer_validated_by_primer3, combinations_validated_by_primer3

TEST_DATA_PATHDIR = Path(__file__).parent.parent.joinpath('tests').joinpath('data')


class TestKmerGenerator(unittest.TestCase):

    def test_generate_kmers(self):
        with open(TEST_DATA_PATHDIR / 'heterochromatin_prueba.bed', 'rb') as regions_fhand, \
                open(str(config.KMERS_GENERATOR_FASTA), 'rb') as genome_fhand:
            heterochromatic_regions = GenomeRegions(regions_fhand)
            genome = parse_fasta(genome_fhand)
            kmer_generator = KmerLocationGenerator(genome, 8, heterochromatic_regions)
            kmer_locations = kmer_generator.generate_kmer_locations()
            kmer_locations = list(kmer_locations)
            assert len(kmer_generator.kmer_counters[True].most_common(60)) == 18

    def test_genome_region(self):
        region1 = GenomeRegion('chrom1', 1, 2)
        region2 = GenomeRegion('chrom2', 1, 2)
        assert region1 < region2

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


class TestPrimerCombinations(unittest.TestCase):

    def test_combination_validated_by_primer3(self):
        kmers_for_pairing = ("CAAGGGAGGAACAAGTCTGC",
                             "ATGCATAAATACAACATCG", "GATCGGATGCAATCGAGCA")
        compatibility_group_primers = [b'CCAAAAGACCCCTGCACTTA']
        for primer in kmers_for_pairing:
            if not combinations_validated_by_primer3(compatibility_group_primers, primer):
                continue
            compatibility_group_primers.append(primer)
        assert compatibility_group_primers == [b"CCAAAAGACCCCTGCACTTA", "CAAGGGAGGAACAAGTCTGC"]

        kmers_for_pairing = ("CAAGGGAGGAACAAGTCTGC", "CCAAAAGACCCCTGCACTTA",
                             "GATCGGATGCAATCGAGCA")
        compatibility_group_primers = [b'ATGCATAAATACAACATCG']
        for primer in kmers_for_pairing:
            if not combinations_validated_by_primer3(compatibility_group_primers, primer):
                continue
            compatibility_group_primers.append(primer)
        assert compatibility_group_primers == [b'ATGCATAAATACAACATCG']

    def test_make_compatibility_groups(self):
        primers = [b'TCATCAAG', b'CAACTCAA', b'GACGGACC', b'TGATGAAG',
                   b'TTGACTTG', b'CATCAACA', b'AAGAAGGA', b'CTTGAAGA',
                   b'TTCCCAAA', b'CAACATCA', b'AAGGAAGA', b'TTCTCCAA',
                   b'TTCTCCTT', b'AAAGAAGG', b'TTTCCCTT', b'TGTGATGA',
                   b'TTTGAGGA', b'AGGAAGAA', b'TTGTGATG', b'TCTTTCTC',
                   b'TGAAGTTG']
        compatibility_primers_groups = get_compatible_groups_of_primers(primers)
        assert compatibility_primers_groups[0] == [b'TCATCAAG', b'CAACTCAA', b'GACGGACC',
                                                   b'TTGACTTG', b'AAGAAGGA', b'TTCCCAAA',
                                                   b'TTCTCCTT', b'TGTGATGA', b'TCTTTCTC',
                                                   b'TGAAGTTG']
        assert compatibility_primers_groups[1] == [b'TGATGAAG', b'CAACTCAA', b'GACGGACC',
                                                   b'TTGACTTG', b'CATCAACA', b'TTCCCAAA',
                                                   b'AAGGAAGA', b'TTCTCCTT', b'TTTGAGGA',
                                                   b'TTGTGATG']


class TestPcrProducts(unittest.TestCase):

    def test_get_pcr_products(self):
        kmer_locations = {b'GGTAGGAT': [KmerAndLocation(seq=b'GGTAGGAT', chrom_location=(b'chrom1', 0), is_heterochromatic=False)],
                          b'ATCCTACC': [KmerAndLocation(seq=b'ATCCTACC', chrom_location=(b'chrom1', 1), is_heterochromatic=True),
                                        KmerAndLocation(seq=b'ATCCTACC', chrom_location=(b'chrom2', 3), is_heterochromatic=True)],
                          b'GTAGGGCA': [KmerAndLocation(seq=b'GTAGGGCA', chrom_location=(b'chrom1', 30000), is_heterochromatic=True),
                                        KmerAndLocation(seq=b'GTAGGGCA', chrom_location=(b'chrom1', 30), is_heterochromatic=True),
                                        KmerAndLocation(seq=b'GTAGGGCA', chrom_location=(b'chrom4', 40), is_heterochromatic=True)],
                          b'TGCCCTAC': [KmerAndLocation(seq=b'TGCCCTAC', chrom_location=(b'chrom2', 40), is_heterochromatic=True)],
                          b'GATGGTAG': [KmerAndLocation(seq=b'GATGGTAG', chrom_location=(b'chrom1', 40), is_heterochromatic=True)],
                          b'CTACCATC': [KmerAndLocation(seq=b'CTACCATC', chrom_location=(b'chrom4', 60), is_heterochromatic=True)]}
        primer_combination = [b'GGTAGGAT', b'GTAGGGCA', b'GATGGTAG']
        pcr_products = get_pcr_products(kmer_locations, primer_combination)
        assert pcr_products == {(b'GGTAGGAT', b'GTAGGGCA'): [(PrimerAndLocation(strand=0, seq=b'GGTAGGAT', chrom_location=(b'chrom1', 0), is_heterochromatic=False),
                                                              PrimerAndLocation(strand=1, seq=b'ATCCTACC', chrom_location=(b'chrom1', 1), is_heterochromatic=True)),
                                                             (PrimerAndLocation(strand=1, seq=b'ATCCTACC', chrom_location=(b'chrom1', 1), is_heterochromatic=True),
                                                              PrimerAndLocation(strand=0, seq=b'GTAGGGCA', chrom_location=(b'chrom1', 30), is_heterochromatic=True))],
                                (b'GGTAGGAT', b'GATGGTAG'): [(PrimerAndLocation(strand=0, seq=b'GGTAGGAT', chrom_location=(b'chrom1', 0), is_heterochromatic=False),
                                                              PrimerAndLocation(strand=1, seq=b'ATCCTACC', chrom_location=(b'chrom1', 1), is_heterochromatic=True)),
                                                             (PrimerAndLocation(strand=1, seq=b'ATCCTACC', chrom_location=(b'chrom1', 1), is_heterochromatic=True),
                                                              PrimerAndLocation(strand=0, seq=b'GATGGTAG', chrom_location=(b'chrom1', 40), is_heterochromatic=True))],
                                (b'GTAGGGCA', b'GATGGTAG'): [(PrimerAndLocation(strand=0, seq=b'GTAGGGCA', chrom_location=(b'chrom4', 40), is_heterochromatic=True),
                                                              PrimerAndLocation(strand=1, seq=b'CTACCATC', chrom_location=(b'chrom4', 60), is_heterochromatic=True))]}

    def test_pcr_products_annotation(self):
        annotations = defaultdict(dict)
        pcr_products = {(b'GGTAGGAT', b'GTAGGGCA'): [(PrimerAndLocation(strand=0, seq=b'GGTAGGAT', chrom_location=(b'chrom1', 0), is_heterochromatic=False),
                                                      PrimerAndLocation(strand=1, seq=b'ATCCTACC', chrom_location=(b'chrom1', 1), is_heterochromatic=False)),
                                                     (PrimerAndLocation(strand=1, seq=b'ATCCTACC', chrom_location=(b'chrom1', 1), is_heterochromatic=True),
                                                      PrimerAndLocation(strand=0, seq=b'GTAGGGCA', chrom_location=(b'chrom1', 30), is_heterochromatic=False))],
                        (b'GGTAGGAT', b'GATGGTAG'): [(PrimerAndLocation(strand=0, seq=b'GGTAGGAT', chrom_location=(b'chrom1', 0), is_heterochromatic=True),
                                                      PrimerAndLocation(strand=1, seq=b'ATCCTACC', chrom_location=(b'chrom1', 900), is_heterochromatic=True)),
                                                     (PrimerAndLocation(strand=1, seq=b'ATCCTACC', chrom_location=(b'chrom1', 1), is_heterochromatic=True),
                                                      PrimerAndLocation(strand=0, seq=b'GATGGTAG', chrom_location=(b'chrom1', 40), is_heterochromatic=True))],
                        (b'GTAGGGCA', b'GATGGTAG'): [(PrimerAndLocation(strand=0, seq=b'GTAGGGCA', chrom_location=(b'chrom4', 0), is_heterochromatic=False),
                                                      PrimerAndLocation(strand=1, seq=b'CTACCATC', chrom_location=(b'chrom4', 350), is_heterochromatic=True))]}

        for primers, products in pcr_products.items():
            annotations[primers]["total_products"] = count_total_products(products)
        assert annotations[(b'GGTAGGAT', b'GTAGGGCA')]["total_products"] == 2
        assert annotations[(b'GGTAGGAT', b'GATGGTAG')]["total_products"] == 2
        assert annotations[(b'GTAGGGCA', b'GATGGTAG')]["total_products"] == 1

        for primers, products in pcr_products.items():
            annotations[primers]["pcr_revcomp_dimers"] = annotate_reverse_complement_dimers(products)
        assert annotations[(b'GGTAGGAT', b'GTAGGGCA')]["pcr_revcomp_dimers"][True] == 1
        assert annotations[(b'GGTAGGAT', b'GTAGGGCA')]["pcr_revcomp_dimers"][False] == 1
        assert annotations[(b'GGTAGGAT', b'GATGGTAG')]["pcr_revcomp_dimers"][True] == 1
        assert annotations[(b'GGTAGGAT', b'GATGGTAG')]["pcr_revcomp_dimers"][False] == 1
        assert annotations[(b'GTAGGGCA', b'GATGGTAG')]["pcr_revcomp_dimers"][True] == 0
        assert annotations[(b'GTAGGGCA', b'GATGGTAG')]["pcr_revcomp_dimers"][False] == 1

        for primers, products in pcr_products.items():
            annotations[primers]["length_viable_products"] = annotate_length_viable_products(products)
        assert annotations[(b'GGTAGGAT', b'GTAGGGCA')]["length_viable_products"]["viable"] == 0
        assert annotations[(b'GGTAGGAT', b'GTAGGGCA')]["length_viable_products"]["inviable"] == 2
        assert annotations[(b'GGTAGGAT', b'GATGGTAG')]["length_viable_products"]["viable"] == 1
        assert annotations[(b'GGTAGGAT', b'GATGGTAG')]["length_viable_products"]["inviable"] == 1
        assert annotations[(b'GTAGGGCA', b'GATGGTAG')]["length_viable_products"]["viable"] == 1
        assert annotations[(b'GTAGGGCA', b'GATGGTAG')]["length_viable_products"]["inviable"] == 0

        for primers, products in pcr_products.items():
            annotations[primers]["euchromatin_products"] = annotate_products_by_euchromatin_region(products)
        assert annotations[(b'GGTAGGAT', b'GTAGGGCA')]["euchromatin_products"]["euchromatic"] == 1
        assert annotations[(b'GGTAGGAT', b'GTAGGGCA')]["euchromatin_products"]["heterochromatic"] == 0
        assert annotations[(b'GGTAGGAT', b'GTAGGGCA')]["euchromatin_products"]["mixed"] == 1
        assert annotations[(b'GGTAGGAT', b'GATGGTAG')]["euchromatin_products"]["euchromatic"] == 0
        assert annotations[(b'GGTAGGAT', b'GATGGTAG')]["euchromatin_products"]["heterochromatic"] == 2
        assert annotations[(b'GGTAGGAT', b'GATGGTAG')]["euchromatin_products"]["mixed"] == 0
        assert annotations[(b'GTAGGGCA', b'GATGGTAG')]["euchromatin_products"]["euchromatic"] == 0
        assert annotations[(b'GTAGGGCA', b'GATGGTAG')]["euchromatin_products"]["heterochromatic"] == 0
        assert annotations[(b'GTAGGGCA', b'GATGGTAG')]["euchromatin_products"]["mixed"] == 1


if __name__ == "__main__":
    unittest.main()
