from primer_xp import (parse_fasta, KmerLocationGenerator, GenomeRegion,
                       get_first_that_complies, GenomeRegions,
                       HETEROCHROMATIN_BED, generate_kmer_locations, KmerAndLocation, filter_kmers_by_heterochromatin_stats,
                       get_top_kmers_by_minimum_abundance)
from primer3 import kmer_validated_by_primer3
import config
import unittest


class TestKmerGenerator(unittest.TestCase):

  def test_generate_kmers(self):
    heterochromatic_regions = GenomeRegions(
        open('test_hetero.bed', 'rb'))
    genome_fhand = open(str(config.KMERS_GENERATOR_FASTA), 'rb')
    genome = parse_fasta(genome_fhand)
    kmer_generator = KmerLocationGenerator(
        genome, 8, heterochromatic_regions)
    kmer_locations = kmer_generator.generate_kmer_locations()
    list(kmer_locations)
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
    assert not GenomeRegion('c', 10, 12).overlaps(
        GenomeRegion('c', 12, 13))
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
             KmerAndLocation(b'GGGGGAGGAACAAGTCTGC',
                             (b'chrom1', 2), False),
             KmerAndLocation(b'TCCCCTTTAAAAAGGCTCAA',
                             (b'chrom1', 3), True),
             KmerAndLocation(b'TTGAGCCTTTTTAAAGGGGA', (b'chrom1', 4), False)]
    primer3_results = [kmer_validated_by_primer3(
        kmer[0].decode()) for kmer in kmers]
    assert primer3_results == [True, True, False, False]


class TestSortByEuchromatinStats(unittest.TestCase):

  def test_sort_by_euchromatin_abundance(self):
    heterochromatic_regions = GenomeRegions(
        open('test_hetero_abundance.bed', 'rb'))
    genome_fhand = open("test_kmers_abundance.fa", 'rb')
    genome = parse_fasta(genome_fhand)
    kmer_generator = KmerLocationGenerator(
        genome, 8, heterochromatic_regions)
    filtered_kmer_locations = [KmerAndLocation(seq=b'GGTAGGAT', chrom_location=(b'test_kmers', 2), is_heterochromatic=True),
                               KmerAndLocation(seq=b'GGTAGGAT', chrom_location=(
                                   b'test_kmers', 10), is_heterochromatic=False),
                               KmerAndLocation(seq=b'GGTAGGAT', chrom_location=(
                                   b'test_kmers', 18), is_heterochromatic=True),
                               KmerAndLocation(seq=b'GGTAGGAT', chrom_location=(
                                   b'test_kmers', 26), is_heterochromatic=False),
                               KmerAndLocation(seq=b'TAGGATGG', chrom_location=(
                                   b'test_kmers', 4), is_heterochromatic=True),
                               KmerAndLocation(seq=b'TAGGATGG', chrom_location=(
                                   b'test_kmers', 12), is_heterochromatic=False),
                               KmerAndLocation(seq=b'TAGGATGG', chrom_location=(
                                   b'test_kmers', 20), is_heterochromatic=True),
                               KmerAndLocation(seq=b'TAGGATGG', chrom_location=(
                                   b'test_kmers', 28), is_heterochromatic=False),
                               KmerAndLocation(seq=b'GTAGGATG', chrom_location=(
                                   b'test_kmers', 3), is_heterochromatic=True),
                               KmerAndLocation(seq=b'GTAGGATG', chrom_location=(
                                   b'test_kmers', 11), is_heterochromatic=False),
                               KmerAndLocation(seq=b'GTAGGATG', chrom_location=(
                                   b'test_kmers', 19), is_heterochromatic=True),
                               KmerAndLocation(seq=b'GTAGGATG', chrom_location=(b'test_kmers', 27), is_heterochromatic=False)]
    filtered_kmers = filter_kmers_by_heterochromatin_stats(
        kmer_generator, criterion="euchromatin abundance", max_num_kmers=1)
    for filtered_kmer in filtered_kmers:
      assert filtered_kmer in filtered_kmer_locations

  def test_sort_by_euchromatin_ratio(self):
    filtered_kmer_locations = [KmerAndLocation(seq=b'TGGTAGGG', chrom_location=(b'test_kmers', 33), is_heterochromatic=False),
                               KmerAndLocation(seq=b'GGTAGGGC', chrom_location=(
                                   b'test_kmers', 34), is_heterochromatic=False),
                               KmerAndLocation(seq=b'GTAGGGCA', chrom_location=(
                                   b'test_kmers', 35), is_heterochromatic=False),
                               KmerAndLocation(seq=b'TAGGGCAT', chrom_location=(
                                   b'test_kmers', 36), is_heterochromatic=False),
                               KmerAndLocation(seq=b'AGGGCATA', chrom_location=(
                                   b'test_kmers', 37), is_heterochromatic=False),
                               KmerAndLocation(seq=b'GGGCATAT', chrom_location=(
                                   b'test_kmers', 38), is_heterochromatic=False),
                               KmerAndLocation(seq=b'GGCATATG', chrom_location=(b'test_kmers', 39), is_heterochromatic=False)]
    heterochromatic_regions = GenomeRegions(
        open('test_hetero_abundance.bed', 'rb'))
    genome_fhand = open("test_kmers_abundance.fa", 'rb')
    genome = parse_fasta(genome_fhand)
    kmer_generator = KmerLocationGenerator(
        genome, 8, heterochromatic_regions)
    filtered_kmers = filter_kmers_by_heterochromatin_stats(
        kmer_generator, criterion="euchromatin ratio", max_num_kmers=7)
    for filtered_kmer_location in filtered_kmer_locations:
      assert filtered_kmer_location in filtered_kmers

  def test_sort_by_minimum_total_abundance(self):
    heterochromatic_regions = GenomeRegions(
        open('test_hetero_abundance.bed', 'rb'))
    genome_fhand = open("test_kmers_abundance.fa", 'rb')
    genome = parse_fasta(genome_fhand)
    kmer_generator = KmerLocationGenerator(
        genome, 8, heterochromatic_regions)
    kmer_locations = kmer_generator.generate_kmer_locations()
    kmer_locations = list(kmer_locations)
    filtered_kmers = get_top_kmers_by_minimum_abundance(kmer_locations, kmer_generator, max_num_kmers=2,
                                                        min_abundance=0)
    filtered_kmer_locations = [KmerAndLocation(seq=b'ATGGTAGG', chrom_location=(b'test_kmers', 0), is_heterochromatic=True),
                               KmerAndLocation(seq=b'ATGGTAGG', chrom_location=(
                                   b'test_kmers', 8), is_heterochromatic=True),
                               KmerAndLocation(seq=b'ATGGTAGG', chrom_location=(
                                   b'test_kmers', 16), is_heterochromatic=True),
                               KmerAndLocation(seq=b'ATGGTAGG', chrom_location=(
                                   b'test_kmers', 24), is_heterochromatic=True),
                               KmerAndLocation(seq=b'ATGGTAGG', chrom_location=(
                                   b'test_kmers', 32), is_heterochromatic=False),
                               KmerAndLocation(seq=b'GCATATGT', chrom_location=(
                                   b'test_kmers', 40), is_heterochromatic=False),
                               KmerAndLocation(seq=b'GCATATGT', chrom_location=(
                                   b'test_kmers', 48), is_heterochromatic=True),
                               KmerAndLocation(seq=b'GCATATGT', chrom_location=(
                                   b'test_kmers', 56), is_heterochromatic=True),
                               KmerAndLocation(seq=b'GCATATGT', chrom_location=(
                                   b'test_kmers', 64), is_heterochromatic=True),
                               KmerAndLocation(seq=b'GCATATGT', chrom_location=(b'test_kmers', 72), is_heterochromatic=True)]
    for kmer in filtered_kmer_locations:
      assert kmer in filtered_kmers

    filtered_kmers = get_top_kmers_by_minimum_abundance(
        kmer_locations, kmer_generator, max_num_kmers=1000, min_abundance=4)
    assert len(filtered_kmers) == 66


if __name__ == "__main__":
  unittest.main()