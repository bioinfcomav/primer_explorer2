import unittest
from regions import GenomeRegion


class TestKmerGenerator(unittest.TestCase):

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


if __name__ == "__main__":
    unittest.main()
