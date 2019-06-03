import unittest
from pathlib import Path

import primer_explorer
from primer_explorer.regions import GenomeRegion, GenomeRegions

TEST_DATA_PATHDIR = Path(primer_explorer.__file__).parent.parent.joinpath('tests').joinpath('data')


class TestKmerGenerator(unittest.TestCase):

    def test_genome_region(self):
        region1 = GenomeRegion('chrom1', 5, 10)
        region2 = GenomeRegion('chrom2', 5, 10)
        assert region1 < region2

        assert region1 < GenomeRegion('chrom1', 6, 13)
        assert not region1 < GenomeRegion('chrom1', 2, 3)
        assert not region1 < GenomeRegion('chrom0', 12, 3)

        assert not GenomeRegion('c', 0, 2).overlaps(GenomeRegion('c', 2, 4))
        assert not GenomeRegion('c', 2, 4).overlaps(GenomeRegion('c', 0, 2))
        assert not GenomeRegion('c', 0, 2).overlaps(GenomeRegion('c', 3, 4))
        assert not GenomeRegion('c', 3, 4).overlaps(GenomeRegion('c', 0, 2))
        assert GenomeRegion('c', 0, 2).overlaps(GenomeRegion('c', 0, 1))
        assert GenomeRegion('c', 10, 12).overlaps(GenomeRegion('c', 9, 11))
        assert GenomeRegion('c', 10, 12).overlaps(GenomeRegion('c', 11, 13))
        assert not GenomeRegion('c', 10, 12).overlaps(GenomeRegion('c', 12, 13))
        assert GenomeRegion('c', 10, 15).overlaps(GenomeRegion('c', 11, 13))

    def test_genome_regions(self):
        with (TEST_DATA_PATHDIR / 'heterochromatin.bed').open() as bed_fhand:
            regions = list(GenomeRegions(bed_fhand=bed_fhand))
            assert len(regions) == 2
            assert regions[0].chrom == 'test_kmers'
            assert regions[0].start == 0
            assert regions[0].stop == 5


if __name__ == "__main__":
    unittest.main()
