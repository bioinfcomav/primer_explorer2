import unittest
from pathlib import Path

import primer_explorer
from primer_explorer.pcr import get_pcr_products, PrimerAndLocation
from primer_explorer.kmer import KmerAndLocation

TEST_DATA_PATHDIR = Path(primer_explorer.__file__).parent.parent.joinpath('tests').joinpath('data')


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
                                                              PrimerAndLocation(strand=1, seq=b'ATCCTACC', chrom_location=(b'chrom1', 1), is_heterochromatic=True))],
                                (b'GGTAGGAT', b'GATGGTAG'): [(PrimerAndLocation(strand=0, seq=b'GGTAGGAT', chrom_location=(b'chrom1', 0), is_heterochromatic=False),
                                                              PrimerAndLocation(strand=1, seq=b'ATCCTACC', chrom_location=(b'chrom1', 1), is_heterochromatic=True))],
                                (b'GTAGGGCA', b'GATGGTAG'): [(PrimerAndLocation(strand=0, seq=b'GTAGGGCA', chrom_location=(b'chrom4', 40), is_heterochromatic=True),
                                                              PrimerAndLocation(strand=1, seq=b'CTACCATC', chrom_location=(b'chrom4', 60), is_heterochromatic=True))]}


if __name__ == "__main__":
    unittest.main()
