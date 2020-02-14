import unittest
from pathlib import Path

import primer_explorer
from primer_explorer.pcr import get_pcr_products, PrimerAndLocation
from primer_explorer.kmer import KmerAndLocation

TEST_DATA_PATHDIR = Path(primer_explorer.__file__).parent.parent.joinpath('tests').joinpath('data')


class TestPcrProducts(unittest.TestCase):

    def test_pcr_prducts(self):
        kmer_locations = {
            b'ATT': [KmerAndLocation(seq=b'ATT', chrom_location=(b'1', 0), is_heterochromatic=True),
                     KmerAndLocation(seq=b'ATT', chrom_location=(b'1', 5), is_heterochromatic=True),
                     KmerAndLocation(seq=b'ATT', chrom_location=(b'1', 50000), is_heterochromatic=True)],
            b'AAT': [KmerAndLocation(seq=b'AAT', chrom_location=(b'1', 100), is_heterochromatic=True),
                     KmerAndLocation(seq=b'AAT', chrom_location=(b'1', 20100), is_heterochromatic=True)],
            b'CGT': [KmerAndLocation(seq=b'CGT', chrom_location=(b'1', 20000), is_heterochromatic=True)],
            b'ACG': [KmerAndLocation(seq=b'ACG', chrom_location=(b'1', 120), is_heterochromatic=True),
                     KmerAndLocation(seq=b'ACG', chrom_location=(b'1', 20050), is_heterochromatic=True),
                     KmerAndLocation(seq=b'ACG', chrom_location=(b'1', 50100), is_heterochromatic=True)]
        }

        primer_combination = [b'ATT', b'CGT']
        pcr_products = get_pcr_products(kmer_locations, primer_combination)
        expected = [
            (PrimerAndLocation(strand=0, seq=b'ATT', chrom_location=(b'1', 0), is_heterochromatic=True),
             PrimerAndLocation(strand=1, seq=b'ACG', chrom_location=(b'1', 120), is_heterochromatic=True)),
            (PrimerAndLocation(strand=0, seq=b'ATT', chrom_location=(b'1', 5), is_heterochromatic=True),
             PrimerAndLocation(strand=1, seq=b'ACG', chrom_location=(b'1', 120), is_heterochromatic=True)),
            (PrimerAndLocation(strand=0, seq=b'ATT', chrom_location=(b'1', 50000), is_heterochromatic=True),
             PrimerAndLocation(strand=1, seq=b'ACG', chrom_location=(b'1', 50100), is_heterochromatic=True)),
            (PrimerAndLocation(strand=0, seq=b'CGT', chrom_location=(b'1', 20000), is_heterochromatic=True),
             PrimerAndLocation(strand=1, seq=b'AAT', chrom_location=(b'1', 20100), is_heterochromatic=True))
        ]

        assert expected == pcr_products[(b'ATT', b'CGT')]


if __name__ == "__main__":
    unittest.main()
