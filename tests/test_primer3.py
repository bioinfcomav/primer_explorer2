import unittest

from primer3.primer3 import combinations_validated_by_primer3
from combinations import get_compatible_groups_of_primers


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


if __name__ == "__main__":
    unittest.main()
