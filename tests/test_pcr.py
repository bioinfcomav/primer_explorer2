import unittest
from collections import defaultdict

from primer_explorer.pcr import (get_pcr_products,
                                 PrimerAndLocation, annotate_reverse_complement_dimers,
                                 annotate_length_viable_products,
                                 annotate_products_by_euchromatin_region,
                                 count_total_products)
from primer_explorer.kmer import KmerAndLocation


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
