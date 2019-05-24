import unittest
from pathlib import Path
import primer_explorer
import pickle
import io
from primer_explorer.stats import get_stats_by_pair_in_sets
from primer_explorer.report import write_detailed_report

TEST_DATA_PATHDIR = Path(primer_explorer.__file__).parent.parent.joinpath('tests').joinpath('data')


class TestReport(unittest.TestCase):

    def test_report_stats(self):
        producs_fpath = TEST_DATA_PATHDIR / 'pcr_products.pickle'
        pcr_products_sets = pickle.load(open(producs_fpath, 'rb'))
        fhand = io.StringIO()
        stats = get_stats_by_pair_in_sets(pcr_products_sets)
        write_detailed_report(fhand, stats)
        # print(fhand.getvalue())


if __name__ == "__main__":
    unittest.main()
