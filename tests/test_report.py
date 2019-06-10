import unittest
from pathlib import Path
import primer_explorer
import pickle
import io
from primer_explorer.stats import get_stats_by_pair_in_sets
from primer_explorer.report import write_detailed_report, write_stats_in_excel

TEST_DATA_PATHDIR = Path(primer_explorer.__file__).parent.parent.joinpath('tests').joinpath('data')


class TestReport(unittest.TestCase):

    def xtest_report_stats(self):
        products_path = TEST_DATA_PATHDIR / 'pcr_products.pickle'
        pcr_products_sets = pickle.load(products_path.open('rb'))
        fhand = io.StringIO()
        stats = get_stats_by_pair_in_sets(pcr_products_sets)
        write_detailed_report(fhand, stats)
        # print(fhand.getvalue())

    def test_excel_write(self):
        products_path = TEST_DATA_PATHDIR / 'pcr_products.pickle'
        pcr_products_sets = pickle.load(products_path.open('rb'))
        stats = get_stats_by_pair_in_sets(pcr_products_sets)
        write_stats_in_excel('/tmp/eee.xlsx', stats)


if __name__ == "__main__":
    unittest.main()
