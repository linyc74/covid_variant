import pandas as pd
from covid_variant.result import ReportResult
from .setup import TestCase


class TestReportResult(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        ReportResult(self.settings).main(
            mutation_df=pd.read_csv(f'{self.indir}/mutations.csv'),
            covid_strain_df=pd.read_csv(f'{self.indir}/strains.csv')
        )
        self.assertFileEqual(
            f'{self.indir}/result.txt',
            f'{self.outdir}/result.txt'
        )
