from covid_variant.covid_variant import CovidVariant
from .setup import TestCase


class TestCovidVariant(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        CovidVariant(self.settings).main(
            gbk=f'{self.indir}/NC_045512.2.gb',
            fq1=f'{self.indir}/54Ct21-NY-23572315_S54_L001_R1.fq.gz',
            fq2=f'{self.indir}/54Ct21-NY-23572315_S54_L001_R2.fq.gz',
            covid_strain_csv=f'{self.indir}/strains.csv'
        )
        self.assertFileEqual(
            f'{self.indir}/result.txt',
            f'{self.outdir}/result.txt'
        )
