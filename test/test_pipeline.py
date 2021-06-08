from covid_variant.pipeline import VariantCallingPipeline
from .setup import TestCase


class TestVariantCallingPipeline(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = VariantCallingPipeline(self.settings).main(
            gbk=f'{self.indir}/NC_045512.2.gb',
            fq1=f'{self.indir}/54Ct21-NY-23572315_S54_L001_R1.fq.gz',
            fq2=f'{self.indir}/54Ct21-NY-23572315_S54_L001_R2.fq.gz'
        )
        expected = f'{self.indir}/raw.vcf'
        with open(expected) as fh1:
            with open(actual) as fh2:
                for line1, line2 in zip(fh1, fh2):
                    if line1.startswith('##bcftools_callCommand'):  # Skip the time stamp in this line
                        continue
                    self.assertEqual(line1, line2)
