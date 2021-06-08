import pandas as pd
from covid_variant.process_vcf import ProcessVcf, RemoveConflictVariants, VcfDfToCdsEditDf, ReadVcf
from .setup import TestCase


class TestProcessVcf(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ProcessVcf(self.settings).main(vcf=f'{self.indir}/{self.__class__.__name__}_in.vcf')
        expected = pd.read_csv(f'{self.indir}/{self.__class__.__name__}_out.csv')
        self.assertDataFrameEqual(expected, actual)


class TestReadVcf(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = ReadVcf(self.settings).main(vcf=f'{self.indir}/{self.__class__.__name__}_in.vcf')
        expected = pd.read_csv(f'{self.indir}/{self.__class__.__name__}_out.csv')
        self.assertDataFrameEqual(expected, actual)


class TestRemoveConflictVariants(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        indf = pd.read_csv(f'{self.indir}/{self.__class__.__name__}_in.csv')
        actual = RemoveConflictVariants(self.settings).main(indf=indf)
        expected = pd.read_csv(f'{self.indir}/{self.__class__.__name__}_out.csv')
        self.assertDataFrameEqual(expected, actual)


class TestVcfDfToCdsEditEf(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        vcf_df = pd.read_csv(f'{self.indir}/{self.__class__.__name__}_in.csv')
        actual = VcfDfToCdsEditDf(self.settings).main(vcf_df=vcf_df)
        expected = pd.read_csv(f'{self.indir}/{self.__class__.__name__}_out.csv')
        self.assertDataFrameEqual(expected, actual)
