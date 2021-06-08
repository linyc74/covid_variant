import pandas as pd
from covid_variant.cds import CDS, Exon
from covid_variant.compare import CompareWtMutantCdses, CompareProteinPair
from .setup import TestCase


class TestCompareWtMutantCdses(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        exon1a = Exon(start=1, end=6, strand='+', sequence='AAACCC')  # KP
        exon1b = Exon(start=1, end=6, strand='+', sequence='ATGCCC')  # MP
        exon2 = Exon(start=10, end=12, strand='+', sequence='GGG')  # G

        wt_cds = CDS(exons=[exon1a, exon2])
        wt_cds.name = 'S'

        mutant_cds = CDS(exons=[exon1b])
        mutant_cds.name = 'S'

        actual = CompareWtMutantCdses(self.settings).main(
            wt_cdses=[wt_cds], mutant_cdses=[mutant_cds]
        )

        expected = pd.DataFrame(data={
            'Protein': ['S', 'S'],
            'Mutation': ['K1M', '3del'],
        })

        self.assertDataFrameEqual(expected, actual)


class TestCompareProteinPair(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        #          123  4 5678
        wt =     '-MTG--G-HKLH-'
        mutant = 'QMAGDKGYH--HE'

        actual = CompareProteinPair(self.settings).main(
            wt=wt.replace('-', ''),
            mutant=mutant.replace('-', '')
        )

        expected = ['1insQ', 'T2A', '4insDK', '5insY', '6del', '7del', '9insE']

        self.assertListEqual(expected, actual)
