from unittest import TestCase
from covid_variant.cds import Exon, CDS


class TestExon(TestCase):

    def test_wrong_nucleotide(self):
        with self.assertRaises(AssertionError):
            Exon(start=1, end=1, strand='+', sequence='N')

    def test_wrong_sequence_length(self):
        with self.assertRaises(AssertionError):
            Exon(start=1, end=5, strand='+', sequence='A')


class TestSubstitute(TestCase):

    STRAND = '-'  # Minus strand should not matter, always operate on the plus strand

    def test_single_base(self):
        exon = Exon(start=11, end=11, strand=self.STRAND, sequence='A')
        exon.substitute(position=11, base='T')
        self.assertEqual('T', exon.sequence)

    def test_first_of_two_bases(self):
        exon = Exon(start=11, end=12, strand=self.STRAND, sequence='AA')
        exon.substitute(position=11, base='T')
        self.assertEqual('TA', exon.sequence)

    def test_second_of_two_bases(self):
        exon = Exon(start=11, end=12, strand=self.STRAND, sequence='AA')
        exon.substitute(position=12, base='T')
        self.assertEqual('AT', exon.sequence)

    def test_first_of_many_bases(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='AAAAAA')
        exon.substitute(position=11, base='T')
        self.assertEqual('TAAAAA', exon.sequence)

    def test_middle_of_many_bases(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='AAAAAA')
        exon.substitute(position=12, base='T')
        self.assertEqual('ATAAAA', exon.sequence)

    def test_last_of_many_bases(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='AAAAAA')
        exon.substitute(position=16, base='T')
        self.assertEqual('AAAAAT', exon.sequence)

    def test_multiple(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='AAAAAA')
        for position, base in [
            (12, 'C'),
            (13, 'T'),
            (14, 'G'),
        ]:
            exon.substitute(position=position, base=base)
        self.assertEqual('ACTGAA', exon.sequence)


class TestDelete(TestCase):

    STRAND = '-'

    def test_first_base(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='ATAAAA')
        exon.delete(position=11)
        self.assertEqual('TAAAA', exon.sequence)

    def test_middle_base(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='ATAAAA')
        exon.delete(position=12)
        self.assertEqual('AAAAA', exon.sequence)

    def test_last_base(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='ATAAAA')
        exon.delete(position=16)
        self.assertEqual('ATAAA', exon.sequence)


class TestInsert(TestCase):

    STRAND = '-'

    def test_out_of_start_bound(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='AAAAAA')
        with self.assertRaises(AssertionError):
            exon.insert(position=11, bases='T')

    def test_out_of_end_bound(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='AAAAAA')
        with self.assertRaises(AssertionError):
            exon.insert(position=17, bases='T')

    def test_second_of_many_bases(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='AAAAAA')
        exon.insert(position=12, bases='TT')
        self.assertEqual('AttAAAAA', exon.sequence)

    def test_last_of_many_bases(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='AAAAAA')
        exon.insert(position=16, bases='TT')
        self.assertEqual('AAAAAttA', exon.sequence)


class TestCompoundMutation(TestCase):

    STRAND = '-'

    def test_substitute_delete_on_same_position(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='ATAAAA')
        exon.substitute(position=12, base='C')
        exon.delete(position=12)
        self.assertEqual('AAAAA', exon.sequence)

    def test_substitute_delete_on_different_positions(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='ATAAAA')
        exon.substitute(position=13, base='C')
        exon.delete(position=12)
        self.assertEqual('ACAAA', exon.sequence)

    def test_delete_insert_on_same_position(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='ATAAAA')
        exon.delete(position=12)
        exon.insert(position=12, bases='C')
        self.assertEqual('AcAAAA', exon.sequence)

    def test_delete_insert_on_different_position(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='ATAAAA')
        exon.delete(position=12)
        exon.insert(position=16, bases='C')
        self.assertEqual('AAAAcA', exon.sequence)

    def test_substitute_insert_on_same_position(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='ATAAAA')
        exon.substitute(position=12, base='C')
        exon.insert(position=12, bases='G')
        self.assertEqual('AgCAAAA', exon.sequence)

    def test_substitute_insert_on_different_position(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='ATAAAA')
        exon.substitute(position=12, base='C')
        exon.insert(position=13, bases='G')
        self.assertEqual('ACgAAAA', exon.sequence)

    def test_substitute_delete_insert(self):
        exon = Exon(start=11, end=16, strand=self.STRAND, sequence='ATAAAA')
        exon.substitute(position=12, base='C')
        exon.delete(position=13)
        exon.insert(position=14, bases='GG')
        self.assertEqual('ACggAAA', exon.sequence)


class TestCDS(TestCase):

    def test_sequence(self):
        exon1 = Exon(start=1, end=7, strand='-', sequence='AAACCCg')
        exon2 = Exon(start=11, end=12, strand='-', sequence='gg')
        cds = CDS(exons=[exon1, exon2])
        self.assertEqual('AAACCCGGG', cds.sequence)

    def test_coding_sequence(self):
        exon1 = Exon(start=1, end=7, strand='-', sequence='AAACCCg')
        exon2 = Exon(start=11, end=12, strand='-', sequence='gg')
        cds = CDS(exons=[exon1, exon2])
        self.assertEqual('CCCGGGTTT', cds.coding_sequence)

    def test_translate_plus_strand(self):
        exon1 = Exon(start=1, end=7, strand='+', sequence='AAACCCg')
        exon2 = Exon(start=11, end=12, strand='+', sequence='gg')
        cds = CDS(exons=[exon1, exon2])
        self.assertEqual('KPG', cds.translate())

    def test_translate_minus_strand(self):
        exon1 = Exon(start=1, end=7, strand='-', sequence='AAACCCg')
        exon2 = Exon(start=11, end=12, strand='-', sequence='gg')
        cds = CDS(exons=[exon1, exon2])
        self.assertEqual('PGF', cds.translate())

    def test_substitute(self):
        exon1 = Exon(start=1, end=7, strand='+', sequence='AAACCCg')
        exon2 = Exon(start=11, end=12, strand='+', sequence='gg')
        cds = CDS(exons=[exon1, exon2])
        cds.substitute(position=2, base='T')
        cds.substitute(position=3, base='G')
        self.assertEqual('ATGCCCGGG', cds.sequence)

    def test_delete_insert(self):
        exon1 = Exon(start=1, end=7, strand='+', sequence='AAACCCg')
        exon2 = Exon(start=11, end=12, strand='+', sequence='gg')
        cds = CDS(exons=[exon1, exon2])
        cds.delete(position=1)
        cds.delete(position=2)
        cds.delete(position=3)
        cds.insert(position=2, bases='ATG')
        self.assertEqual('atgCCCGGG', cds.sequence)
