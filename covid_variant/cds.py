import sys
from typing import List, Dict
from ngslite import translate, rev_comp


class Exon:

    DEL_CHAR = ' '

    start: int
    end: int
    strand: str
    __sequence: str  # LENGTH should never be changed to maintain positional information

    substitutions: Dict[int, str]
    insertions: Dict[int, str]

    def __init__(
            self,
            start: int,
            end: int,
            strand: str,
            sequence: str):
        """
        start:
            position on the + strand, 1-based inclusive

        end:
            position on the + strand, 1-based inclusive

        strand:
            '+' or '-'

        sequence:
            + strand genomic sequence
        """

        for base in set(sequence.upper()):
            assert base in ['A', 'C', 'G', 'T']
        assert end - start + 1 == len(sequence)
        assert strand in ('+', '-')

        self.start = start
        self.end = end
        self.strand = strand
        self.__sequence = sequence.upper()
        self.substitutions = dict()
        self.insertions = dict()

    @property
    def sequence(self) -> str:
        s = self.__substitute(self.__sequence)
        s = self.__insert(s)
        return s.replace(self.DEL_CHAR, '')

    def __substitute(self, sequence: str) -> str:
        s = list(sequence)
        for position, base in self.substitutions.items():
            p = position - self.start  # zero-based position on + strand
            s[p] = base
        return ''.join(s)

    def __insert(self, sequence: str) -> str:
        s = list(sequence)
        for position, bases in self.insertions.items():
            p = position - self.start  # zero-based position on + strand
            s[p] = bases + s[p]
        return ''.join(s)

    def substitute(self, position: int, base: str):
        assert self.start <= position <= self.end
        assert base.upper() in ['A', 'C', 'G', 'T']
        self.substitutions[position] = base.upper()

    def delete(self, position: int):
        assert self.start <= position <= self.end
        self.substitutions[position] = self.DEL_CHAR

    def insert(self, position: int, bases: str):
        assert self.start + 1 <= position <= self.end
        for b in set(bases):
            assert b.lower() in ['a', 'c', 'g', 't'], f'b = "{b}"'
        self.insertions[position] = bases.lower()


class CDS:

    exons: List[Exon]
    start: int
    end: int
    strand: str

    name: str

    def __init__(self, exons: List[Exon]):
        self.exons = exons
        self.start = min(e.start for e in exons)
        self.end = max(e.end for e in exons)
        self.strand = exons[0].strand
        self.check_length()

    def __str__(self):
        s = []
        for e in self.exons:
            s.append(f'{e.start}..{e.end}')
        s = ', '.join(s)
        return f"CDS ({s}) '{self.strand}'"

    @property
    def sequence(self) -> str:
        s = ''
        for e in self.exons:
            s += e.sequence
        return s

    @property
    def coding_sequence(self) -> str:
        return self.sequence if self.strand == '+' else rev_comp(self.sequence)

    def check_length(self):
        length = len(self.sequence)
        if length % 3 != 0:
            print(
                f'WARNING! {self}: coding sequence length {length}, not multiple of 3',
                file=sys.stderr,
                flush=True)

    def translate(self) -> str:
        return translate(self.coding_sequence)

    def substitute(self, position: int, base: str):
        for e in self.exons:
            if e.start <= position <= e.end:
                e.substitute(position=position, base=base)

    def delete(self, position: int):
        for e in self.exons:
            if e.start <= position <= e.end:
                e.delete(position=position)

    def insert(self, position: int, bases: str):
        for e in self.exons:
            if e.start + 1 <= position <= e.end:
                e.insert(position=position, bases=bases)
