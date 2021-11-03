import pandas as pd
from collections import namedtuple
from typing import List, Dict, Tuple
from Bio.pairwise2 import align, format_alignment
from .cds import CDS
from .template import Processor, Settings


class CompareWtMutantCdses(Processor):

    COLUMNS_OUT = [
        'Protein',
        'Mutation',
    ]

    wt_cdses: List[CDS]
    mutant_cdses: List[CDS]

    wt_protein_dict: Dict[str, str]
    mutant_protein_dict: Dict[str, str]
    outdf: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)
        self.compare_protein_pair = CompareProteinPair(self.settings).main

    def main(
            self,
            wt_cdses: List[CDS],
            mutant_cdses: List[CDS]) -> pd.DataFrame:

        self.wt_cdses = wt_cdses
        self.mutant_cdses = mutant_cdses

        self.set_protein_dict()
        self.init_outdf()
        self.compare_cdses()

        return self.outdf

    def set_protein_dict(self):
        self.wt_protein_dict = {
            cds.name: cds.translate() for cds in self.wt_cdses
        }
        self.mutant_protein_dict = {
            cds.name: cds.translate() for cds in self.mutant_cdses
        }

    def init_outdf(self):
        self.outdf = pd.DataFrame(columns=self.COLUMNS_OUT)

    def compare_cdses(self):
        for name in self.wt_protein_dict.keys():

            wt = self.wt_protein_dict.get(name)
            mutant = self.mutant_protein_dict.get(name, None)

            if mutant is None:
                continue

            mutations = self.compare_protein_pair(wt=wt, mutant=mutant)

            df = pd.DataFrame(data={
                'Protein': [name] * len(mutations),
                'Mutation': mutations,
            })

            self.outdf = self.outdf.append(df, ignore_index=True)


class CompareProteinPair(Processor):

    wt: str
    mutant: str

    mutations: List[str]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)
        self.__align = Align(self.settings).main
        self.__aligned_pair_to_mutations = AlignedPairToMutations(self.settings).main

    def main(self, wt: str, mutant: str) -> List[str]:
        self.wt = wt
        self.mutant = mutant

        self.align()
        self.aligned_pair_to_mutations()

        return self.mutations

    def align(self):
        self.wt, self.mutant = self.__align(wt=self.wt, mutant=self.mutant)

    def aligned_pair_to_mutations(self):
        self.mutations = self.__aligned_pair_to_mutations(
            wt=self.wt, mutant=self.mutant)


class Align(Processor):

    MATCH = 1.
    MISMATCH = -1.
    GAP_OPEN = -1.
    GAP_EXTEND = -1.

    wt: str
    mutant: str

    alignment: namedtuple
    aligned_wt: str
    aligned_mutant: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, wt: str, mutant: str) -> Tuple[str, str]:

        self.wt = wt
        self.mutant = mutant

        self.create_alignment()
        self.set_aligned_sequences()

        return self.aligned_wt, self.aligned_mutant

    def create_alignment(self):
        self.alignment = align.globalms(
            self.wt,
            self.mutant,
            self.MATCH,
            self.MISMATCH,
            self.GAP_OPEN,
            self.GAP_EXTEND)[0]

    def set_aligned_sequences(self):
        lines = format_alignment(*self.alignment).splitlines()
        self.aligned_wt = lines[0].strip()
        self.aligned_mutant = lines[2].strip()


class AlignedPairToMutations(Processor):

    wt: str
    mutant: str

    mutations: List[str]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, wt: str, mutant: str) -> List[str]:
        self.wt = wt
        self.mutant = mutant

        self.assert_equal_length()
        self.set_mutations()
        self.merge_insertions()

        return self.mutations

    def assert_equal_length(self):
        assert len(self.wt) == len(self.mutant)

    def set_mutations(self):
        self.mutations = []
        p = 1
        for wt, mut in zip(self.wt, self.mutant):
            if wt == '-':
                self.mutations.append(f'{p}ins{mut}')
            elif mut == '-':
                self.mutations.append(f'{p}del')
            elif wt != mut:
                self.mutations.append(f'{wt}{p}{mut}')

            if wt != '-':
                p += 1

    def merge_insertions(self):
        merged = []
        for mutation in self.mutations:
            if len(merged) == 0:
                merged.append(mutation)
                continue

            last = merged[-1]

            is_insertion = 'ins' in last
            if not is_insertion:
                merged.append(mutation)
                continue

            same_position = last[:4] == mutation[:4]
            if not same_position:
                merged.append(mutation)
                continue

            merged[-1] = last + mutation[-1]

        self.mutations = merged


