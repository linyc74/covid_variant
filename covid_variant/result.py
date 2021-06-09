import pandas as pd
from typing import List
from .template import Processor, Settings


class ReportResult(Processor):

    SEP = ','

    mutation_df: pd.DataFrame
    covid_strain_df: pd.DataFrame
    tolerate_missing: float

    spike_mutations: List[str]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            mutation_df: pd.DataFrame,
            covid_strain_df: pd.DataFrame,
            tolerate_missing: float):

        self.mutation_df = mutation_df
        self.covid_strain_df = covid_strain_df
        self.tolerate_missing = tolerate_missing

        self.set_spike_mutations()
        self.print_spike_mutations()
        self.add_match_column()
        self.print_matched_strains()

    def set_spike_mutations(self):
        df = self.mutation_df
        is_spike = df['Protein'] == 'S'
        self.spike_mutations = list(df.loc[is_spike, 'Mutation'])

    def print_spike_mutations(self):
        s = ', '.join(self.spike_mutations)
        self.print_write(f'Spike protein mutations: {s}')

    def add_match_column(self):
        self.covid_strain_df['Matched'] = \
            self.covid_strain_df['Spike Protein Substitutions'].apply(self.match)

    def match(self, s: str) -> bool:
        mutations = s.replace(' ', '').split(self.SEP)
        return list_1_in_list_2(
                list_1=mutations,
                list_2=self.spike_mutations,
                tolerated_list_1_missing_fraction=self.tolerate_missing)

    def print_matched_strains(self):
        df = self.covid_strain_df
        df = df.loc[df['Matched']]

        items = []
        for i, row in df.iterrows():
            name = row['Name']
            geo_loc = row['First Detected']
            items.append(f'{name} [{geo_loc}]')

        m = ', '.join(items)
        self.print_write(f'Match: {m}')

    def print_write(self, msg: str):
        print(msg, flush=True)
        with open(f'{self.outdir}/result.txt', 'a') as writer:
            writer.write(msg + '\n')


def list_1_in_list_2(
        list_1: List[str],
        list_2: List[str],
        tolerated_list_1_missing_fraction: float) -> bool:

    common = set(list_1).intersection(set(list_2))
    missing_fraction = abs(len(common) - len(list_1)) / len(list_1)

    return missing_fraction <= tolerated_list_1_missing_fraction
