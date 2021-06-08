import pandas as pd
from typing import List
from copy import deepcopy
from ngslite import read_genbank, Chromosome, GenericFeature
from .cds import CDS, Exon
from .template import Processor, Settings


class ReadGbk(Processor):

    gbk: str

    chromosome: Chromosome
    cdses: List[CDS]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, gbk: str) -> List[CDS]:
        self.gbk = gbk
        self.read_gbk()
        self.features_to_cdses()
        return self.cdses

    def read_gbk(self):
        self.chromosome = read_genbank(file=self.gbk)[0]

    def features_to_cdses(self):
        self.cdses = []
        features = [f for f in self.chromosome.features if f.type == 'CDS']
        self.cdses = list(map(self.to_cds, features))

    def to_cds(self, feature: GenericFeature) -> CDS:
        exons = []
        for start, end, strand in feature.regions:
            exon = Exon(
                start=start,
                end=end,
                strand=strand,
                sequence=self.chromosome.sequence[start - 1:end]
            )
            exons.append(exon)
        cds = CDS(exons=exons)
        cds.name = feature.get_attribute(key='gene')
        return cds


TYPE = 'Type'
POSITION = 'Position'
BASE = 'Base'


class Mutate(Processor):

    cdses: List[CDS]
    cds_edit_df: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)
        self.mutate_one_cds = MutateOneCDS(self.settings).main

    def main(
            self,
            cdses: List[CDS],
            cds_edit_df: pd.DataFrame) -> List[CDS]:

        self.cdses = deepcopy(cdses)
        self.cds_edit_df = cds_edit_df

        for cds in self.cdses:
            self.mutate_one_cds(cds=cds, cds_edit_df=cds_edit_df)

        return self.cdses


class MutateOneCDS(Processor):

    cds: CDS
    cds_edit_df: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            cds: CDS,
            cds_edit_df: pd.DataFrame):

        self.cds = cds
        self.cds_edit_df = cds_edit_df

        self.subset_cds_edit_df()
        self.substitute()
        self.delete()
        self.insert()

        return self.cds

    def subset_cds_edit_df(self):
        within_cds = \
            (self.cds_edit_df[POSITION] >= self.cds.start) & \
            (self.cds_edit_df[POSITION] <= self.cds.end)
        self.cds_edit_df = self.cds_edit_df.loc[within_cds]

    def substitute(self):
        df = self.cds_edit_df
        for i, row in df[df[TYPE] == 'substitute'].iterrows():
            self.cds.substitute(position=row[POSITION], base=row[BASE])

    def delete(self):
        df = self.cds_edit_df
        for i, row in df[df[TYPE] == 'delete'].iterrows():
            self.cds.delete(position=row[POSITION])

    def insert(self):
        df = self.cds_edit_df
        for i, row in df[df[TYPE] == 'insert'].iterrows():
            first_base = row[POSITION] == self.cds.start
            if first_base:
                continue  # by definition insertion cannot take place before the first base
            self.cds.insert(position=row[POSITION], bases=row[BASE])
