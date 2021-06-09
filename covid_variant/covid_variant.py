import pandas as pd
from typing import List, Optional
from .cds import CDS
from .result import ReportResult
from .process_vcf import ProcessVcf
from .compare import CompareWtMutantCdses
from .template import Processor, Settings
from .pipeline import VariantCallingPipeline
from .read_gbk_mutate import ReadGbk, Mutate


class CovidVariant(Processor):

    gbk: str
    fq1: str
    fq2: Optional[str]
    covid_variant_csv: str
    tolerate_missing: float

    vcf: str
    cds_edit_df: pd.DataFrame
    wt_cdses: List[CDS]
    mutant_cdses: List[CDS]
    mutation_df: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            gbk: str,
            fq1: str,
            fq2: Optional[str],
            covid_variant_csv: str,
            tolerate_missing: float):

        self.gbk = gbk
        self.fq1 = fq1
        self.fq2 = fq2
        self.covid_variant_csv = covid_variant_csv
        self.tolerate_missing = tolerate_missing

        self.variant_calling_pipeline()
        self.process_vcf()
        self.read_gbk()
        self.mutate()
        self.compare_wt_and_mutant_cdses()
        self.report_result()

    def variant_calling_pipeline(self):
        self.vcf = VariantCallingPipeline(self.settings).main(
            gbk=self.gbk, fq1=self.fq1, fq2=self.fq2)

    def process_vcf(self):
        self.cds_edit_df = ProcessVcf(self.settings).main(vcf=self.vcf)
        self.cds_edit_df.to_csv(f'{self.outdir}/cds_edit.csv', index=False)

    def read_gbk(self):
        self.wt_cdses = ReadGbk(self.settings).main(gbk=self.gbk)

    def mutate(self):
        self.mutant_cdses = Mutate(self.settings).main(
            cdses=self.wt_cdses,
            cds_edit_df=self.cds_edit_df)

    def compare_wt_and_mutant_cdses(self):
        self.mutation_df = CompareWtMutantCdses(self.settings).main(
            wt_cdses=self.wt_cdses,
            mutant_cdses=self.mutant_cdses)
        self.mutation_df.to_csv(f'{self.outdir}/mutations.csv', index=False)

    def report_result(self):
        ReportResult(self.settings).main(
            mutation_df=self.mutation_df,
            covid_variant_df=pd.read_csv(self.covid_variant_csv),
            tolerate_missing=self.tolerate_missing)
