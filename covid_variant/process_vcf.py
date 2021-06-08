import pandas as pd
from typing import Union
from .template import Processor, Settings


class ProcessVcf(Processor):

    vcf: str

    vcf_df: pd.DataFrame
    cds_edit_df: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, vcf: str) -> pd.DataFrame:
        self.vcf = vcf

        self.read_vcf()
        self.remove_conflict_variants()
        self.vcf_df_to_cds_edit_df()

        return self.cds_edit_df

    def read_vcf(self):
        self.vcf_df = ReadVcf(self.settings).main(vcf=self.vcf)

    def remove_conflict_variants(self):
        self.vcf_df = RemoveConflictVariants(self.settings).main(indf=self.vcf_df)

    def vcf_df_to_cds_edit_df(self):
        self.cds_edit_df = VcfDfToCdsEditDf(self.settings).main(vcf_df=self.vcf_df)


class ReadVcf(Processor):

    vcf: str

    temp_tsv: str
    outdf: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, vcf: str) -> pd.DataFrame:
        self.vcf = vcf

        self.write_temp_tsv()
        self.read_temp_tsv()

        return self.outdf

    def write_temp_tsv(self):
        self.temp_tsv = f'{self.workdir}/temp.tsv'
        with open(self.vcf) as reader:
            with open(self.temp_tsv, 'w') as writer:
                for line in reader:
                    if line.startswith('##'):
                        continue
                    elif line.startswith('#CHROM'):
                        writer.write(line[1:])
                    else:
                        writer.write(line)

    def read_temp_tsv(self):
        self.outdf = pd.read_csv(self.temp_tsv, sep='\t')


class RemoveConflictVariants(Processor):

    COLUMNS_IN = [
        'POS',
        'REF',
        'ALT',
        'QUAL',
    ]
    COLUMNS_OUT = COLUMNS_IN

    indf: pd.DataFrame
    snv_df: pd.DataFrame
    del_df: pd.DataFrame
    ins_df: pd.DataFrame
    outdf: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, indf: pd.DataFrame) -> pd.DataFrame:
        self.indf = indf

        self.separate_three_dfs()
        self.remove_snv_conflict()
        self.remove_deletion_conflict()
        self.remove_insertion_conflict()
        self.merge_back()

        return self.outdf

    def separate_three_dfs(self):
        df = self.indf
        ref_len = df['REF'].apply(len)
        alt_len = df['ALT'].apply(len)
        self.snv_df = df[ref_len == alt_len].reset_index(drop=True)
        self.del_df = df[ref_len > alt_len].reset_index(drop=True)
        self.ins_df = df[ref_len < alt_len].reset_index(drop=True)

    def remove_snv_conflict(self):
        self.snv_df = RemoveSnvConflict(self.settings).main(indf=self.snv_df)

    def remove_deletion_conflict(self):
        if len(self.del_df) > 0:
            self.del_df = RemoveDeletionConflict(self.settings).main(indf=self.del_df)

    def remove_insertion_conflict(self):
        if len(self.ins_df) > 0:
            self.ins_df = RemoveInsertionConflict(self.settings).main(indf=self.ins_df)

    def merge_back(self):
        self.outdf = self.snv_df.append(
            self.del_df, ignore_index=True
        ).append(
            self.ins_df, ignore_index=True
        ).sort_values(
            by='POS',
            ascending=True
        ).reset_index(
            drop=True
        )


class RemoveSnvConflict(Processor):

    indf: pd.DataFrame
    outdf: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, indf: pd.DataFrame) -> pd.DataFrame:
        self.indf = indf

        self.outdf = self.indf.sort_values(
            by='QUAL',
            ascending=False
        ).drop_duplicates(
            subset='POS',
            keep='first'
        ).sort_values(
            by='POS',
            ascending=True
        ).reset_index(
            drop=True
        )

        return self.outdf


class RemoveDeletionConflict(Processor):

    COLUMNS_IN = [
        'POS',
        'REF',
        'ALT',
        'QUAL',
    ]
    COLUMNS_OUT = COLUMNS_IN

    indf: pd.DataFrame
    outdf: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, indf: pd.DataFrame) -> pd.DataFrame:
        self.indf = indf
        self.outdf = pd.DataFrame(columns=self.COLUMNS_OUT)

        self.sort_by_position()
        self.add_deletion_range()
        self.remove_overlapped()
        self.drop_deletion_range()

        return self.outdf

    def sort_by_position(self):
        self.indf = self.indf.sort_values(
            by='POS',
            ascending=True
        ).reset_index(
            drop=True
        )

    def add_deletion_range(self):
        df = self.indf
        for i, row in df.iterrows():
            pos, ref, alt = row['POS'], row['REF'], row['ALT']
            df.loc[i, 'start'] = pos + len(alt)
            df.loc[i, 'end'] = pos + len(ref) - 1  # -1 because the end is inclusive

    def remove_overlapped(self):
        df = self.indf
        prev = None
        for i, row in df.iterrows():
            this = row

            if prev is None:
                prev = this
                continue

            if overlap(prev, this):
                if this['QUAL'] > prev['QUAL']:
                    prev = this
                else:
                    pass  # Just skip this deletion (lower quality), and keep the previous one
            else:
                self.outdf = self.outdf.append(prev, ignore_index=True)
                prev = this

        self.outdf = self.outdf.append(prev, ignore_index=True)  # Append the final one

    def drop_deletion_range(self):
        self.outdf = self.outdf.drop(columns=['start', 'end'])


def overlap(row1: pd.Series, row2: pd.Series) -> bool:
    return (row2['start'] <= row1['end']) and (row1['start'] <= row2['end'])


class RemoveInsertionConflict(Processor):

    COLUMNS_IN = [
        'POS',
        'REF',
        'ALT',
        'QUAL',
    ]
    COLUMNS_OUT = COLUMNS_IN

    indf: pd.DataFrame
    outdf: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, indf: pd.DataFrame) -> pd.DataFrame:
        self.indf = indf
        self.outdf = pd.DataFrame(columns=self.COLUMNS_OUT)

        self.sort_by_position()
        self.add_insertion_position()
        self.remove_conflict()
        self.drop_insertion_position()

        return self.outdf

    def sort_by_position(self):
        self.indf = self.indf.sort_values(
            by='POS',
            ascending=True
        ).reset_index(
            drop=True
        )

    def add_insertion_position(self):
        df = self.indf
        for i, row in df.iterrows():
            pos, ref = row['POS'], row['REF']
            df.loc[i, 'position'] = pos + len(ref)

    def remove_conflict(self):
        self.outdf = self.indf.sort_values(
            by='QUAL',
            ascending=False
        ).drop_duplicates(
            subset='position',
            keep='first'
        ).sort_values(
            by='POS',
            ascending=True
        ).reset_index(
            drop=True
        )

    def drop_insertion_position(self):
        self.outdf = self.outdf.drop(columns='position')


class VcfDfToCdsEditDf(Processor):

    COLUMNS_IN = [
        'POS',
        'REF',
        'ALT',
    ]
    COLUMNS_OUT = [
        'Position',
        'Type',
        'Base',
    ]

    vcf_df: pd.DataFrame
    cds_edit_df: pd.DataFrame

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        self.vcf_df = vcf_df
        self.cds_edit_df = pd.DataFrame(columns=self.COLUMNS_OUT)

        self.sort_by_position()
        for i, row in self.vcf_df.iterrows():
            self.convert_one(row)

        return self.cds_edit_df

    def sort_by_position(self):
        self.vcf_df = self.vcf_df.sort_values(
            by='POS',
            ascending=True
        ).reset_index(
            drop=True
        )

    def convert_one(self, row: pd.Series):
        pos, ref, alt = row['POS'], row['REF'], row['ALT']
        is_snv = len(ref) == len(alt)
        is_deletion = len(ref) > len(alt)
        if is_snv:
            self.convert_snv(row)
        elif is_deletion:
            self.convert_deletion(row)
        else:
            self.convert_insertion(row)

    def convert_snv(self, row: pd.Series):
        self.append({
            'Position': row['POS'],
            'Type': 'substitute',
            'Base': row['ALT'],
        })

    def convert_deletion(self, row: pd.Series):
        pos, ref, alt = row['POS'], row['REF'], row['ALT']

        start = pos + len(alt)
        length = len(ref) - len(alt)
        positions = [(start + i) for i in range(length)]

        new = pd.DataFrame(data={
            'Position': positions,
            'Type': ['delete'] * len(positions),
        })
        self.append(new)

    def convert_insertion(self, row: pd.Series):
        pos, ref, alt = row['POS'], row['REF'], row['ALT']

        insertion_site = pos + len(ref)
        bases = alt[len(ref):]

        self.append({
            'Position': insertion_site,
            'Type': 'insert',
            'Base': bases,
        })

    def append(self, new: Union[dict, pd.DataFrame]):
        self.cds_edit_df = self.cds_edit_df.append(new, ignore_index=True)
