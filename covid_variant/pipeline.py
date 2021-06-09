import os.path
from typing import Tuple, Optional
from ngslite import read_genbank, write_fasta
from .template import Processor, Settings


LOG_FILENAME = 'variant_calling_pipeline.log'


class VariantCallingPipeline(Processor):

    gbk: str
    fq1: str
    fq2: Optional[str]

    fna: str
    bam: str
    vcf: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, gbk: str, fq1: str, fq2: Optional[str]) -> str:

        self.gbk = gbk
        self.fq1 = fq1
        self.fq2 = fq2

        self.write_fna()
        self.trimming()
        self.mapping()
        self.variant_calling()

        return self.vcf

    def write_fna(self):
        self.fna = f'{self.workdir}/genome.fna'
        chromsome = read_genbank(file=self.gbk)[0]
        data = {chromsome.seqname: chromsome.sequence}
        write_fasta(data=data, file=self.fna)

    def trimming(self):
        if self.fq2 is None:
            self.fq1 = UnpairedTrimming(self.settings).main(
                fq=self.fq1)
        else:
            self.fq1, self.fq2 = PairedTrimming(self.settings).main(
                fq1=self.fq1, fq2=self.fq2)

    def mapping(self):
        if self.fq2 is None:
            self.bam = UnpairedMapping(self.settings).main(
                fna=self.fna, fq=self.fq1)
        else:
            self.bam = PairedMapping(self.settings).main(
                fna=self.fna, fq1=self.fq1, fq2=self.fq2)

    def variant_calling(self):
        self.vcf = VariantCalling(self.settings).main(
            fna=self.fna, bam=self.bam)


class Trimming(Processor):

    QUALITY: int = 20
    LENGTH: int = 20
    FQ_SUFFIXES = [
        '.fq',
        '.fq.gz',
        '.fastq',
        '.fastq.gz',
    ]

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)


class UnpairedTrimming(Trimming):

    fq: str

    cmd: str
    trimmed_fq: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, fq: str) -> str:
        self.fq = fq

        self.set_trimmed_fq()
        self.set_cmd()
        self.call(self.cmd)

        return self.trimmed_fq

    def set_trimmed_fq(self):
        fq = os.path.basename(self.fq)

        for suffix in self.FQ_SUFFIXES:
            if fq.endswith(suffix):
                fq = fq[:-len(suffix)]

        self.trimmed_fq = f'{self.workdir}/{fq}_trimmed.fq.gz'

    def set_cmd(self):

        self.cmd = f'''\
trim_galore \\
--quality {self.QUALITY} \\
--phred33 \\
--fastqc \\
--illumina \\
--gzip \\
--length {self.LENGTH} \\
--max_n 0 \\
--trim-n \\
--cores {self.threads} \\
--output_dir {self.workdir} \\
{self.fq} \\
1>> {self.workdir}/{LOG_FILENAME} \\
2>> {self.workdir}/{LOG_FILENAME}'''


class PairedTrimming(Trimming):

    fq1: str
    fq2: str

    cmd: str
    trimmed_fq1: str
    trimmed_fq2: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, fq1: str, fq2: str) -> Tuple[str, str]:

        self.fq1 = fq1
        self.fq2 = fq2

        self.set_trimmed_fq1_fq2()
        self.set_cmd()
        self.call(self.cmd)

        return self.trimmed_fq1, self.trimmed_fq2

    def set_trimmed_fq1_fq2(self):
        fq1 = os.path.basename(self.fq1)
        fq2 = os.path.basename(self.fq2)

        for suffix in self.FQ_SUFFIXES:
            if fq1.endswith(suffix):
                assert fq2.endswith(suffix)
                fq1 = fq1[:-len(suffix)]
                fq2 = fq2[:-len(suffix)]

        self.trimmed_fq1 = f'{self.workdir}/{fq1}_val_1.fq.gz'
        self.trimmed_fq2 = f'{self.workdir}/{fq2}_val_2.fq.gz'

    def set_cmd(self):
        self.cmd = f'''\
trim_galore \\
--paired \\
--quality {self.QUALITY} \\
--phred33 \\
--fastqc \\
--illumina \\
--gzip \\
--length {self.LENGTH} \\
--max_n 0 \\
--trim-n \\
--retain_unpaired \\
--cores {self.threads} \\
--output_dir {self.workdir} \\
{self.fq1} \\
{self.fq2} \\
1>> {self.workdir}/{LOG_FILENAME} \\
2>> {self.workdir}/{LOG_FILENAME}'''


class Mapping(Processor):

    LINE_BREAK = ' \\\n'

    fna: str
    bowtie2_index: str
    sam: str
    bam: str
    sorted_bam: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def indexing(self):
        self.bowtie2_index = f'{self.workdir}/{os.path.basename(self.fna)}'
        args = [
            'bowtie2-build',
            f'--threads {self.threads}',
            self.fna,
            self.bowtie2_index,
            f'1>> {self.workdir}/{LOG_FILENAME}',
            f'2>> {self.workdir}/{LOG_FILENAME}',
        ]
        cmd = self.LINE_BREAK.join(args)
        self.call(cmd)

    def sam_to_bam(self):
        self.bam = f'{self.workdir}/aligned.bam'
        args = [
            'samtools view -b -h',
            f'-o {self.bam}',
            self.sam
        ]
        cmd = self.LINE_BREAK.join(args)
        self.call(cmd)

    def sort_bam(self):
        self.sorted_bam = f'{self.outdir}/aligned_sorted.bam'
        args = [
            'samtools sort',
            f'-o {self.sorted_bam}',
            self.bam
        ]
        cmd = self.LINE_BREAK.join(args)
        self.call(cmd)


class UnpairedMapping(Mapping):

    fq: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, fna: str, fq: str) -> str:

        self.fna = fna
        self.fq = fq

        self.indexing()
        self.mapping()
        self.sam_to_bam()
        self.sort_bam()

        return self.sorted_bam

    def mapping(self):
        self.sam = f'{self.workdir}/aligned.sam'
        args = [
            'bowtie2',
            f'--threads {self.threads}',
            f'-x {self.bowtie2_index}',
            f'-U {self.fq}',
            f'-S {self.sam}',
            f'1>> {self.workdir}/{LOG_FILENAME}',
            f'2>> {self.workdir}/{LOG_FILENAME}',
        ]
        cmd = self.LINE_BREAK.join(args)
        self.call(cmd)


class PairedMapping(Mapping):

    fq1: str
    fq2: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, fna: str, fq1: str, fq2: str) -> str:

        self.fna = fna
        self.fq1 = fq1
        self.fq2 = fq2

        self.indexing()
        self.mapping()
        self.sam_to_bam()
        self.sort_bam()

        return self.sorted_bam

    def mapping(self):
        self.sam = f'{self.workdir}/aligned.sam'
        args = [
            'bowtie2',
            f'--threads {self.threads}',
            f'-x {self.bowtie2_index}',
            f'-1 {self.fq1}',
            f'-2 {self.fq2}',
            f'-S {self.sam}',
            f'1>> {self.workdir}/{LOG_FILENAME}',
            f'2>> {self.workdir}/{LOG_FILENAME}',
        ]
        cmd = self.LINE_BREAK.join(args)
        self.call(cmd)


class VariantCalling(Processor):

    LINE_BREAK = ' \\\n'

    fna: str
    bam: str

    vcf: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(self, fna: str, bam: str) -> str:

        self.fna = fna
        self.bam = bam

        self.variant_calling()

        return self.vcf

    def variant_calling(self):
        self.vcf = f'{self.outdir}/raw.vcf'
        args = [
            'bcftools mpileup',
            f'--threads {self.threads}',
            f'--output-type u',  # uncompressed BCF
            f'--fasta-ref {self.fna}',
            self.bam,
            f'2>> {self.workdir}/{LOG_FILENAME}',
            '|',
            'bcftools call',
            f'--threads {self.threads}',
            '--consensus-caller',
            '--variants-only',
            '--output-type v',  # uncompressed VCF
            f'-o {self.vcf}',
            f'2>> {self.workdir}/{LOG_FILENAME}'
        ]
        cmd = self.LINE_BREAK.join(args)
        self.call(cmd)
