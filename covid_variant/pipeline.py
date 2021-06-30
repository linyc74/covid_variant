import gzip
import random
import os.path
from typing import Tuple, Optional
from ngslite import read_genbank, write_fasta
from .template import Processor, Settings


LOG_FILENAME = 'variant_calling_pipeline.log'


class VariantCallingPipeline(Processor):

    gbk: str
    fq1: str
    fq2: Optional[str]
    target_coverage: float

    fna: str
    bam: str
    vcf: str

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            gbk: str,
            fq1: str,
            fq2: Optional[str],
            target_coverage: float) -> str:

        self.gbk = gbk
        self.fq1 = fq1
        self.fq2 = fq2
        self.target_coverage = target_coverage

        self.write_fna()
        self.trimming()
        self.sampling()
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
            self.fq1 = TrimmingUnpaired(self.settings).main(
                fq=self.fq1)
        else:
            self.fq1, self.fq2 = TrimmingPaired(self.settings).main(
                fq1=self.fq1, fq2=self.fq2)

    def sampling(self):
        if self.fq2 is None:
            self.fq1 = SamplingUnpaired(self.settings).main(
                gbk=self.gbk,
                fq=self.fq1,
                target_coverage=self.target_coverage)
        else:
            self.fq1, self.fq2 = SamplingPaired(self.settings).main(
                gbk=self.gbk,
                fq1=self.fq1,
                fq2=self.fq2,
                target_coverage=self.target_coverage)

    def mapping(self):
        if self.fq2 is None:
            self.bam = MappingUnpaired(self.settings).main(
                fna=self.fna, fq=self.fq1)
        else:
            self.bam = MappingPaired(self.settings).main(
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


class TrimmingUnpaired(Trimming):

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


class TrimmingPaired(Trimming):

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


class Sampling(Processor):

    RANDOM_SEED = 1

    gbk: str
    target_coverage: float

    genome_size: int
    total_read_bases: int
    fraction: float

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)
        random.seed(self.RANDOM_SEED)

    def set_genome_size(self):
        self.genome_size = len(read_genbank(self.gbk)[0].sequence)

    def set_fraction(self):
        sample_coverage = self.total_read_bases / self.genome_size
        self.fraction = self.target_coverage / sample_coverage


class SamplingUnpaired(Sampling):

    fq: str
    sub_fq: str

    __fq: gzip.GzipFile
    __sub_fq: gzip.GzipFile

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            gbk: str,
            fq: str,
            target_coverage: float) -> str:

        self.gbk = gbk
        self.fq = fq
        self.target_coverage = target_coverage

        self.set_genome_size()
        self.set_total_read_bases()
        self.set_fraction()
        if self.fraction >= 1.:
            return self.fq
        self.set_sub_fq()
        self.random_sampling()

        return self.sub_fq

    def set_total_read_bases(self):
        t = 0
        with gzip.open(self.fq) as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    t += len(line.strip())
        self.total_read_bases = t

    def set_sub_fq(self):
        self.sub_fq = f'{self.workdir}/subsampled.fq.gz'

    def random_sampling(self):
        self.__open_files()

        end = False
        while True:
            write = random.random() <= self.fraction
            for i in range(4):
                line = self.__fq.readline()
                if write:
                    self.__sub_fq.write(line)
                if line == b'':
                    end = True
                    break
            if end:
                break

        self.__close_files()

    def __open_files(self):
        self.__fq = gzip.open(self.fq)
        self.__sub_fq = gzip.open(self.sub_fq, mode='wb')

    def __close_files(self):
        for f in [self.__fq, self.__sub_fq]:
            f.close()


class SamplingPaired(Sampling):

    fq1: str
    fq2: str
    sub_fq1: str
    sub_fq2: str

    __fq1: gzip.GzipFile
    __fq2: gzip.GzipFile
    __sub_fq1: gzip.GzipFile
    __sub_fq2: gzip.GzipFile

    def __init__(self, settings: Settings):
        super().__init__(settings=settings)

    def main(
            self,
            gbk: str,
            fq1: str,
            fq2: Optional[str],
            target_coverage: float) -> Tuple[str, str]:

        self.gbk = gbk
        self.fq1 = fq1
        self.fq2 = fq2
        self.target_coverage = target_coverage

        self.set_genome_size()
        self.set_total_read_bases()
        self.set_fraction()
        if self.fraction >= 1.:
            return self.fq1, self.fq2
        self.set_sub_fq1_fq2()
        self.random_sampling()

        return self.sub_fq1, self.sub_fq2

    def set_total_read_bases(self):
        t = 0
        for fq in [self.fq1, self.fq2]:
            with gzip.open(fq) as fh:
                for i, line in enumerate(fh):
                    if i % 4 == 1:
                        t += len(line.strip())
        self.total_read_bases = t

    def set_sub_fq1_fq2(self):
        self.sub_fq1 = f'{self.workdir}/subsampled.1.fq.gz'
        self.sub_fq2 = f'{self.workdir}/subsampled.2.fq.gz'

    def random_sampling(self):
        self.__open_files()

        end = False
        while True:
            write = random.random() <= self.fraction
            for i in range(4):
                line1 = self.__fq1.readline()
                line2 = self.__fq2.readline()
                if write:
                    self.__sub_fq1.write(line1)
                    self.__sub_fq2.write(line2)
                if line1 == b'' or line2 == b'':
                    end = True
                    break
            if end:
                break

        self.__close_files()

    def __open_files(self):
        self.__fq1 = gzip.open(self.fq1)
        self.__fq2 = gzip.open(self.fq2)
        self.__sub_fq1 = gzip.open(self.sub_fq1, mode='wb')
        self.__sub_fq2 = gzip.open(self.sub_fq2, mode='wb')

    def __close_files(self):
        for f in [self.__fq1, self.__fq2, self.__sub_fq1, self.__sub_fq2]:
            f.close()


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


class MappingUnpaired(Mapping):

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


class MappingPaired(Mapping):

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
        """
        For some reason, even though --ploidy was set to 1,
        --consensus-caller (old method) still gave me diploid calling result, e.g. "A,C"

        I set the calling method to --multiallelic-caller (new method), and the problem got resolved
        """

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
            '--multiallelic-caller',
            '--variants-only',
            '--ploidy 1',
            '--output-type v',  # uncompressed VCF
            f'-o {self.vcf}',
            f'2>> {self.workdir}/{LOG_FILENAME}'
        ]
        cmd = self.LINE_BREAK.join(args)
        self.call(cmd)
