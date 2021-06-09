from os import makedirs
from shutil import rmtree
from os.path import exists, dirname
from .template import Settings
from .covid_variant import CovidVariant


def get_temp_path(prefix: str = 'temp', suffix: str = '') -> str:
    i = 0
    while True:
        path = f'{prefix}_{i:06d}{suffix}'
        if not exists(path):
            break
        i += 1
    return path


class Main:

    fq1: str
    fq2: str
    outdir: str
    threads: int
    debug: bool

    settings: Settings
    gbk: str
    covid_strain_csv: str

    def main(
            self,
            fq1: str,
            fq2: str,
            outdir: str,
            threads: int,
            debug: bool):

        self.fq1 = fq1
        self.fq2 = fq2
        self.outdir = outdir
        self.threads = threads
        self.debug = debug

        self.set_settings()
        self.makedirs()
        self.set_reference_paths()
        self.execute()
        self.clean_up()

    def set_settings(self):
        self.settings = Settings(
            workdir=get_temp_path(prefix='workdir'),
            outdir=self.outdir,
            threads=self.threads,
            debug=self.debug,
            mock=False)

    def makedirs(self):
        for d in [self.settings.workdir, self.settings.outdir]:
            makedirs(d, exist_ok=True)

    def set_reference_paths(self):
        ref_dir = f'{dirname(dirname(__file__))}/reference'
        self.gbk = f'{ref_dir}/NC_045512.2.gb'
        self.covid_strain_csv = f'{ref_dir}/strains.csv'

    def execute(self):
        CovidVariant(self.settings).main(
            gbk=self.gbk,
            fq1=self.fq1,
            fq2=self.fq2,
            covid_strain_csv=self.covid_strain_csv)

    def clean_up(self):
        if not self.debug:
            rmtree(self.settings.workdir)


def main(
        fq1: str,
        fq2: str,
        outdir: str,
        threads: int,
        debug: bool):

    Main().main(
        fq1=fq1,
        fq2=fq2,
        outdir=outdir,
        threads=threads,
        debug=debug)
