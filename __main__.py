import argparse
import covid_variant


__version__ = '1.0.2'


class EntryPoint:

    parser: argparse.ArgumentParser

    def main(self):
        self.set_parser()
        self.add_required_arguments()
        self.add_optional_arguments()
        args = self.parser.parse_args()
        covid_variant.main(
            fq1=args.fq1,
            fq2=args.fq2,
            outdir=args.outdir,
            threads=args.threads,
            debug=args.debug)

    def set_parser(self):
        prog = 'python covid_variant'

        description = f'Covid Variant (v{__version__}) by Yu-Cheng Lin (yclin.python@gmail.com)'

        dependencies = '\n  '.join([
            'python (>=3.6)',
            'pandas (0.25.1)',
            'ngslite (1.1.2)',
            'trim_galore (0.6.6)',
            'bowtie2 (2.4.1)',
            'samtools (1.10)',
            'bcftools (1.9)',
        ])

        epilog = f'dependency:\n  {dependencies}'

        self.parser = argparse.ArgumentParser(
            prog=prog,
            description=description,
            epilog=epilog,
            add_help=False,
            formatter_class=argparse.RawTextHelpFormatter)

    def add_required_arguments(self):
        group = self.parser.add_argument_group('required arguments')

        group.add_argument(
            '-1', '--fq1', type=str, required=True,
            help='path to read 1 fastq file')

        group.add_argument(
            '-2', '--fq2', type=str, required=True,
            help='path to read 2 fastq file')

    def add_optional_arguments(self):
        group = self.parser.add_argument_group('optional arguments')
        default = '(default: %(default)s)'

        group.add_argument(
            '-o', '--outdir', type=str, required=False, default='covid_variant_outdir',
            help=f'path to the output directory {default}')

        group.add_argument(
            '-t', '--threads', type=int, required=False, default=4,
            help=f'number of CPU threads {default}')

        group.add_argument(
            '-d', '--debug', action='store_true', help='debug mode')

        group.add_argument(
            '-h', '--help', action='help',
            help='show this help message and exit')


if __name__ == '__main__':
    EntryPoint().main()
