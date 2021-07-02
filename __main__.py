import argparse
import covid_variant


__version__ = '1.1.1-beta'


PROG = 'python covid_variant'
DESCRIPTION = f'Covid Variant (version {__version__}) by Yu-Cheng Lin (linyc74@gmail.com)'
REQUIRED = [
    {
        'keys': ['-1', '--fq1'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to read 1 fastq file',
        }
    },
]
OPTIONAL = [
    {
        'keys': ['-2', '--fq2'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'path to read 2 fastq file (default: %(default)s)',
        }
    },
    {
        'keys': ['-o', '--outdir'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'covid_variant_outdir',
            'help': 'path to the output directory (default: %(default)s)',
        }
    },
    {
        'keys': ['--tolerate-missing'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.1,
            'help': 'fraction of missing mutations to be tolerated when matching to known variants (default: %(default)s)',
        }
    },
    {
        'keys': ['--target-coverage'],
        'properties': {
            'type': float,
            'required': False,
            'default': float('inf'),
            'help': 'subsample fastq reads to the target coverage (default: %(default)s)',
        }
    },
    {
        'keys': ['-t', '--threads'],
        'properties': {
            'type': int,
            'required': False,
            'default': 4,
            'help': 'number of CPU threads (default: %(default)s)',
        }
    },
    {
        'keys': ['-d', '--debug'],
        'properties': {
            'action': 'store_true',
            'help': 'debug mode',
        }
    },
    {
        'keys': ['-h', '--help'],
        'properties': {
            'action': 'help',
            'help': 'show this help message',
        }
    },
    {
        'keys': ['-v', '--version'],
        'properties': {
            'action': 'version',
            'version': __version__,
            'help': 'show version',
        }
    },
]


class EntryPoint:

    parser: argparse.ArgumentParser

    def main(self):
        self.set_parser()
        self.add_required_arguments()
        self.add_optional_arguments()
        self.run()

    def set_parser(self):
        self.parser = argparse.ArgumentParser(
            prog=PROG,
            description=DESCRIPTION,
            add_help=False,
            formatter_class=argparse.RawTextHelpFormatter)

    def add_required_arguments(self):
        group = self.parser.add_argument_group('required arguments')
        for item in REQUIRED:
            group.add_argument(*item['keys'], **item['properties'])

    def add_optional_arguments(self):
        group = self.parser.add_argument_group('optional arguments')
        for item in OPTIONAL:
            group.add_argument(*item['keys'], **item['properties'])

    def run(self):
        args = self.parser.parse_args()
        covid_variant.main(
            fq1=args.fq1,
            fq2=args.fq2,
            outdir=args.outdir,
            tolerate_missing=args.tolerate_missing,
            target_coverage=args.target_coverage,
            threads=args.threads,
            debug=args.debug)


if __name__ == '__main__':
    EntryPoint().main()
