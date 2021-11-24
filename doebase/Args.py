from os import getcwd, path
from argparse  import ArgumentParser
from ._version import __version__

DEFAULT_libsize = 32
DEFAULT_get_sequences = True
DEFAULT_backtranslate = True
DEFAULT_condon_table = 'Eecoli.cut'
DEFAULT_parts_file = path.join(
    path.dirname(path.abspath(__file__)),
    'data',
    'ref_parts.csv'
    )

def build_args_parser(
    prog: str,
    description: str = '',
    epilog: str = ''
) -> ArgumentParser:

    parser = ArgumentParser(
        prog = prog,
        description = description,
        epilog = epilog
    )

    # Build Parser
    parser = add_arguments(parser)

    return parser

def add_arguments(parser: ArgumentParser) -> ArgumentParser:
    parser.add_argument(
        '--func',
        required=True,
        type=str,
        choices=[
            # 'makeDoeOptDes',
            # 'callDoE',
            # 'evaldes',
            # 'doeRequest',
            # 'getDoe',
            # 'mainDoe',
            'doeGetSBOL',
        ],
        help='function to process'
    )
    parser.add_argument(
        'genes_file',
        type=str,
        help='Genes File'
    )
    parser.add_argument(
        'outfile',
        type=str,
        help='Ouput (SBOL) file'
    )
    parser.add_argument(
        '--ref_parts_file',
        type=str,
        default=DEFAULT_parts_file,
        help=f'Parts File (defualt: {DEFAULT_parts_file})'
    )
    parser.add_argument(
        '--sbol_file',
        type=str,
        help='SBOL file containing optimised versions of the genes (RBS, etc)'
    )
    parser.add_argument(
        '--libsize',
        type=int,
        default=DEFAULT_libsize,
        help=f'(default: {DEFAULT_libsize})'
    )
    parser.add_argument(
        '--get_sequences',
        type=bool,
        default=DEFAULT_get_sequences,
        help=f'(default: {DEFAULT_get_sequences})'
    )
    parser.add_argument(
        '--backtranslate',
        type=bool,
        default=DEFAULT_backtranslate,
        help=f'(default: {DEFAULT_backtranslate})'
    )
    parser.add_argument(
        '--codon_table',
        type=str,
        default=DEFAULT_condon_table,
        help=f'(default: {DEFAULT_condon_table})'
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {}'.format(__version__),
        help='show the version number and exit'
    )
    return parser
