from logging import Logger
from .Args import build_args_parser
from .synbioParts import doeGetSBOL


def entry_point():
    parser = build_args_parser(
        prog = 'partsgenie_client',
        description='Requests a input (SBOL) file to a PartsGenie server'
    )
    args = parser.parse_args()

    if args.func == 'doeGetSBOL':
        diagnostics = doeGetSBOL(
            pfile=args.ref_parts_file,
            gfile=args.genes_file, 
            gsbol=args.sbol_file,
            libsize=args.libsize,
            getSequences=args.get_sequences,
            backtranslate=args.backtranslate,
            codontable=args.codon_table
        )
        open(args.outfile, 'w').write(diagnostics['sbol'])


if __name__ == '__main__':
    entry_point()
