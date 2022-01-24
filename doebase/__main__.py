from .Args import build_args_parser
from .synbioParts import doeGetSBOL


def entry_point():
    parser = build_args_parser(
        prog = 'doebase',
        description='An optimal design of experiments (DoE) base package for synthetic biology'
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
