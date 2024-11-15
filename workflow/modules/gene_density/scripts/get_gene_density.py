import gffpandas.gffpandas as gffpd
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import sys

def get_gene_density(gff_path:str, window_size:int, tbl_path:str):
    """ Measure gene density by window

    Args:
        gff_path (str): path to .gff
        tbl_path (str): path and prfx to density table
    """
    window_size = int(window_size)
    annotation_gff = gffpd.read_gff3(gff_path)
    genes_gff = annotation_gff.filter_feature_of_type(['gene']).df
    print(genes_gff)
    genes_gff.loc[:, 'start'] = genes_gff.loc[:, 'start'] // window_size
    gene_density = genes_gff.groupby('start').size().reset_index(name='density')
    gene_density["start"] = gene_density["start"] * window_size
    gene_density.insert(1, "end", gene_density["start"] + window_size)

    gene_density.to_csv(tbl_path, sep="\t", header=False, index=False)


def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input', type=str,
        help="Path to gff"
    )
    parser.add_argument('-w', '--window', type=str,
        help="Path to gff"
    )
    parser.add_argument('-o', '--output', type=str,
        help="Path to density table"
    )


    ## if no args then return help message
    if len(sys.argv) == 2:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()

    return args


if __name__ == "__main__":
 
    args = parse_command_line()

    get_gene_density(args.input, args.window, args.output)