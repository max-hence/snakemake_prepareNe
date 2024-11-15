from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pandas import read_csv
import sys

def rescale_map(input_bed:str, input_fai:str, output_fai:str):
    """Determine the correct chr length based on kept bed regions

    Args:
        input_bed (str): _description_
        input_fai (str): _description_
        output_fai (str): _description_
    """
    bed_table = read_csv(input_bed, header=None, delimiter="\t", index_col=None, names=["chrom", "chromStart", "chromEnd"])
    bed_table["difference"] = bed_table["chromEnd"] - bed_table["chromStart"]
    correct_size = bed_table['difference'].sum()

    with open(output_fai, "w") as output:
        with open(input_fai, "r") as fai:
            for line in fai:
                chr_id = line.split("\t")[0]
                end_line = '\t'.join(line.split('\t')[2:])

                output.write(f"{chr_id}\t{correct_size}\t{end_line}")


def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input',type=str,
        help="Specifies the path to vcf stats (if --snp) or bed_file (if --bed)"
    )
    parser.add_argument('-f', '--fai',type=str,
        help="Specifies the path to fai"
    )
    parser.add_argument('-o', '--output',type=str,
        help="Specifies the path to modified fai"
    )
    
    ## if no args then return help message
    if len(sys.argv) == 4:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    return args


if __name__== "__main__":

    args = parse_command_line()
    rescale_map(args.input, args.fai, args.output)