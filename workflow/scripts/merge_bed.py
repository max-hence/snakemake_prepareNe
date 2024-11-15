from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pandas import read_csv
import sys

def merge_beds(bed_map:str, bed_call:str, output_bed:str):
    """Determine the correct chr length based on kept bed regions

    Args:
        input_bed (str): _description_
        input_fai (str): _description_
        output_fai (str): _description_
    """
    bed_map_df = read_csv(bed_map,  header=None, delimiter="\t", index_col=None, names=["chrom", "chromStart", "chromEnd"])
    map_idx = 0
    map_start = int(bed_map_df.iloc[map_idx]["chromStart"])
    map_end = int(bed_map_df.iloc[map_idx]["chromEnd"])

    with open(output_bed, "w") as output:
        with open(bed_call, "r") as bed:
            for row in bed:
                end = int(row.split("\t")[2])
                if end <= map_start:
                    continue
                else:
                    if end > map_end:
                        map_idx += 1
                        if map_idx == bed_map_df.shape[0]: break
                        map_start = int(bed_map_df.iloc[map_idx]["chromStart"])
                        map_end = int(bed_map_df.iloc[map_idx]["chromEnd"])
                    else:
                        output.write(row)


def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--bed-map',type=str,
        help="Specifies the path to map.bed"
    )
    parser.add_argument('-b', '--bed-call',type=str,
        help="Specifies the path to callable.bed"
    )
    parser.add_argument('-o', '--output',type=str,
        help="Specifies the path to merged beds"
    )
    
    ## if no args then return help message
    if len(sys.argv) == 4:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    return args


if __name__== "__main__":

    args = parse_command_line()
    merge_beds(args.bed_map, args.bed_call, args.output)