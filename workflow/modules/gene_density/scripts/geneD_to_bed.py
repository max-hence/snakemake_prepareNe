from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pandas import read_csv
import sys

def geneD2bed(geneD_path:str, output_bed:str, low_range:int=0, up_range:int=100, plot:bool=False):
    """Writes a bed file that includes position within a given gene density interval
    run chr by chr
    Args:
        geneD_path (str): path/to/file.rec (with space divided columns)
        bed_path (str): output_path/to/file.bed
        low_range (int, optional): Defaults to 0.
        high_range (int, optional): Defaults to 100.
        ex : low_range = 50, high_range = 100 will keep only the regions with the 50% highest rec rate
        plot (bool, optional): if you want to plot (.png or .svg)
    """

    geneD_df = read_csv(geneD_path, delimiter = "\t", header = None, names = ["seq_id", "start", "end", "density"], index_col=None)
    quantiles = [int(q)/100 for q in [low_range, up_range]]

    # defines treshold by chromosome

    density_quantiles = geneD_df["density"].quantile(q=quantiles)
    density_quantiles = [density_quantiles[q] for q in quantiles] # dict of low and upper range for each chrom

    # writes bed file
    with open(output_bed, "w") as bedfile:
        start_intervall = None
        for _, row in geneD_df.iterrows():
            new_start, new_end, density = int(round(row["start"])), int(row["end"]), int(row["density"])
            # if position in rec rate range and not NA, interval is written
            if density >= density_quantiles[0] and density < density_quantiles[1]:
                if start_intervall == None: # to merge contiguous intervals
                    start_intervall = new_start
                previous_end = new_end
            else:
                if start_intervall != None:
                    bedfile.write(f"{row["seq_id"]}\t{start_intervall}\t{previous_end}\n")
                start_intervall = None

        if start_intervall != None: bedfile.write(f"{row["seq_id"]}\t{start_intervall}\t{new_end}\n")


def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input', type=str,
        help="Path to file.gendist"
    )
    parser.add_argument('-o', '--output', type=str,
        help="Path to file.bed"
    )
    parser.add_argument('-l', '--low', type=int, default=0,
        help="lower range"
    )
    parser.add_argument('-u', '--up', type=int, default=100,
        help="upper range"
    )
    
    ## if no args then return help message
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if not (args.low or args.up):
        parser.error("At least one of -l (low) or -u (up) must be provided.")

    if args.low < 0 or args.low > 100 or args.up < 0 or args.low > 100:
        parser.error("low_range and up_range must be pourcentages btw 0 and 100")
    return args


if __name__=="__main__":

    args = parse_command_line()
    geneD2bed(args.input, args.output, args.low, args.up)