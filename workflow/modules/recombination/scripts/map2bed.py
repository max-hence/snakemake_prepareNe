from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pandas import read_csv
import sys
from math import isnan


def map2bed(recMap_path:str, output_path:str, low_range:int=0, up_range:int=100):
    """Writes a bed file that includes position with a given recombination rate

    Args:
        recMap_path (str): path/to/file.rec (with space divided columns)
        output_path (str): output_path/to/file.bed
        low_range (int, optional): Defaults to 0.
        high_range (int, optional): Defaults to 100.
        ex : low_range = 50, high_range = 100 will keep only the regions with the 50% highest rec rate
    """

    rec_map = read_csv(recMap_path, delimiter = " ", header = 0, names = ["set", "map", "start", "end", "recRate", "upperRecRate", "lowerRecRate"])
    quantiles = [int(q)/100 for q in [low_range, up_range]]

    # defines treshold by chromosome

    rec_quantiles = rec_map["recRate"].quantile(q=quantiles)
    rec_quantiles = [rec_quantiles[q] for q in quantiles] # dict of low and upper range for each chrom

    # writes bed file
    with open(output_path, "w") as bedfile:
        # bedfile.write("chrom\tchromStart\tchromEnd\n") # header
        start_intervall = None
        for _, row in rec_map.iterrows():
            new_start, new_end, recRate = int(round(row["start"])), int(row["end"]), row["recRate"]

            # if position in rec rate range and not NA, interval is written
            if recRate >= rec_quantiles[0] and recRate <= rec_quantiles[1] and not isnan(recRate):
                if start_intervall == None: # to merge contiguous intervals
                    start_intervall = new_start
                previous_end = new_end
            else:
                if start_intervall != None:
                    bedfile.write(f"{row['map']}\t{start_intervall}\t{previous_end}\n")
                start_intervall = None

        if start_intervall != None: bedfile.write(f"{row['map']}\t{start_intervall}\t{new_end}\n")

    # print(f"{row['map']} : {round(rec_quantiles[0], 12)}, {round(rec_quantiles[1], 12)}")


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
    map2bed(args.input, args.output, args.low, args.up)