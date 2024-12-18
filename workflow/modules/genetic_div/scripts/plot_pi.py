from pandas import read_csv
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def plot_pi(infiles:list, window:bool, plot_path:str):
    """
    Args:
        infiles (str): list of .pi files
        plot_path (str): _description_
    """

    fig, ax = plt.subplots(figsize=(10, 4))
    color = ["#DE464F", "#69BBAD"]
    for i, input_pi in enumerate(infiles):
        if window:
            df_pi = read_csv(input_pi, sep='\t') #, header=None, names=["CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "PI"])
            # Central position for each gene
            df_pi['center'] = (df_pi['BIN_START'] + df_pi['BIN_END']) / 2
            ax.plot(df_pi['center'], df_pi['PI'], color=color[i % len(color)], label=input_pi.split("/")[-1]) #, s=5) # scaterplot for window < 50000
            ax.set_xlim(df_pi['BIN_START'][2], df_pi['BIN_END'][len(df_pi)-3])
        else:
            df_pi = read_csv(input_pi, sep='\t')#, header=None, names=['CHR', 'POS', 'PI'])
            ax.scatter(df_pi['POS'], df_pi['PI'], color=color[i % len(color)], label="pi", s=5)

    # Configurations supplémentaires
    ax.set_xlabel(f"{df_pi['CHROM'][0]}")
    ax.set_ylabel("pi")
    ax.legend()
    plt.savefig(plot_path)


def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input', nargs="+", required=True,
        help="List of pi files"
    )
    parser.add_argument('--window', action="store_true",
        help="if format is 'chrom', 'start', 'end', 'pi' instead of 'chrom', 'pos', 'pi'"
    )
    parser.add_argument('-o', '--output',type=str,
        help="Path to png"
    )
    
    ## if no args then return help message
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_command_line()
    plot_pi(args.input, args.window, args.output)