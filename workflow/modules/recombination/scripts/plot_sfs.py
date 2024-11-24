from argparse import ArgumentParser, RawDescriptionHelpFormatter
import matplotlib.pyplot as plt
import numpy as np
import sys


def plot_sfs(sfs_path:str, output_path:str):
    """Plot sfs

    Args:
        sfs_path (str): _description_
        output_path (str): _description_
    """
    with open(sfs_path, 'r') as file:
        sfs = np.array([int(value) for value in file.read().split(' ')])

    table = np.column_stack((np.arange(1, len(sfs) + 1), sfs))

    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.bar(table[:,0], table[:,1], align='center', alpha=0.7)
    plt.title('SFS')
    plt.xlabel('Frequency in Population')
    plt.ylabel('Proportion of SNPs')

    plt.savefig(output_path)

def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input',type=str,
        help="Specifies the path to sfs"
    )
    parser.add_argument('-o', '--output',type=str,
        help="Specifies the path to plot"
    )
    
    ## if no args then return help message
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    return args


if __name__=="__main__":

    args = parse_command_line()
    
    plot_sfs(args.input, args.output)