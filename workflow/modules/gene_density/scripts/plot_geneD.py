from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pandas import read_csv
import matplotlib.pyplot as plt
import sys

def cut_chr_by_geneD(geneD_df):
    first_idx = geneD_df.index[0]
    start = first_idx
    color = geneD_df["color"][start]
    chr_chunks = []
    for i, row in geneD_df.iterrows():
        if row["color"] != color:
            chunk = geneD_df.iloc[start-first_idx:i-first_idx, :]
            chr_chunks.append(chunk)
            start = i
            color = row["color"]
    if start < i:
        chr_chunks.append(geneD_df.iloc[start-first_idx:i-first_idx, :])
    return chr_chunks

def plot_geneD(geneD_path:str, thresholds:list, plot_path:str):

    geneD_df = read_csv(geneD_path, delimiter = "\t", header = None, names = ["seq_id", "start", "end", "density"], index_col=None)
    
    # defines tresholds by chromosome
    quantiles = [int(q)/100 for q in thresholds]
    density_quantiles = geneD_df["density"].quantile(q=quantiles)
    density_quantiles = [density_quantiles[q] for q in quantiles]

    if len(density_quantiles) == 1: density_quantiles.append(density_quantiles[0])

    geneD_df["color"] = geneD_df['density'].apply(lambda x: '#D44B53' if x >= density_quantiles[1] else '#009e73' if x < density_quantiles[0] else '#7F7F7F') # red, blue, grey
    
    plt.figure(figsize=(14, 8))
    plt.plot(geneD_df['start'], geneD_df['density'], linestyle='-', color="#7F7F7F", lw=2)

    chr_chunks = cut_chr_by_geneD(geneD_df)
    for chunk in chr_chunks:
        plt.plot(chunk['start'], chunk['density'], color=chunk['color'].unique()[0], lw=2) # degue mais ça passe
    
    plt.savefig(plot_path)
    # need to divide density map by chunk of same density category for goof visualisation
    

def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input', type=str,
        help="Path to file.gendist"
    )
    parser.add_argument('-t', '--threshold', nargs="+", type=int,
        help="One threshold or two to split density map in two, in %"
        # either one value and one color for values below and other for values above threshold
        # or two threshold and one color for values below first and other of values above 2nd
    )
    parser.add_argument('-o', '--output', type=str,
        help="Path to file.bed"
    )

    args = parser.parse_args()

    args.threshold = [int(n) for n in args.threshold]
    if len(args.threshold) > 2:
        raise ValueError("""--threshold arg must be one or two values (in %)\n
        - Either one threshold then one color for values below and another for values above (e.g. 50)\n
        - Or two threshold and one color for values below first and other of values above 2nd (e.g. 40 60)
    """)
    return args


if __name__=="__main__":

    args = parse_command_line()
    plot_geneD(args.input, args.threshold, args.output)