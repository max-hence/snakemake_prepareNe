# Plot gene density of a chromosome colored by categories

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
    """
        geneD_path:
        thresholds:
        plot_path:
    """
    geneD_df = read_csv(geneD_path, delimiter = "\t", header = None, names = ["seq_id", "start", "end", "density"], index_col=None)
    
    # defines tresholds by chromosome
    quantiles = [int(q)/100 for q in thresholds]
    density_quantiles = geneD_df["density"].quantile(q=quantiles)
    density_quantiles = [density_quantiles[q] for q in quantiles]
    
    print(density_quantiles)
    if len(density_quantiles) == 1:
        geneD_df["color"] = geneD_df['density'].apply(lambda x: '#009e73' if x < density_quantiles[0] else '#D44B53') # red, green, grey
    else:
        geneD_df["color"] = geneD_df['density'].apply(
            lambda x: "#E69F00" if density_quantiles[0] < x <= density_quantiles[1] 
            else '#009E73' if x <= density_quantiles[0] 
            else '#D44B53'
        )
    plt.figure(figsize=(14, 8))
    plt.scatter(geneD_df['start'], geneD_df['density'], color=geneD_df["color"])

    plt.savefig(plot_path)

    # for curves but not good and harsh
    # need to divide density map by chunk of same density category for goof visualisation
    # chr_chunks = cut_chr_by_geneD(geneD_df)
    # for chunk in chr_chunks:
    #   plt.plot(chunk['start'], chunk['density'], color=chunk['color'].unique()[0], lw=2) # degue mais ça passe

def plot_recmap_by_geneD(geneD_path:str, thresholds:list, recmap_path:str, plot_path:str):
    """
        Color rec map by gene density
    """

    recmap_df = read_csv(recmap_path, delimiter = " ", header = None, names = ["rec_id", "seq_id", "start", "end", "rec_rate", "conf_intervall_low", "conf_intervall_up"], index_col=None)

    geneD_df = read_csv(geneD_path, delimiter = "\t", header = None, names = ["seq_id", "start", "end", "density"], index_col=None)
    
    # defines tresholds by chromosome
    quantiles = [int(q)/100 for q in thresholds]
    density_quantiles = geneD_df["density"].quantile(q=quantiles)
    density_quantiles = [density_quantiles[q] for q in quantiles]

    if len(density_quantiles) == 1:
        geneD_df["color"] = geneD_df['density'].apply(lambda x: '#009E73' if x < density_quantiles[0] else '#D44B53') # red, green, grey
    else:
        geneD_df["color"] = geneD_df['density'].apply(
            lambda x: "#E69F00" if density_quantiles[0] < x <= density_quantiles[1] 
            else '#009E73' if x <= density_quantiles[0] 
            else '#D44B53'
        )

    geneD_df['rec_rate'] = geneD_df.apply(lambda row: find_rec_rate(row, recmap_df), axis=1)
    geneD_df = geneD_df.dropna()

    plt.figure(figsize=(14, 8))
    plt.scatter(geneD_df["start"], geneD_df["rec_rate"], color=geneD_df["color"], s=25,)
    plt.savefig(plot_path)

def find_rec_rate(row, recmap_df):
    match = recmap_df[
        (recmap_df['seq_id'] == row['seq_id']) & 
        (recmap_df['start'] <= row['end']) & 
        (recmap_df['end'] >= row['start'])
    ]
    return match['rec_rate'].iloc[0] if not match.empty else None


def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input', type=str,
        help="Path to file.geneD"
    )
    parser.add_argument('-t', '--threshold', nargs="+", type=int,
        help="One threshold or two to split density map in two, in %"
        # either one value and one color for values below and other for values above threshold
        # or two threshold and one color for values below first and other of values above 2nd
    )
    parser.add_argument('-r', '--recmap', type=str,
        help="Path to recmap.rec"
    )
    parser.add_argument('-o', '--output', type=str,
        help="Path to file.png"
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
    plot_recmap_by_geneD(args.input, args.threshold, args.recmap, args.output)