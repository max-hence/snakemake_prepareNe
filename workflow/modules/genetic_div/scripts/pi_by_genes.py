from pandas import read_csv
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import sys

def barre_chargement(progression, total, longueur=50):
    pourcentage = progression / total
    bar_length = int(longueur * pourcentage)
    bar = '█' * bar_length + '-' * (longueur - bar_length)
    sys.stdout.write(f"\r|{bar}| {round(pourcentage * 100, 3)}%")
    sys.stdout.flush()

def p0p4_by_genes(pi_by_sites:str, dgnracy_bed:str, dgnracy_summary:str, output_path:str):
    """Measure p0p4 ratios by genes

    Args:
        pi_by_sites (str): _description_
        gene_bed (str): _description_
        output_path (str): _description_
    """
    header_bed = ["seq_id", "start", "end", "rna_id", "folding", "anc", "aa", "alt_aa"]
    dgnracy_df = read_csv(dgnracy_bed, delimiter="\t", names=header_bed, index_col=None)

    header_summ = ["seq_id", "start", "end", "rna_id", "cds_length", "mrna_length", "is_longest", "f0",	"f2", "f3", "f4"]
    summ_df = read_csv(dgnracy_summary, delimiter="\t", names=header_summ, index_col=None)
    len_summ = summ_df.shape[0]


    with open(pi_by_sites, 'r') as pi:
        with open(output_path, "w") as p0p4_by_gene:
            summ_idx = 0
            sum_p0, sum_p4 = 0, 0
            barre_chargement(summ_idx, len_summ)
            gene_start, gene_end = summ_df.iloc[summ_idx]["start"], summ_df.iloc[summ_idx]["end"]
            for row in pi:
                chrom = row.split("\t")[0]
                if chrom == "CHROM": continue
                pos, pi = int(row.split("\t")[1]), float(row.strip().split("\t")[2])

                # if pos before gene start
                if pos < gene_start: continue

                # if pos after gene end
                while pos >= gene_end and summ_idx < len_summ:
                    if sum_p0 * sum_p4 != 0:
                        ratio = p0p4_ratio(sum_p0, sum_p4, int(summ_df.iloc[summ_idx]["f0"]), int(summ_df.iloc[summ_idx]["f4"]))
                        p0p4_by_gene.write(f"{chrom}\t{gene_start}\t{gene_end}\t{ratio}\n")
                        barre_chargement(summ_idx, len_summ)
                    sum_p0, sum_p4 = 0, 0
                    summ_idx += 1
                    gene_start, gene_end = summ_df.iloc[summ_idx]["start"], summ_df.iloc[summ_idx]["end"]

                # if pos inside genes
                if pos in dgnracy_df["start"].values:
                    folding = dgnracy_df[dgnracy_df["start"] == pos]["folding"].values[0]
                    if folding == 4: sum_p4 += pi
                    elif folding == 0: sum_p0 += pi


def pi_by_genes(pi_by_sites:str, gene_bed:str, output_path:str):
    """Measure pi by genes

    Args:
        pi_by_sites (str): _description_
        gene_bed (str): bed file with gene intervals
        output_path (str): _description_
    """

    genes_df = read_csv(gene_bed, delimiter="\t", names=["seq_id", "start", "end", "rna_id"], index_col=None)
    len_genes = genes_df.shape[0]
    with open(pi_by_sites, 'r') as pi:
        with open(output_path, "w") as pi_by_gene:
            gene_idx = 0
            sum_pi = 0
            barre_chargement(gene_idx, len_genes)
            gene_start, gene_end = int(genes_df.iloc[gene_idx]["start"]), int(genes_df.iloc[gene_idx]["end"])
            for row in pi:
                chrom = row.split("\t")[0]
                if chrom == "CHROM": continue
                pos, pi = int(row.split("\t")[1]), float(row.strip().split("\t")[2])

                # if after gene
                while pos >= gene_end and gene_idx != genes_df.shape[0]:
                    if sum_pi != 0: # if the vcf skipped a gene
                        mean_pi = round(sum_pi / (gene_end - gene_start), 10)
                        pi_by_gene.write(f"{chrom}\t{gene_start}\t{gene_end}\t{mean_pi}\n")
                        barre_chargement(gene_idx, len_genes)
                    
                    sum_pi = 0
                    gene_idx += 1
                    gene_start, gene_end = genes_df.iloc[gene_idx]["start"], genes_df.iloc[gene_idx]["end"]
                    
                
                # if before gene
                if pos < gene_start:  continue

                # if inside
                else: sum_pi += pi


def p0p4_ratio(sum_p0:int, sum_p4:int, n_f0:int, n_f4:int):
    """_summary_

    Args:
        sum_p0 (int): _description_
        sum_p4 (int): _description_
        n_f0 (int): _description_
        n_f4 (int): _description_
    """

    return (sum_p0/n_f0)/(sum_p4/n_f4)


def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input',type=str,
        help="Path to pi by sites .pi"
    )
    parser.add_argument('--p0p4', action='store_true',
        help="For p0p4 ratio"
    )
    parser.add_argument('--pi', action='store_true',
        help="For pi"
    )
    parser.add_argument('-d', '--degeneracy', type=str,
        help="Path to gene positions .bed"
    )
    parser.add_argument('-s', '--summary', type=str,
        help="Path to gene positions .bed"
    )
    parser.add_argument('-o', '--output',type=str,
        help="Path to bed"
    )
    
    ## if no args then return help message
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()

    return args


if __name__ == "__main__":
 
    args = parse_command_line()

    if args.p0p4:
        p0p4_by_genes(args.input, args.degeneracy, args.summary, args.output)

    elif args.pi:
        pi_by_genes(args.input, args.summary, args.output)