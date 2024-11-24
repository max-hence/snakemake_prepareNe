from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pandas import read_csv
import sys

def rescale_snps(snps_stats:str, raw_stats:str, input_fai:str, output_fai:str):
    """ Rescale chr length based on removed variants

    Args:
        snps_stats (str): stats file from bcftools index givin nbr of snps for one chr
        raw_stats (str): stats file of vcf before filtering
        input_fai (str): previous fai
        output_fai (str): new fai
    """
    with open(snps_stats, "r") as snps:
        chr_id_snps, chr_len_snps, variants_snps = snps.readline().strip().split("\t")
    
    with open(raw_stats, "r") as raw:
        chr_id_raw, chr_len_raw, variants_raw = raw.readline().strip().split("\t")

    if chr_id_snps != chr_id_raw: raise ValueError("AH CA VA PAS DU TOUT")

    with open(input_fai, "r") as fai:
        with open(output_fai, "w") as new_fai:
            for row in fai:
                chr_id, chr_len = row.strip().split('\t')
                if chr_id == chr_id_snps:
                    new_len = int((int(chr_len) * int(variants_snps)) / int(variants_raw))
                    new_fai.write(f"{chr_id}\t{new_len}\n")


def rescale_bed(input_bed:str, input_fai:str, output_fai:str):
    """Determine the correct chr length based on kept bed regions

    Args:
        input_bed (str): _description_
        input_fai (str): _description_
        output_fai (str): _description_
    """
    bed_table = read_csv(input_bed, header=None, delimiter="\t", index_col=None, names=["chrom", "chromStart", "chromEnd", "n_samples","samples"])
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
    parser.add_argument('-r', '--raw', type=str,
        help="Stats of vcf before filtering"    
    )
    parser.add_argument('-f', '--fai',type=str,
        help="Specifies the path to fai"
    )
    parser.add_argument('-o', '--output',type=str,
        help="Specifies the path to modified fai"
    )
    parser.add_argument('--method', type=str,
        help="Tells which rescaling method to apply (snp or bed)"
    )
    
    ## if no args then return help message
    if len(sys.argv) == 4:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    return args


if __name__== "__main__":

    args = parse_command_line()

    if args.method == "snp":
        rescale_snps(args.input, args.raw, args.fai, args.output)

    elif args.method == "bed":
        rescale_bed(args.input, args.fai, args.output)
    else:
        raise ValueError("method arg is either snp or bed")