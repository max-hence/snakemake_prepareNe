from argparse import ArgumentParser, RawDescriptionHelpFormatter

def snp2bed(infile, outfile):
    """Make a bed file from a list of positions

    Args:
        infile (_type_): _description_
        outfile (_type_): _description_
    """
    first_row=True
    with open(outfile, "w") as bed: 
        with open(infile, "r") as snps:
            for row in snps:
                chr_id, new_start = row.strip().split("\t")
                if first_row: 
                    start, end = new_start, int(new_start) + 1
                    first_row=False
                else:
                    if int(new_start) == end: end += 1
                    else:
                        bed.write(f"{chr_id}\t{start}\t{end}\n")
                        start, end = new_start, int(new_start) + 1
        bed.write(f"{chr_id}\t{start}\t{end}\n")    

def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input',type=str,
        help=""
    )
    parser.add_argument('-o', '--output',type=str,
        help=""
    )

    args = parser.parse_args()

    return args

if __name__=="__main__":

    args = parse_command_line()

    snp2bed(args.input, args.output)
