from argparse import ArgumentParser, RawDescriptionHelpFormatter
import sys

def reverse_bed(bed:str, last_pos:str, mask:str):
    """Reverse bed file to show regions to exclude
    e.g. 
        BED : 
            1	0	30000	3	msp_0,msp_1,msp_2
            1	50000	1000000	5	msp_0,msp_1,msp_2,msp_3,msp_4
        MASK:
            1	30000	50000

    Args:
        bed (str): _description_
        last_pos: last chrEnd (can be absent if last bed's lines have been removed)
        mask (str): _description_
    """
    start = "0"
    with open(mask, 'w') as reverse_bed:
        with open(bed, 'r') as file:
            for row in file:
                chrom, new_start, new_end = row.split("\t")[:3]
                #print('Avant : ', new_end)
                if new_start == start:
                    start = new_end
                    continue
                else:
                    reverse_bed.write(f"{chrom}\t{start}\t{new_start}\n")
                    start = new_end
            #print("Apres : ", new_end)
            if int(new_end) < int(last_pos):
                reverse_bed.write(f"{chrom}\t{new_end}\t{last_pos}\n")


def parse_command_line():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""\n
    """)
    parser.add_argument('-i', '--input',type=str,
        help="Specifies the path to bed file"
    )
    parser.add_argument('-l', '--last-pos',type=str,
        help="Last chrEnd of raw bed"
    )
    parser.add_argument('-o', '--output',type=str,
        help="Specifies the path to mask"
    )
    
    ## if no args then return help message
    if len(sys.argv) == 3:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    return args

if __name__=="__main__":

    args = parse_command_line()
    reverse_bed(args.input, args.last_pos, args.output)