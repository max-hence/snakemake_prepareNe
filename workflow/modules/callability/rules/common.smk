
def get_chr_list(fai_path:str):
    
    return [row.split('\t')[0] for row in open(fai_path, "r")]

def get_previews(wc):

    return expand("results/callability/sfs/{{prefix}}.SNPS.NA.preview.{{chr}}.txt")