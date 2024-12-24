
def get_chr_list(fai_path:str):
    
    return [row.split('\t')[0] for row in open(fai_path, "r")]

def get_previews(wc):

    return expand("results/geneD/sfs/{{prefix}}.{density}.rdmSNP.preview.{chr}.txt", chr=get_chr_list(config["fai_path"]), density=["lowD", "highD"])