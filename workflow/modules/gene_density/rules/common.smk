
def get_chr_list(fai_path:str):
    
    return [row.split('\t')[0] for row in open(fai_path, "r")]

def get_output():
    final_prefix = config["final_prefix"]
    
    if final_prefix == "":
        raise (WorkflowError("'final_prefix' is not set in config."))
    
    out: list  = [] # output list
    chromosomes: list = get_chr_list(config["fai_path"])
    density_cat = ["lowD","highD"]

    out.extend(expand("results/geneD/stats/{prefix}.{chr}.geneD.png", prefix=final_prefix, chr=chromosomes))

    # Stairway Plot inputs
    out.extend(expand("results/geneD/ne_inference/strway_plt/{density}/{prefix}.{density}.{chr}.blueprint", prefix=final_prefix, density=density_cat, chr=chromosomes))
    # SMC++ inputs
    out.extend(expand("results/geneD/ne_inference/smcpp/{density}/{prefix}.{density}.{chr}.smc.gz", prefix=final_prefix, density=density_cat, chr=chromosomes))

    return out

def get_previews(wc):

    return expand("results/geneD/sfs/{{prefix}}.{{density}}.rdmSNP.preview.{chr}.txt", chr=get_chr_list(config["fai_path"]))