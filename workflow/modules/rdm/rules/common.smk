
def get_chr_list(fai_path:str):
    """Return the list of scf that will analyzed

    Args:
        fai_path (str): tab separated table (fai format). 
                        Genome info, with only scf that will be analyzed
    Returns:
        list of scaffolds of interest
    """
    
    return [row.split('\t')[0] for row in open(fai_path, "r")]


def get_output():
    chromosomes: list = get_chr_list(config["fai_path"])
    rec: list = ["rec1","rec2","rec3"]
    out: list  = []

    if config["final_prefix"] == "":
            raise (WorkflowError("'final_prefix' is not set in config."))

    out.extend(expand("results/sfs/rdm/{prefix}.rdm.{chr}.png", prefix=config["final_prefix"], chr=chromosomes))
    out.extend(expand("results/ne_inference/strway_plt/rdm/{prefix}.rdm.{chr}.blueprint", prefix=config["final_prefix"], chr=chromosomes))
    out.extend(expand("results/ne_inference/smcpp/rdm/{prefix}.rdm.{chr}.smc.gz", prefix=config["final_prefix"], chr=chromosomes))
    
    return out


def get_previews(vcf_path:str):
    chromosomes = get_chr_list(config["fai_path"])
    return expand(vcf_path, prefix=config["final_prefix"], chr=chromosomes)