
def get_chr_list(fai_path:str):
    
    return [row.split('\t')[0] for row in open(fai_path, "r")]


def get_output():
    chromosomes: list = get_chr_list(config["fai_path"])

    if config["final_prefix"] == "":
            raise (WorkflowError("'final_prefix' is not set in config."))

    # vcf filtered to keep only bi-allelic SNPS
    out.extend(expand("results/geneD/{prefix}.geneD", prefix=config["final_prefix"], chr=chromosomes))
    
    # splitted vcf
    out.extend(expand("results/geneD/ne_inference/strway_plot/{prefix}.lowD.{chr}.vcf.gz", prefix=config["final_prefix"], chr=chromosomes))
    out.extend(expand("results/geneD/ne_inference/strway_plot/{prefix}.highD.{chr}.vcf.gz.tbi", prefix=config["final_prefix"], chr=chromosomes))

#    out.extend(expand("results/geneD/ne_inference/smcpp/{prefix}.geneD.{chr}.smc.gz", prefix=config["final_prefix"], chr=chromosomes))
#    out.extend(expand("results/geneD/ne_inference/strway_plot/{prefix}.geneD.{chr}.blueprint", prefix=config["final_prefix"], chr=chromosomes))