
def check_fai_format(fai_path:str):
    with open(fai_path, "r") as fai:
        if len(fai_path.readline().split("\t")) != 2:
            raise (WorkflowError(" fai must contain only two tab-separated columns : chr_id\tchr_length"))


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
    final_prefix = config["final_prefix"]
    chromosomes: list = get_chr_list(config["fai_path"])
    rec: list = ["rec1","rec2","rec3"]
    sfs_params_method: list = ["ml", "max"]
    sampling_method: list = ["small","max","ml","strict"]
    out: list  = []

    if final_prefix == "":
        raise (WorkflowError("'final_prefix' is not set in config."))

    # vcf filtered to keep only bi-allelic SNPS
    out.extend(expand("results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz", prefix=final_prefix, chr=chromosomes))
    out.extend(expand("results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi", prefix=final_prefix, chr=chromosomes))

    if not config["stop"]:      
        # SFS plots
        # out.extend(expand("results/sfs/snps/{subsample}/{prefix}.{subsample}.{chr}.png", prefix=final_prefix, chr=chromosomes, subsample=sampling_method))
        # out.extend(expand("results/sfs/snps_na/{prefix}.{chr}.png", prefix=final_prefix, chr=chromosomes))
        
        # Stairway Plot inputs
        out.extend(expand("results/ne_inference/strway_plt/{subsample}/{prefix}.SNPS.{subsample}.{chr}.blueprint", prefix=final_prefix, chr=chromosomes, subsample=sampling_method))
        out.extend(expand("results/ne_inference/strway_plt/snps_na/{prefix}.SNPS.NA.{chr}.blueprint", prefix=final_prefix, chr=chromosomes))

        # SMC++ inputs
        out.extend(expand("results/ne_inference/smcpp/{subsample}/{prefix}.SNPS.{subsample}.{chr}.smc.gz", prefix=final_prefix, chr=chromosomes, subsample=sampling_method))
        out.extend(expand("results/ne_inference/smcpp/full/{prefix}.SNPS.full.{chr}.smc.gz", prefix=final_prefix, chr=chromosomes))
        out.extend(expand("results/ne_inference/smcpp/snps_na/{prefix}.SNPS.NA.{chr}.smc.gz", prefix=final_prefix, chr=chromosomes))

    if config["density"]:
        out.extend(rules.density_all.input)
    if config["recombination"]:
        out.extend(rules.rec_all.input)
    # if config["genetic_div"]: out.append(rules.genetic_div_all.input)
    return out


def get_previews(preview_path:str):
    chromosomes = get_chr_list(config["fai_path"])
    return expand(preview_path, prefix=config["final_prefix"], chr=chromosomes)