
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

    # vcf filtered to keep only bi-allelic SNPS
    out.extend(expand("results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz", prefix=config["final_prefix"], chr=chromosomes))
    out.extend(expand("results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi", prefix=config["final_prefix"], chr=chromosomes))

    # subsampled vcf (after correction on callability)
    out.extend(expand("results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz", prefix=config["final_prefix"], chr=chromosomes))
    out.extend(expand("results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz.tbi", prefix=config["final_prefix"], chr=chromosomes))
    
    if not ["stop"]:       
        # SFS plots
        out.extend(expand("results/sfs/snps/small/{prefix}.small.{chr}.png", prefix=config["final_prefix"], chr=chromosomes))
        out.extend(expand("results/sfs/snps/max/{prefix}.max.{chr}.png", prefix=config["final_prefix"], chr=chromosomes))
        out.extend(expand("results/sfs/snps/strict/{prefix}.strict.{chr}.png", prefix=config["final_prefix"], chr=chromosomes))
        out.extend(expand("results/sfs/snps_na/{prefix}.{chr}.png", prefix=config["final_prefix"], chr=chromosomes))
        
        # Stairway Plot inputs
        out.extend(expand("results/ne_inference/strway_plt/snps/small/{prefix}.SNPS.small.{chr}.blueprint", prefix=config["final_prefix"], chr=chromosomes))
        out.extend(expand("results/ne_inference/strway_plt/snps/max/{prefix}.SNPS.max.{chr}.blueprint", prefix=config["final_prefix"], chr=chromosomes))
        out.extend(expand("results/ne_inference/strway_plt/snps/ml/{prefix}.SNPS.ml.{chr}.blueprint", prefix=config["final_prefix"], chr=chromosomes))
        out.extend(expand("results/ne_inference/strway_plt/snps/strict/{prefix}.SNPS.strict.{chr}.blueprint", prefix=config["final_prefix"], chr=chromosomes))
        out.extend(expand("results/ne_inference/strway_plt/snps_na/{prefix}.SNPS.NA.{chr}.blueprint", prefix=config["final_prefix"], chr=chromosomes))

        # SMC++ inputs
        out.extend(expand("results/ne_inference/smcpp/snps/full/{prefix}.SNPS.full.{chr}.smc.gz", prefix=config["final_prefix"], chr=chromosomes))
        out.extend(expand("results/ne_inference/smcpp/snps/small/{prefix}.SNPS.small.{chr}.smc.gz", prefix=config["final_prefix"], chr=chromosomes))
        out.extend(expand("results/ne_inference/smcpp/snps/max/{prefix}.SNPS.max.{chr}.smc.gz", prefix=config["final_prefix"], chr=chromosomes))
        out.extend(expand("results/ne_inference/smcpp/snps/ml/{prefix}.SNPS.ml.{chr}.smc.gz", prefix=config["final_prefix"], chr=chromosomes))
        out.extend(expand("results/ne_inference/smcpp/snps/strict/{prefix}.SNPS.strict.{chr}.smc.gz", prefix=config["final_prefix"], chr=chromosomes))
        out.extend(expand("results/ne_inference/smcpp/snps_na/{prefix}.SNPS.NA.{chr}.smc.gz", prefix=config["final_prefix"], chr=chromosomes))


    # Recombination
    if config["recombination"]:
        out.extend(expand("results/bed/rec/max/{prefix}.{rec_quant}.{chr}.callable.bed", prefix=config["final_prefix"], rec_quant=rec, chr=chromosomes))
        out.extend(expand("results/stats/rec/ml/{prefix}.{rec_quant}.320kSNP.{chr}.fai", prefix=config["final_prefix"], rec_quant=rec, chr=chromosomes))
        out.extend(expand("results/sfs/rec/{rec_quant}/{prefix}.{rec_quant}.{chr}.png", prefix=config["final_prefix"], rec_quant=rec, chr=chromosomes))
        out.extend(expand("results/ne_inference/smcpp/rec/{rec_quant}/{prefix}.{rec_quant}.{chr}.smc.gz", prefix=config["final_prefix"], rec_quant=rec, chr=chromosomes))
        out.extend(expand("results/ne_inference/strway_plt/rec/{rec_quant}/{prefix}.{rec_quant}.{chr}.blueprint", prefix=config["final_prefix"], rec_quant=rec, chr=chromosomes))
    
    return out


def get_previews(vcf_path:str):
    chromosomes = get_chr_list(config["fai_path"])
    return expand(vcf_path, prefix=config["final_prefix"], chr=chromosomes)