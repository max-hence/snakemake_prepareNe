def get_chr_list(fai_path:str):
    
    return [row.split('\t')[0] for row in open(fai_path, "r")]


def get_output():
    chromosomes: list = get_chr_list(config["fai_path"])
    rec_quant = ["rec1","rec2","rec3"]
    final_prefix = config["final_prefix"]
    out: list = []

    if final_prefix == "":
            raise (WorkflowError("'final_prefix' is not set in config."))

    out.extend(expand("results/rec/vcf/{prefix}.{rec}.{chr}.vcf", prefix=final_prefix, rec=rec_quant, chr=chromosomes))
    out.extend(expand("results/rec/bed/{prefix}.{rec}.{chr}.callable.bed", prefix=final_prefix, rec=rec_quant, chr=chromosomes))    
    out.extend(expand("results/rec/stats/{prefix}.{rec}.rdmSNP.{chr}.fai", prefix=final_prefix, rec=rec_quant, chr=chromosomes))
    out.extend(expand("results/rec/sfs/{rec}/{prefix}.{rec}.{chr}.png", prefix=final_prefix, rec=rec_quant, chr=chromosomes))
    out.extend(expand("results/rec/ne_inference/smcpp/{rec}/{prefix}.{rec}.{chr}.smc.gz", prefix=final_prefix, rec=rec_quant, chr=chromosomes))
    out.extend(expand("results/rec/ne_inference/strway_plt/{rec}/{prefix}.{rec}.{chr}.blueprint", prefix=final_prefix, rec=rec_quant, chr=chromosomes))

    return out


def get_previews(preview_path:str):
    chromosomes: list = get_chr_list(config["fai_path"])
    rec_quant = ["rec1","rec2","rec3"]
    return expand(preview_path, prefix=config["final_prefix"], rec=rec_quant, chr=chromosomes)