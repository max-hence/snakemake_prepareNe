def get_chr_list(fai_path:str):
    
    return [row.split('\t')[0] for row in open(fai_path, "r")]


def get_previews(preview_path:str):
    chromosomes: list = get_chr_list(config["fai_path"])
    rec_quant = ["rec1","rec2","rec3"]
    return expand(preview_path, prefix=config["final_prefix"], rec=rec_quant, chr=chromosomes)