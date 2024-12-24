def get_chr_list(fai_path:str):
    
    return [row.split('\t')[0] for row in open(fai_path, "r")]

def get_previews(preview_path:str):
    # Gather all previews files from easySFS projects
    return expand("results/rec/sfs/{{prefix}}.{rec}.preview.{chr}.txt", chr=get_chr_list(config["fai_path"]), rec=["rec1","rec2","rec3"])