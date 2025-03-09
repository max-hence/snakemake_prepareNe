# Prepare inputs for MSMC2

mu = config["mutation_rate"]
generation = config["generation_time"]
time_segments = "4*1+25*1+4*1+6*1"
chromosomes: list = get_chr_list(config["fai_path"])

rule split_sample:
    """ Split vcf by samples and chromosomes and remove Nas """
    input:
        vcf = "results/snps/vcf/{prefix}.SNPS.{chr}.vcf.gz"
    output:
        vcf_by_sample = "results/snps/vcf/{prefix}.SNPS.{sample}.{chr}.vcf.gz",
        vcf_by_sample_idx = "results/snps/vcf/{prefix}.SNPS.{sample}.{chr}.vcf.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.{sample}.{chr}.log"
    params:
        sample="{wildcards.sample}"
    shell:
        """
        bcftools view -s {params.sample} {input.vcf} -Oz -o {output.vcf_by_sample}
        tabix -p vcf {output.vcf_by_sample}
        """

rule remove_na:
    """ Remove all NA from vcf and keep positions in a tsv file"""
    input:
        vcf_by_sample = "results/snps/vcf/{prefix}.SNPS.{sample}.{chr}.vcf.gz",
        vcf_by_sample_idx = "results/snps/vcf/{prefix}.SNPS.{sample}.{chr}.vcf.gz.tbi"
    output:
        vcf_no_na = "results/snps/vcf/{prefix}.SNPS.{sample}.{chr}.no_na.vcf.gz",
        vcf_no_na_idx = "results/snps/vcf/{prefix}.SNPS.{sample}.{chr}.no_na.vcf.gz.tbi",
        na_pos = "results/snps/bed/{prefix}.{sample}.{chr}.no_na.tsv"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.{sample}.{chr}.log"
    shell:
        """
        bcftools view -H -i 'GT~"\."' {input.vcf_by_sample} | \
            cut -f1,2 > {output.na_pos}

        bcftools view -e 'GT~"\."' -Oz -o {output.vcf_no_na} {input.vcf_by_sample}
        tabix -p vcf {input.vcf_by_sample}
        """
        
rule na2bed:
    """ Convert na positions as bed file"""
    input:
        na_pos = "results/snps/bed/{prefix}.{sample}.{chr}.no_na.tsv",
        script = workflow.source_path("../scripts/snp2bed.py")
    output:
        na_bed = "results/snps/bed/{prefix}.{sample}.{chr}.no_na.bed"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.{sample}.{chr}.log"
    shell:
        """
        python3 {input.script} -i {input.na_pos} -o {output.na_bed}
        """

rule split_callable:
    """ Split callability bed by samples """
    input:
        callability = "results/raw/bed/{prefix}.raw.{chr}.callable.bed"
    output:
        callability_by_sample = "results/raw/bed/{prefix}.raw.{sample}.{chr}.callable.bed"
    params:
        sample="{wildcards.sample}"
    log:
        "logs/{prefix}.{sample}.{chr}.log"
    shell:
        """
        grep {params.sample} {input.callability} | cut -f1,2,3 | bgzip > {output.callability_by_sample}
        """

rule mutlihetsep:
    """ Built input file for MSMC2 """
    input:
        vcf = "results/snps/vcf/{prefix}.SNPS.{sample}.{chr}.vcf.gz",
        na_mask = "results/snps/bed/{prefix}.{sample}.{chr}.no_na.bed",
        callability_mask = "results/raw/bed/{prefix}.raw.{sample}.{chr}.callable.bed",
        msmc2 = config["msmc2_dir"]
    output:
        multihetsep = "results/snps/ne_inference/msmc2/multihetsep/{prefix}.SNPS.{sample}.{chr}.multihetsep.txt"
    conda:
        "../envs/ne.yml"
    log:
        "logs/{prefix}.{sample}.{chr}.log"
    shell:
        """
        python3 {input.msmc2}/generate_multihetsep.py \
            --mask {input.callability_mask} \
            --negative_mask {input.na_mask} \
            {input.vcf} \
            > {output.multihetsep}
        """

rule run_msmc2:
    """ Run MSMC2 """
    input:
        multihetsep = expand("results/snps/ne_inference/msmc2/multihetsep/{{prefix}}.SNPS.{{sample}}.{chr}.multihetsep.txt", chr=chromosomes)
    output:
        msmc2_results = "results/snps/ne_inference/msmc2/inference/{prefix}.SNPS.{sample}"
    conda:
        "../envs/ne.yml"
    shell:
        """
        echo {input.multihetsep}
        msmc2_Linux -p {time_segments} -o {output.msmc2_results} {input.multihetsep}
        """