# Correct genotypes based on callability bed file 
# Run easy_SFS to determine the best indiv number to take to maximize SNPs and minimize NA
# Then resample the bed file

include: "common.smk"
rule prepare_pop:
    """
        Writes pop file which will be used on a lot of following rules
    """
    input:
        vcf = config["vcf_path"]
    output:
        pop_path = "results/{prefix}.pop"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.log"
    shell:
        """ 
            bcftools query -l {input.vcf} | awk '{{print $0 "\tpop1"}}' > {output.pop_path}
        """

rule split_by_chr:
    """
        split vcf by chr
    """
    input:
        vcf = config["vcf_path"]
    output:
        splitted_vcf = temp("results/vcf/raw/{prefix}.raw.{chr}.vcf.gz"),
        splitted_vcf_idx = temp("results/vcf/raw/{prefix}.raw.{chr}.vcf.gz.tbi"),
        stats = "results/stats/raw/{prefix}.raw.{chr}.stats"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.{chr}.log"
    shell:
        """
            bcftools view -O z -r {wildcards.chr} -o {output.splitted_vcf} {input.vcf}
            tabix -p vcf {output.splitted_vcf}
            bcftools index -s {output.splitted_vcf} > {output.stats}
        """

rule split_bed:
    input:
        bed = config["bed_path"]
    output:
        splitted_bed = temp("results/bed/raw/{prefix}.raw.{chr}.callable.bed")
    log:
        "logs/{prefix}.{chr}.log"
    shell:
        """
            awk -v k="{wildcards.chr}" '$1 == k' "{input.bed}" > {output.splitted_bed}
        """

    #################################
    ### Filter on Bi-allelic SNPs ###
    #################################

rule filter_snps:
    """
        Remove Multiallelic, indels and MNP
    """
    input:
        splitted_vcf = "results/vcf/raw/{prefix}.raw.{chr}.vcf.gz",
        splitted_vcf_idx = "results/vcf/raw/{prefix}.raw.{chr}.vcf.gz.tbi"
    output:
        vcf_snps = temp("results/vcf/snps/{prefix}.SNPS.{chr}.vcf"),
        vcf_snps_gz = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz",
        vcf_snps_idx = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi",
        snps_stats = "results/stats/snps/full/{prefix}.SNPS.{chr}.stats"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.{chr}.log"
    shell:
        """
        bcftools view -Oz -m2 -M2 -v snps -o {output.vcf_snps} {input.splitted_vcf}
        bgzip -c {output.vcf_snps} > {output.vcf_snps_gz}
        tabix -p vcf {output.vcf_snps_gz}
        bcftools index -s {output.vcf_snps_gz} > {output.snps_stats}
        """

        ######################
        ### VCF resampling ###
        ######################
 
# Resample and rescale vcf with different criteria

    ### Max snps
rule sfs_projection:
    """
        Run easySFS to run SFS projection for each sample size (2:2n)
    """
    input:
        vcf_snps = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf",
        pop_path = "results/{prefix}.pop",
        easySFS = config["easySFS_path"],
    output:
        preview = "results/sfs/snps/{prefix}.SNPS.preview.{chr}.txt",
    conda:
        "../envs/easySFS.yml"
    log:
        "logs/{prefix}.{chr}.log"
    shell:
        """
            python3 {input.easySFS} -i {input.vcf_snps} -p {input.pop_path} \
            --preview -v -a > {output.preview}
        """

rule get_best_params:
    """
        merge all results from easySFS run by chr and find best params based on two methods
        - one maximizes nbr of snps
        - one maximizes a likelihood value (including n_snps and n_indiv as params)
    """
    input:
        previews = get_previews("results/sfs/snps/{prefix}.SNPS.preview.{chr}.txt"),
        script = workflow.source_path("../scripts/get_sfs_param.py")
    output:
        best_sample_size = "results/sfs/snps/{prefix}.SNPS.best_sample_{sfs_params_method}.txt"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.{sfs_params_method}.log"
    params:
        method = lambda wildcards: wildcards.sfs_params_method
    shell:
        """
            python3 {input.script} -i {input.previews} --method {params.method} -o {output.best_sample_size}
        """

    # Max SNPs
rule trim_bed:
    """
        Remove regions in bed where less than <n> samples have been well called
    """
    input:
        raw_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed"
    output:
        trimmed_bed = "results/bed/snps/{subsample}/{prefix}.SNPS.{subsample}.{chr}.callable.bed"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.{subsample}.{chr}.log"
    shell:
        """
        total_sample=$(wc -l {input.pop})
        sampling_size=$(({wildcards.subsample}*$total_sample/100))
        awk -v n=$sampling_size '$4 >= n' {input.raw_bed} > {output.trimmed_bed}
        """

rule trim_chr:
    """
    Trim chromosome
    """
    input:
        vcf = 
        bed = "results/bed/snps/{subsample}/{prefix}.SNPS.{subsample}.{chr}.callable.bed"
    output:
        trimmed_bed =
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.{subsample}.{chr}.log"
    shell:
        """
        bcftools view ...
        """

rule resize_chr:
    """
    Change chr size based on new bed
    """
    input:
        trimmed_bed = "results/bed/snps/{subsample}/{prefix}.SNPS.{subsample}.{chr}.callable.bed",
        fai = config["fai_path"],
        script = workflow.source_path("../scripts/rescale_genlen.py")
    output:
        rescaled_fai = "results/stats/snps/{subsample}/{prefix}.SNPS.{subsample}.{chr}.fai"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.{subsample}.{chr}.log"
    shell:
        """
        python3 {input.script} -i {input.trimmed_bed} -f {input.fai} -o {output.rescaled_fai} --method bed
        """