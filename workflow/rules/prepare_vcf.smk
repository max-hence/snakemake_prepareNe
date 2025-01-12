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
        splitted_vcf = temp("results/raw/vcf/{prefix}.raw.{chr}.vcf.gz"),
        splitted_vcf_idx = temp("results/raw/vcf/{prefix}.raw.{chr}.vcf.gz.tbi"),
        stats = "results/raw/stats/{prefix}.raw.{chr}.stats"
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
        splitted_bed = temp("results/raw/bed/{prefix}.raw.{chr}.callable.bed")
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
        splitted_vcf = "results/raw/vcf/{prefix}.raw.{chr}.vcf.gz",
        splitted_vcf_idx = "results/raw/vcf/{prefix}.raw.{chr}.vcf.gz.tbi"
    output:
        vcf_snps = temp("results/snps/vcf/{prefix}.SNPS.{chr}.vcf"),
        vcf_snps_gz = "results/snps/vcf/{prefix}.SNPS.{chr}.vcf.gz",
        vcf_snps_idx = "results/snps/vcf/{prefix}.SNPS.{chr}.vcf.gz.tbi",
        snps_stats = "results/snps/stats/{prefix}.SNPS.{chr}.stats"
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

    ###########
    ### SFS ###
    ###########

rule sfs_projection:
    """
        Run easySFS to run SFS projection for each sample size (2:2n)
    """
    input:
        vcf_snps = "results/snps/vcf/{prefix}.SNPS.{chr}.vcf",
        pop_path = "results/{prefix}.pop",
        easySFS = config["easySFS_path"],
    output:
        preview = "results/snps/sfs/{prefix}.SNPS.preview.{chr}.txt",
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
        previews = get_previews("results/snps/sfs/{prefix}.SNPS.preview.{chr}.txt"),
        script = workflow.source_path("../scripts/get_sfs_param.py")
    output:
        best_sample_size = "results/snps/sfs/{prefix}.SNPS.best_sample.txt"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.log"
    shell:
        """
            python3 {input.script} -i {input.previews} --method ml -o {output.best_sample_size}
        """


rule trim_bed:
    """
        Remove regions in bed where less than <best_sample> indiv have been well called
        Resize chr length based on trimmed bed
    """
        input:
            best_sample = "results/snps/sfs/{prefix}.SNPS.best_sample.txt",
            fai = config["fai_path"],
            bed = "results/raw/bed/{prefix}.raw.{chr}.callable.bed",
            script = workflow.source_path("../scripts/rescale_genlen.py")
        output:
            trimmed_bed = "results/snps/bed/{prefix}.SNPS.{chr}.callable.bed",
            resized_fai = "results/snps/stats/{prefix}.SNPS.resized.{chr}.fai"
        conda:
            "../envs/vcf_processing.yml"
        shell:
            """
                sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 ) / 2 ))
                awk -v n=$sampling_size '$4 >= n' {input.bed} > {output.trimmed_bed}

                python3 {input.script} -i {output.trimmed_bed} -f {input.fai} \
                -o {output.resized_fai} --method "bed"
            """