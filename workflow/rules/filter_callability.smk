include: "common.smk"

    #############################
    ### Filter on Callability ###
    #############################

rule correct_genotype:
    """
        Changes wrongly called genotype into NA
    """
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf",
        splitted_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        script = workflow.source_path("../scripts/correct_genotype.py")
    output:
        corrected_vcf = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf",
        corrected_vcf_gz = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz",
        corrected_vcf_idx = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        python3 {input.script} -i {input.vcf} -b {input.splitted_bed} -o {output.corrected_vcf}
        bgzip < {output.corrected_vcf} > {output.corrected_vcf_gz} && tabix -p vcf {output.corrected_vcf_gz}
        """


rule sfs_projection_na:
    """
        Run easySFS to run SFS projection to find the best sample/snps ratio
    """
        input:
            vcf_na = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf",
            pop_path = "results/{prefix}.pop",
            easySFS = config["easySFS_path"],
        output:
            preview = "results/sfs/snps_na/{prefix}.SNPS.NA.preview.{chr}.txt"
        conda:
            "../envs/easySFS.yml"
        shell:
            """
                {input.easySFS} -i {input.vcf_na} -p {input.pop_path} \
                --preview -v -a > {output.preview}
            """

rule get_best_params_na:
    """
        merge all results from easySFS run by chr and measure best snp/ratio with ML
    """
    input:
        previews = get_previews("results/sfs/snps_na/{prefix}.SNPS.NA.preview.{chr}.txt"),
        script = workflow.source_path("../scripts/get_sfs_param.py")
    output:
        best_sample = "results/sfs/snps_na/{prefix}.SNPS.NA.best_sample.txt"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """ 
            python3 {input.script} -i {input.previews} -m "ml" -o {output.best_sample}
        """


rule trim_bed_na:
    """
        Remove regions in bed where less than <best_sample> indiv have been well called
        Resize chr length based on trimmed bed
    """
        input:
            best_sample = "results/sfs/snps_na/{prefix}.SNPS.NA.best_sample.txt",
            corrected_vcf = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz",
            corrected_vcf_idx = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz.tbi",
            rescaled_fai = "results/stats/snps/full/{prefix}.SNPS.full.{chr}.fai",
            splitted_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
            script = workflow.source_path("../scripts/rescale_genlen.py")
        output:
            trimmed_bed = "results/bed/snps_na/{prefix}.SNPS.NA.{chr}.callable.bed",
            resized_fai = "results/stats/snps_na/{prefix}.SNPS.NA.{chr}.fai"
        conda:
            "../envs/vcf_processing.yml"
        shell:
            """
                sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 ) / 2 ))
                bed_length=$(tail -1 {input.splitted_bed} | cut -f3)
                awk -v n=$sampling_size '$4 >= n' {input.splitted_bed} > {output.trimmed_bed}

                python3 {input.script} -i {output.trimmed_bed} -f {input.rescaled_fai} \
                -o {output.resized_fai} --method "bed"
            """