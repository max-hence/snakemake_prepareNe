include: "common.smk"

    #############################
    ### Filter on Callability ###
    #############################

rule correct_genotype:
    """
        Changes wrongly called genotype into NA
    """
    input:
        vcf = "results/snps/vcf/{prefix}.SNPS.{chr}.vcf",
        splitted_bed = "results/raw/bed/{prefix}.raw.{chr}.callable.bed",
        script = workflow.source_path("../scripts/correct_genotype.py")
    output:
        corrected_vcf = temp("results/callability/vcf/{prefix}.SNPS.NA.{chr}.vcf"),
        corrected_vcf_gz = "results/callability/vcf/{prefix}.SNPS.NA.{chr}.vcf.gz",
        corrected_vcf_idx = "results/callability/vcf/{prefix}.SNPS.NA.{chr}.vcf.gz.tbi"
    conda:
        "../envs/callability.yml"
    shell:
        """
        python3 {input.script} -i {input.vcf} -b {input.splitted_bed} -o {output.corrected_vcf}
        bgzip < {output.corrected_vcf} > {output.corrected_vcf_gz} && tabix -p vcf {output.corrected_vcf_gz}
        """

    ###########
    ### SFS ###
    ###########

rule sfs_projection:
    """
        Run easySFS to run SFS projection to find the best sample/snps ratio
    """
        input:
            vcf_na = "results/callability/vcf/{prefix}.SNPS.NA.{chr}.vcf",
            pop_path = "results/{prefix}.pop",
            easySFS = config["easySFS_path"],
        output:
            preview = "results/callability/sfs/{prefix}.SNPS.NA.preview.{chr}.txt"
        conda:
            "../envs/callability.yml"
        log:
            "logs/na/{prefix}.{chr}"
        shell:
            """
                {input.easySFS} -i {input.vcf_na} -p {input.pop_path} \
                --preview -v -a > {output.preview}
            """

rule get_best_params:
    """
        merge all results from easySFS run by chr and measure best snp/ratio with ML
    """
    input:
        previews = get_previews,
        script = workflow.source_path("../scripts/get_sfs_param.py")
    output:
        best_sample = "results/callability/sfs/{prefix}.SNPS.NA.best_sample.txt"
    conda:
        "../envs/callability.yml"
    shell:
        """ 
            python3 {input.script} -i {input.previews} -m "ml" -o {output.best_sample}
        """


rule trim_bed_ml:
    """
        Remove regions in bed where less than <best_sample> indiv have been well called
        Resize chr length based on trimmed bed
    """
    input:
        best_sample = "results/callability/sfs/{prefix}.SNPS.NA.best_sample.txt",
        splitted_bed = "results/raw/bed/{prefix}.raw.{chr}.callable.bed",
        script = workflow.source_path("../scripts/rescale_genlen.py")
    output:
        trimmed_bed = "results/callability/bed/ml/{prefix}.SNPS.NA.ml.{chr}.callable.bed",
    conda:
        "../envs/callability.yml"
    shell:
        """
            sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 ) / 2 ))
            bed_length=$(tail -1 {input.splitted_bed} | cut -f3)
            awk -v n=$sampling_size '$4 >= n' {input.splitted_bed} > {output.trimmed_bed}
        """

rule trim_bed_strict:
    """
        Remove regions in bed where not every indiv have been well called
    """
    input:
        splitted_bed = "results/raw/bed/{prefix}.raw.{chr}.callable.bed",
        pop_path = "results/{prefix}.pop",
        script = workflow.source_path("../scripts/rescale_genlen.py")
    output:
        trimmed_bed = "results/callability/bed/strict/{prefix}.SNPS.NA.strict.{chr}.callable.bed",
    conda:
        "../envs/callability.yml"
    shell:
        """
            sampling_size=$(wc -l < {input.pop_path})
            awk -v n=$sampling_size '$4 >= n' {input.splitted_bed} > {output.trimmed_bed}
        """

rule trim_bed_small:
    """
        Remove regions in bed where less than 10 indiv have been well called
    """
    input:
        splitted_bed = "results/raw/bed/{prefix}.raw.{chr}.callable.bed",
        script = workflow.source_path("../scripts/rescale_genlen.py")
    output:
        trimmed_bed = "results/callability/bed/small/{prefix}.SNPS.NA.small.{chr}.callable.bed",
    conda:
        "../envs/callability.yml"
    shell:
        """
            sampling_size=10
            awk -v n=$sampling_size '$4 >= n' {input.splitted_bed} > {output.trimmed_bed}
        """

rule resize_chr:
    """
        Resize chr length based on trimmed bed
    """
    input:
        trimmed_bed = "results/callability/bed/{call_filter}/{prefix}.SNPS.NA.{call_filter}.{chr}.callable.bed",
        fai = "results/snps/stats/{prefix}.SNPS.resized.{chr}.fai",
        script = workflow.source_path("../scripts/rescale_genlen.py")
    output:
        resized_fai = "results/callability/stats/{prefix}.SNPS.NA.{call_filter}.{chr}.fai"
    conda:
        "../envs/callability.yml"
    shell:
        """
        python3 {input.script} -i {input.trimmed_bed} -f {input.fai} \
            -o {output.resized_fai} --method "bed"
        """