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
    shell:
        """ 
            bcftools query -l {input.vcf} | awk '{{print $0 "\tpop1"}}' > {output.pop_path}
        """

rule split_vcf:
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
        splitted_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed"
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
        vcf_snps_idx = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        bcftools view -Oz -m2 -M2 -v snps -o {output.vcf_snps} {input.splitted_vcf}
        bgzip -c {output.vcf_snps} > {output.vcf_snps_gz}
        tabix -p vcf {output.vcf_snps_gz}
        """


rule resize_vcf:
    """
    Rescale chr length based on lost variants when filtering for bi allelic snps
    Writes the updated length in vcf header
    """
    input:
        fai = config["fai_path"],
        raw_stats = "results/stats/raw/{prefix}.raw.{chr}.stats",
        vcf_snps = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz",
        vcf_snps_idx = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi"
    output:
        rescaled_fai = "results/stats/snps/full/{prefix}.SNPS.full.{chr}.fai",
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """ 
        snps=$(bcftools index -s {input.vcf_snps} | cut -f3)
        new_length=$(($(cut -f2 {input.raw_stats}) - $(cut -f3 {input.raw_stats}) + snps))
        awk '$1 == "{wildcards.chr}"' {input.fai} | awk -v val=$new_length 'BEGIN {{OFS="\t"}} {{ $2=val; print }}' > {output.rescaled_fai}
        """

        ######################
        ### VCF resampling ###
        ######################
 
# Resample and rescale vcf with different criteria

    ### Max snps
rule sfs_projection:
    """
        Run easySFS to run SFS projection to find the best sample/snps ratio
    """
        input:
            vcf_snps = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf",
            pop_path = "results/{prefix}.pop",
            easySFS = config["easySFS_path"],
        output:
            preview = "results/sfs/snps/{prefix}.SNPS.preview.{chr}.txt",
        conda:
            "../envs/easySFS.yml"
        shell:
            """
                python3 {input.easySFS} -i {input.vcf_snps} -p {input.pop_path} \
                --preview -v -a > {output.preview}
            """


    ### 10 indiv
rule trim_bed_small:
    """
        Remove regions in bed where less than 10 indiv have been well called
    """
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz",
        vcf_idx = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi",
        rescaled_fai = "results/stats/snps/full/{prefix}.SNPS.full.{chr}.fai",
        raw_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/rescale_genlen.py"
    output:
        trimmed_bed = "results/bed/snps/small/{prefix}.SNPS.small.{chr}.callable.bed",
        resized_fai = "results/stats/snps/small/{prefix}.SNPS.small.{chr}.fai"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        sampling_size=5
        bed_length=$(tail -1 {input.raw_bed} | cut -f3)
        awk -v n=$sampling_size '$4 >= n' {input.raw_bed} > {output.trimmed_bed}
        
        python3 {input.script} -i {output.trimmed_bed} -f {input.rescaled_fai} \
        -o {output.resized_fai}
        """

rule get_best_params:
    """
        merge all results from easySFS run by chr and find best params based on two methods
        - one maximizes nbr of snps
        - one maximizes a likelihood value (including n_snps and n_indiv as params)
    """
    input:
        previews = get_previews("results/sfs/snps/{prefix}.SNPS.preview.{chr}.txt"),
        get_sfs_param = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/get_sfs_param.py"

    output:
        best_sample_max = "results/sfs/snps/{prefix}.SNPS.best_sample_max.txt",
        best_sample_ml = "results/sfs/snps/{prefix}.SNPS.best_sample_ml.txt",
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        previews=$(echo {input.previews} | tr ' ' ',')
        python3 {input.get_sfs_param} -i $previews -m "max" -o {output.best_sample_max}
        python3 {input.get_sfs_param} -i $previews -m "ml" -o {output.best_sample_ml}
        """

rule trim_bed_max:
    """
        Remove regions in bed where less than <best_sample> indiv have been well called
        Resize chr length based on trimmed bed
    """
    input:
        best_sample = "results/sfs/snps/{prefix}.SNPS.best_sample_max.txt",
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz",
        vcf_idx = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi",
        rescaled_fai = "results/stats/snps/full/{prefix}.SNPS.full.{chr}.fai",
        raw_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/rescale_genlen.py"
    output:
        trimmed_bed = "results/bed/snps/max/{prefix}.SNPS.max.{chr}.callable.bed",
        resized_fai = "results/stats/snps/max/{prefix}.SNPS.max.{chr}.fai"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 ) / 2 ))
        bed_length=$(tail -1 {input.raw_bed} | cut -f3)
        awk -v n=$sampling_size '$4 >= n' {input.raw_bed} > {output.trimmed_bed}
        
        python3 {input.script} -i {output.trimmed_bed} -f {input.rescaled_fai} \
        -o {output.resized_fai}
        """

    ### all indiv
rule trim_bed_strict:
    """
        Remove regions in bed where less than <all_samples> indivs have been well called
    """
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz",
        vcf_idx = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi",
        rescaled_fai = "results/stats/snps/full/{prefix}.SNPS.full.{chr}.fai",
        raw_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        pop_path = "results/{prefix}.pop",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/rescale_genlen.py"
    output:
        trimmed_bed = "results/bed/snps/strict/{prefix}.SNPS.strict.{chr}.callable.bed",
        resized_fai = "results/stats/snps/strict/{prefix}.SNPS.strict.{chr}.fai"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        sampling_size=$(wc -l < {input.pop_path})
        bed_length=$(tail -1 {input.raw_bed} | cut -f3)
        awk -v n=$sampling_size '$4 >= n' {input.raw_bed} > {output.trimmed_bed}

        python3 {input.script} -i {output.trimmed_bed} -f {input.rescaled_fai} \
        -o {output.resized_fai}
        """

rule trim_bed_ml:
    """
        Remove regions in bed where less than <best_sample> indiv have been well called
        best sample is measure as the number of indiv that maximizes a log-likelihood function (that gives the sfs with the highest signal)
        Resize chr length based on trimmed bed
    """
    input:
        best_sample = "results/sfs/snps/{prefix}.SNPS.best_sample_ml.txt",
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz",
        vcf_idx = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi",
        rescaled_fai = "results/stats/snps/full/{prefix}.SNPS.full.{chr}.fai",
        raw_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/rescale_genlen.py"
    output:
        trimmed_bed = "results/bed/snps/ml/{prefix}.SNPS.ml.{chr}.callable.bed",
        resized_fai = "results/stats/snps/ml/{prefix}.SNPS.ml.{chr}.fai"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 ) / 2 ))
        bed_length=$(tail -1 {input.raw_bed} | cut -f3)
        awk -v n=$sampling_size '$4 >= n' {input.raw_bed} > {output.trimmed_bed}
        
        python3 {input.script} -i {output.trimmed_bed} -f {input.rescaled_fai} -o {output.resized_fai}
        """