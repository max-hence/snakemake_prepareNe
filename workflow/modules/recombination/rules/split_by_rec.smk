include: "common.smk"

rule split_recmap:
    "Divides rec map by chromsome"

    input:
        recmap = config["recmap"]
    output:
        recmap_by_chrom = "results/map/{prefix}.{chr}.rec"
    shell:
        """
        awk -v k={wildcards.chr} '$2 == k' {input.recmap} > {output.recmap_by_chrom}
        """

rule split_by_rec:
    "Split recombination map by quantiles of rec rate"

    input:
        recmap = "results/map/{prefix}.{chr}.rec",
        script = workflow.source_path("../scripts/map2bed.py")
    output:
        bed_1 = "results/rec/bed/{prefix}.rec1.{chr}.bed", # lower than 33%
        bed_2 = "results/rec/bed/{prefix}.rec2.{chr}.bed", # between 33 and 66%
        bed_3 = "results/rec/bed/{prefix}.rec3.{chr}.bed", # higher than 66%
    conda:
        "../envs/recombination.yml"
    shell:
        """
        python3 {input.script} -i {input.recmap} -l 0 -u 33 -o {output.bed_1}
        python3 {input.script} -i {input.recmap} -l 33 -u 66 -o {output.bed_2}
        python3 {input.script} -i {input.recmap} -l 66 -u 100 -o {output.bed_3}
        """

rule resize_chr:
        input:
            bed = "results/rec/bed/{prefix}.{rec}.{chr}.bed",
            fai = "results/stats/snps_na/{prefix}.SNPS.NA.{chr}.fai",
            script = workflow.source_path("../scripts/rescale_genlen.py")
        output:
            fai = "results/rec/stats/{prefix}.{rec}.{chr}.fai",
        conda:
            "../envs/recombination.yml"
        shell:
            """
            python3 {input.script} -i {input.bed} -f {input.fai} -o {output.fai} --method "bed"
            """

rule split_vcf:
    "Remove regions outside rec threshold to rerun easySFS on that subtable"

    input:
        vcf = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz", # vcf corrected by callability
        vcf_idx = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz.tbi",
        bed = "results/rec/bed/{prefix}.{rec}.{chr}.bed"
    output:
        splitted_vcf = temp("results/rec/vcf/{prefix}.{rec}.{chr}.vcf"),
        splitted_vcf_gz = "results/rec/vcf/{prefix}.{rec}.{chr}.vcf.gz",
        splitted_vcf_idx = "results/rec/vcf/{prefix}.{rec}.{chr}.vcf.gz.tbi",
        stats = "results/rec/stats/{prefix}.{rec}.{chr}.stats"
    conda:
        "../envs/recombination.yml"
    log:
        "logs/rec/{prefix}.{rec}.{chr}.log"
    shell:
        """
        bcftools view -R {input.bed} {input.vcf} -o {output.splitted_vcf}
        bgzip -c {output.splitted_vcf} > {output.splitted_vcf_gz}
        tabix -p vcf {output.splitted_vcf_gz}

        bcftools index -s {output.splitted_vcf_gz} > {output.stats}
        """

rule rdm_sample_vcf:
    " Randomly sampling SNPs "
    input:
        vcf = "results/rec/vcf/{prefix}.{rec}.{chr}.vcf.gz",
        vcf_idx = "results/rec/vcf/{prefix}.{rec}.{chr}.vcf.gz.tbi",
        fai = "results/rec/stats/{prefix}.{rec}.{chr}.fai",
        snps_stats = "results/rec/stats/{prefix}.{rec}.{chr}.stats",
        script = workflow.source_path("../scripts/rescale_genlen.py")
    output:
        unsorted_vcf = temp("results/rec/vcf/{prefix}.{rec}.rdmSNP.unsorted.{chr}.vcf.gz"),
        rdm_vcf = temp("results/rec/vcf/{prefix}.{rec}.rdmSNP.{chr}.vcf"),
        rdm_vcf_gz = "results/rec/vcf/{prefix}.{rec}.rdmSNP.{chr}.vcf.gz",
        rdm_vcf_idx = "results/rec/vcf/{prefix}.{rec}.rdmSNP.{chr}.vcf.gz.tbi",
        rdm_stats = "results/rec/stats/{prefix}.{rec}.rdmSNP.{chr}.stats",
        rdm_fai = "results/rec/stats/{prefix}.{rec}.rdmSNP.{chr}.fai",
    conda:
        "../envs/recombination.yml"
    params:
      nsnps=250000
    shell:
        """
        bcftools view --header-only {input.vcf} > {output.unsorted_vcf}
        bcftools view --no-header {input.vcf} | \
            awk '{{ printf("%f\\t%s\\n",rand(),$0);}}' | ( sort -t $'\\t'  -T . -k1,1g || true) | \
            head -n {params.nsnps} | cut -f 2-  >> {output.unsorted_vcf}
            bcftools sort -o {output.rdm_vcf_gz}  {output.unsorted_vcf}  ## ensure output is sorted by position
        tabix -p vcf {output.rdm_vcf_gz}
        bgzip -cd {output.rdm_vcf_gz} > {output.rdm_vcf}

        # Chromosome rescaling
        bcftools index -s {output.rdm_vcf_gz} > {output.rdm_stats}
        python3 {input.script} -i {input.snps_stats} -r {output.rdm_stats} -f {input.fai} -o {output.rdm_fai} --method snp
        """

    ###################
    ### Subsampling ###
    ###################

rule sfs_projection:
    "Run easySFS on the sub vcf"
    input:
        subsampled_vcf = "results/rec/vcf/{prefix}.{rec}.rdmSNP.{chr}.vcf",
        pop_path = "results/{prefix}.pop",
        easySFS = config["easySFS_path"],
    output:
        preview = "results/rec/sfs/{prefix}.{rec}.preview.{chr}.txt",
    conda:
        "../envs/recombination.yml"
    shell:
        """
        python3 {input.easySFS} -i {input.subsampled_vcf} -p {input.pop_path} \
        --preview -v -a > {output.preview}
        """

rule get_best_params:
    """
        merge all results from easySFS run by chr
    """
    input:
        previews = get_previews("results/rec/sfs/{prefix}.{rec}.preview.{chr}.txt"),
        get_sfs_param = workflow.source_path("../scripts/get_sfs_param.py")
    output:
        best_sample = "results/rec/sfs/{prefix}.{rec}.best_sample.txt"
    conda:
        "../envs/recombination.yml"
    shell:
        """
        python3 {input.get_sfs_param} -i {input.previews} -m "ml" -o {output.best_sample}
        """

rule intersect_beds:
    """
    Trim callable bed to keep only regions inside the given rec threshold (merging btw 2 beds)
    """
    input:
        bed = "results/rec/bed/{prefix}.{rec}.{chr}.bed",
        callable_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        script = workflow.source_path("../scripts/merge_bed.py")
    output:
        intersect_bed = "results/rec/bed/{prefix}.{rec}.intersect.{chr}.bed"
    conda:
        "../envs/recombination.yml"
    log:
        "logs/rec/{prefix}.{rec}.{chr}.log"
    shell:
        """
        python3 {input.script} -i {input.bed} -b {input.callable_bed} -o {output.intersect_bed}
        """

rule trim_bed:
    """
        Remove regions in bed where less than <best_sample> indiv have been well called
        Resize chr length based on trimmed bed
    """
    input:
        best_sample = "results/rec/sfs/{prefix}.{rec}.best_sample.txt",
        rescaled_fai = "results/rec/stats/{prefix}.{rec}.rdmSNP.{chr}.fai",
        intersect_bed = "results/rec/bed/{prefix}.{rec}.intersect.{chr}.bed",
        script = workflow.source_path("../scripts/rescale_genlen.py")
    output:
        trimmed_bed = "results/rec/bed/{prefix}.{rec}.{chr}.callable.bed",
        resized_fai = "results/rec/stats/{prefix}.{rec}.rdmSNP.resized.{chr}.fai"
    conda:
        "../envs/recombination.yml"
    shell:
        """
        sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 ) / 2 ))
        awk -v n=$sampling_size '$4 >= n' {input.intersect_bed} > {output.trimmed_bed}
        
        python3 {input.script} -i {output.trimmed_bed} -f {input.rescaled_fai} \
        -o {output.resized_fai} --method bed
        """
