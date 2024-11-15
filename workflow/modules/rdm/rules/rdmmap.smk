### Rdm sampling for control

include: "common.smk"

rule rdm_map:
    input:
        recmap = config["recmap"],
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/modules/rdm/scripts/rdm_recmask.py"
    output:
        tmp_map = "results/tmp/map.tmp"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        python3 {input.script} -i {input.recmap} -o {output.tmp_map}
        """

rule divide_map:
    input:
        tmp_map = "results/tmp/map.tmp",
    output:
        rdm_map_by_chr = "results/map/rdm/{prefix}.rdm.{chr}.rec"
    shell:
        """
        awk -v k="{wildcards.chr}" '$2 == k' "{input.tmp_map}" > {output.rdm_map_by_chr}
        """


rule map2bed:
    " Transform map to bed and rescale "
    input:
        rdm_map = "results/map/rdm/{prefix}.rdm.{chr}.rec",
        fai = "results/stats/snps_na/{prefix}.SNPS.NA.{chr}.fai",
        script_fai = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/rescale_maplen.py"
    output:
        rdm_bed = "results/bed/rdm/{prefix}.rdm.{chr}.bed",
        rescaled_fai = "results/stats/rdm/{prefix}.rdm.{chr}.fai"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        cut -d ' ' -f 2-4 {input.rdm_map} | tr ' ' '\t' > {output.rdm_bed}
        python3 {input.script_fai} -i {output.rdm_bed} -f {input.fai} -o {output.rescaled_fai}
        """

rule trim_vcf_rdm:
    input:
        vcf = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz", # vcf corrected by callability
        vcf_idx = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz.tbi",
        bed = "results/bed/rdm/{prefix}.rdm.{chr}.bed",
    output:
        trimmed_vcf = "results/vcf/rdm/{prefix}.rdm.{chr}.vcf",
        trimmed_vcf_gz = "results/vcf/rdm/{prefix}.rdm.{chr}.vcf.gz",
        trimmed_vcf_gz_tbi = "results/vcf/rdm/{prefix}.rdm.{chr}.vcf.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        bcftools view -R {input.bed} {input.vcf} -o {output.trimmed_vcf}
        bgzip -c {output.trimmed_vcf} > {output.trimmed_vcf_gz}
        tabix -p vcf {output.trimmed_vcf_gz}
        """

rule sfs_projection_rdm:
    "Run easySFS on the sub vcf"
    input:
        trimmed_vcf = "results/vcf/rdm/{prefix}.rdm.{chr}.vcf",
        pop_path = "results/{prefix}.pop",
        easySFS = config["easySFS_path"],
    output:
        preview = "results/sfs/rdm/{prefix}.rdm.{chr}.preview",
    conda:
        "../envs/easySFS.yml"
    shell:
        """
        python3 {input.easySFS} -i {input.trimmed_vcf} -p {input.pop_path} \
        --preview -v -a > {output.preview}
        """

rule get_best_params_rdm:
    """
        merge all results from easySFS run by chr
    """
    input:
        previews = get_previews("results/sfs/rdm/{prefix}.rdm.{chr}.preview"),
        get_sfs_param = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/get_sfs_param.py"
    output:
        best_sample = "results/sfs/rdm/{prefix}.rdm.best_sample.txt"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        python3 {input.get_sfs_param} -i $(echo {input.previews} | tr ' ' ',') -m "max" -o {output.best_sample}
        """

rule intersect_beds_rdm:
    """
    Trim callable bed to keep only regions inside the given rec threshold
    """
    input:
        bed = "results/bed/rdm/{prefix}.rdm.{chr}.bed",
        callable_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/merge_bed.py"
    output:
        intersect_bed = "results/bed/rdm/{prefix}.rdm.intersect.{chr}.bed"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        python3 {input.script} -i {input.bed} -b {input.callable_bed} -o {output.intersect_bed}
        """

rule trim_bed_rdm:
    """
        Remove regions in bed where less than <best_sample> indiv have been well called
        Resize chr length based on trimmed bed
    """
    input:
        best_sample = "results/sfs/rdm/{prefix}.rdm.best_sample.txt",
        rescaled_fai = "results/stats/rdm/{prefix}.rdm.{chr}.fai",
        intersect_bed = "results/bed/rdm/{prefix}.rdm.intersect.{chr}.bed",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/rescale_genlen.py"
    output:
        trimmed_bed = "results/bed/rdm/{prefix}.rdm.{chr}.callable.bed",
        resized_fai = "results/stats/rdm/{prefix}.rdm.rescaled.{chr}.fai"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 ) / 2 ))
        awk -v n=$sampling_size '$4 >= n' {input.intersect_bed} > {output.trimmed_bed}

        python3 {input.script} -i {output.trimmed_bed} -f {input.rescaled_fai} \
        -o {output.resized_fai}
        """