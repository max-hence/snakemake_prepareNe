include: "common.smk"

rule split_recmap:
    "Divides rec map by chromsome"

    input:
        recmap = config["recmap"]
    output:
        recmap_by_chrom = "results/map/{prefix}.{chr}.rec"
    shell:
        """
            awk -v k="{wildcards.chr}" '$2 == k' "{input.recmap}" > {output.recmap_by_chrom}
        """

rule exclude_regions:
    "Creates bed file to remove regions upper or lower a given recombination threshold"

    input:
        recmap = "results/map/{prefix}.{chr}.rec",
        fai = "results/stats/snps_na/{prefix}.SNPS.NA.{chr}.fai",
        script_bed = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/map2bed.py",
        script_fai = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/rescale_maplen.py"
    output:
        bed_1 = "results/bed/rec/{prefix}.rec1.{chr}.bed", # lower than 33%
        bed_2 = "results/bed/rec/{prefix}.rec2.{chr}.bed", # between 33 and 66%
        bed_3 = "results/bed/rec/{prefix}.rec3.{chr}.bed", # higher than 66%
        fai_1 = "results/stats/rec/{prefix}.rec1.{chr}.fai",
        fai_2 = "results/stats/rec/{prefix}.rec2.{chr}.fai",
        fai_3 = "results/stats/rec/{prefix}.rec3.{chr}.fai"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
            python3 {input.script_bed} -i {input.recmap} -l 0 -u 33 -o {output.bed_1}
            python3 {input.script_bed} -i {input.recmap} -l 33 -u 66 -o {output.bed_2}
            python3 {input.script_bed} -i {input.recmap} -l 66 -u 100 -o {output.bed_3}

            python3 {input.script_fai} -i {output.bed_1} -f {input.fai} -o {output.fai_1}
            python3 {input.script_fai} -i {output.bed_2} -f {input.fai} -o {output.fai_2}
            python3 {input.script_fai} -i {output.bed_3} -f {input.fai} -o {output.fai_3}
        """

rule trim_vcf:
    """Remove regions outside rec threshold to rerun easySFS on that subtable
    """

    input:
        vcf = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz", # vcf corrected by callability
        vcf_idx = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz.tbi",
        bed_1 = "results/bed/rec/{prefix}.rec1.{chr}.bed", # lower than 33%
        bed_2 = "results/bed/rec/{prefix}.rec2.{chr}.bed", # between 33 and 66%
        bed_3 = "results/bed/rec/{prefix}.rec3.{chr}.bed", # higher than 66%
    output:
        trimmed_vcf_1 = temp("results/vcf/rec/{prefix}.rec1.{chr}.vcf"),
        trimmed_vcf_gz_1 = "results/vcf/rec/{prefix}.rec1.{chr}.vcf.gz",
        trimmed_vcf_idx_1 = "results/vcf/rec/{prefix}.rec1.{chr}.vcf.gz.tbi",
        trimmed_vcf_2 = temp("results/vcf/rec/{prefix}.rec2.{chr}.vcf"),
        trimmed_vcf_gz_2 = "results/vcf/rec/{prefix}.rec2.{chr}.vcf.gz",
        trimmed_vcf_idx_2 = "results/vcf/rec/{prefix}.rec2.{chr}.vcf.gz.tbi",
        trimmed_vcf_3 = temp("results/vcf/rec/{prefix}.rec3.{chr}.vcf"),
        trimmed_vcf_gz_3 = "results/vcf/rec/{prefix}.rec3.{chr}.vcf.gz",
        trimmed_vcf_idx_3 = "results/vcf/rec/{prefix}.rec3.{chr}.vcf.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        bcftools view -R {input.bed_1} {input.vcf} -o {output.trimmed_vcf_1}
        bgzip -c {output.trimmed_vcf_1} > {output.trimmed_vcf_gz_1}
        tabix -p vcf {output.trimmed_vcf_gz_1}

        bcftools view -R {input.bed_2} {input.vcf} -o {output.trimmed_vcf_2}
        bgzip -c {output.trimmed_vcf_2} > {output.trimmed_vcf_gz_2}
        tabix -p vcf {output.trimmed_vcf_gz_2}

        bcftools view -R {input.bed_3} {input.vcf} -o {output.trimmed_vcf_3}
        bgzip -c {output.trimmed_vcf_3} > {output.trimmed_vcf_gz_3}
        tabix -p vcf {output.trimmed_vcf_gz_3}
        """

rule rdm_sample_vcf:
    " Randomly sampling SNPs "
    input: 
        vcf = "results/vcf/rec/{prefix}.{rec_quant}.{chr}.vcf.gz",
        vcf_idx = "results/vcf/rec/{prefix}.{rec_quant}.{chr}.vcf.gz.tbi",
        fai = "results/stats/rec/{prefix}.{rec_quant}.{chr}.fai",
    output:
        unsorted_vcf = temp("results/vcf/rec/{prefix}.{rec_quant}.320kSNP.unsorted.{chr}.vcf.gz"),
        rdm_vcf = temp("results/vcf/rec/{prefix}.{rec_quant}.320kSNP.{chr}.vcf"),
        rdm_vcf_gz = "results/vcf/rec/{prefix}.{rec_quant}.320kSNP.{chr}.vcf.gz",
        rdm_vcf_idx = "results/vcf/rec/{prefix}.{rec_quant}.320kSNP.{chr}.vcf.gz.tbi",
        rdm_fai = "results/stats/rec/{prefix}.{rec_quant}.320kSNP.{chr}.fai",
    conda:
        "../envs/vcf_processing.yml"
    params: 
      nsnps=320000
    shell:
        """
        bcftools view --header-only {input.vcf} > {output.unsorted_vcf}
        bcftools view --no-header {input.vcf} | \
            awk '{{ printf("%f\\t%s\\n",rand(),$0);}}' | ( sort -t $'\\t'  -T . -k1,1g || true) | \
            head -n {params.nsnps} | cut -f 2-  >> {output.unsorted_vcf}
            bcftools sort -o {output.rdm_vcf_gz}  {output.unsorted_vcf}  ## ensure output is sorted by position
        tabix -p vcf {output.rdm_vcf_gz}
        bgzip -cd {output.rdm_vcf_gz} > {output.rdm_vcf}

        # resize (total length - lost snps)
        nsnps=$(bcftools index -s {input.vcf} | cut -f3)
        lost_snps=$((nsnps - {params.nsnps}))
        awk '{{ $2 = $2 - lost_snps; OFS="\\t"; print }}' {input.fai} > {output.rdm_fai}
        """


rule sfs_projection_rec:
    "Run easySFS on the sub vcf"
    input:
        subsampled_vcf = "results/vcf/rec/{prefix}.{rec_quant}.320kSNP.{chr}.vcf",
        pop_path = "results/{prefix}.pop",
        easySFS = config["easySFS_path"],
    output:
        preview = "results/sfs/rec/{prefix}.{rec_quant}.preview.{chr}.txt",
    conda:
        "../envs/easySFS.yml"
    shell:
        """
        python3 {input.easySFS} -i {input.subsampled_vcf} -p {input.pop_path} \
        --preview -v -a > {output.preview}
        """

rule get_best_params_rec:
    """
        merge all results from easySFS run by chr
    """
    input:
        previews_1 = get_previews("results/sfs/rec/{prefix}.rec1.preview.{chr}.txt"),
        previews_2 = get_previews("results/sfs/rec/{prefix}.rec2.preview.{chr}.txt"),
        previews_3 = get_previews("results/sfs/rec/{prefix}.rec3.preview.{chr}.txt"),
        get_sfs_param = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/get_sfs_param.py"
    output:
        best_sample_1 = "results/sfs/rec/{prefix}.rec1.best_sample.txt",
        best_sample_2 = "results/sfs/rec/{prefix}.rec2.best_sample.txt",
        best_sample_3 = "results/sfs/rec/{prefix}.rec3.best_sample.txt"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        previews_1=$(echo {input.previews_1} | tr ' ' ',')
        previews_2=$(echo {input.previews_2} | tr ' ' ',')
        previews_3=$(echo {input.previews_3} | tr ' ' ',')
        
        python3 {input.get_sfs_param} -i $previews_1 -m "ml" -o {output.best_sample_1}
        python3 {input.get_sfs_param} -i $previews_2 -m "ml" -o {output.best_sample_2}
        python3 {input.get_sfs_param} -i $previews_3 -m "ml" -o {output.best_sample_3}
        """

rule intersect_beds:
    """
    Trim callable bed to keep only regions inside the given rec threshold (merging btw 2 beds)
    """
    input:
        bed = "results/bed/rec/{prefix}.{rec_quant}.{chr}.bed",
        callable_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/merge_bed.py"
    output:
        intersect_bed = "results/bed/rec/{prefix}.{rec_quant}.intersect.{chr}.bed"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        for rec in $(echo "rec1 rec2 rec3"); do
            python3 {input.script} -i "results/bed/rec/{wildcards.prefix}.$rec.{wildcards.chr}.bed" -b {input.callable_bed} \
            -o "results/bed/rec/{wildcards.prefix}.$rec.intersect.{wildcards.chr}.bed"
        done
        """


rule trim_bed_max_rec:
    """
        Remove regions in bed where less than <best_sample> indiv have been well called
        Resize chr length based on trimmed bed
    """
    input:
        #TODO ça sert plus à rien de mettre ça sous forme de fai, ça servait pour changer le header mais du coup juste chr + len is ok
        best_sample = "results/sfs/rec/{prefix}.{rec_quant}.best_sample.txt",
        rescaled_fai = "results/stats/rec/{prefix}.{rec_quant}.320kSNP.{chr}.fai",
        intersect_bed = "results/bed/rec/{prefix}.{rec_quant}.intersect.{chr}.bed",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/rescale_genlen.py"
    output:
        trimmed_bed = "results/bed/rec/ml/{prefix}.{rec_quant}.{chr}.callable.bed",
        resized_fai = "results/stats/rec/ml/{prefix}.{rec_quant}.320kSNP.{chr}.fai"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 ) / 2 ))
        awk -v n=$sampling_size '$4 >= n' {input.intersect_bed} > {output.trimmed_bed}
        
        python3 {input.script} -i {output.trimmed_bed} -f {input.rescaled_fai} \
        -o {output.resized_fai}
        """
