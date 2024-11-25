# Split vcf by gene density based on genome annotation

rule split_gff:
    " Split gff by chr"

    input:
        gff = config["annotation"]
    output:
        gff_by_chr = temp("results/geneD/{prefix}.{chr}.gff")
    shell:
        """
        awk -v chr={wildcards.chr} '$1 == chr' {input.gff} > {output.gff_by_chr}
        """

rule gene_density:
    " Measure gene density by window "

    input:
        gff = "results/geneD/{prefix}.{chr}.gff",
        script = workflow.source_path("../scripts/get_gene_density.py")
    output:
        density = "results/geneD/{prefix}.{chr}.geneD"
    conda:
        "../envs/gene_density.yml"
    params:
        window = config["window_size"]
    shell:
        """
        python3 {input.script} -i {input.gff} -w {params.window} -o {output.density}
        """

rule geneD_to_bed:
    input:
        density = "results/geneD/{prefix}.{chr}.geneD",
        script = workflow.source_path("../scripts/geneD_to_bed.py")
    output:
        low_bed = "results/geneD/{prefix}.lowD.{chr}.bed",
        high_bed = "results/geneD/{prefix}.highD.{chr}.bed"
    conda:
        "../envs/gene_density.yml"
    shell:
        """
        python3 {input.script} -i {input.density} -o {output.low_bed} -l 0 -u 45
        python3 {input.script} -i {input.density} -o {output.low_bed} -l 55 -u 100
        """

rule resize_chr:
    "Measure chr size from new bed file"
    input:
        bed = "results/geneD/{prefix}.{density_cat}.{chr}.bed",
        fai = "results/stats/snps_na/{prefix}.SNPS.NA.{chr}.fai",
        script = workflow.source_path("../scripts/rescale_genlen.py")
    output:
        fai = "results/stats/rec/{prefix}.{density_cat}.{chr}.fai",
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        python3 {input.script} -i {input.bed} -f {input.fai} -o {output.fai} --method "bed"
        """


rule split_vcf_by_geneD:
    input:
        vcf = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz",
        vcf_idx = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz.tbi",
        bed = "results/geneD/{prefix}.{density_cat}.{chr}.bed"
    output:
        splitted_vcf = "results/geneD/vcf/{prefix}.{density_cat}.{chr}.vcf.gz",
        splitted_vcf_idx = "results/geneD/vcf/{prefix}.{density_cat}.{chr}.vcf.gz.tbi",
        stats = "results/stats/rec/{prefix}.{density_cat}.{chr}.stats"
    conda: 
        "../envs/gene_density.yml"
    shell:
        """
        bcftools view -R {input.bed} {input.vcf} -o {output.splitted_vcf}
        tabix -p vcf {output.splitted_vcf}

        bcftools index -s {output.splitted_vcf} > {output.stats}
        """

rule rdm_sample_vcf:
    " Random SNPs sampling "
    input: 
        vcf = "results/geneD/vcf/{prefix}.{density_cat}.{chr}.vcf.gz",
        vcf_idx = "results/geneD/vcf/{prefix}.{density_cat}.{chr}.vcf.gz.tbi",
        fai = "results/stats/rec/{prefix}.{density_cat}.{chr}.fai",
        stats = "results/stats/rec/{prefix}.{density_cat}.{chr}.stats"
    output:
        unsorted_vcf = temp("results/geneD/vcf/{prefix}.{density_cat}.rdmSNP.unsorted.{chr}.vcf.gz"),
        rdm_vcf = temp("results/geneD/vcf/{prefix}.{density_cat}.rdmSNP.{chr}.vcf"),
        rdm_vcf_gz = "results/geneD/vcf/{prefix}.{density_cat}.rdmSNP.{chr}.vcf.gz",
        rdm_vcf_idx = "results/geneD/vcf/{prefix}.{density_cat}.rdmSNP.{chr}.vcf.gz.tbi",
        rdm_stats = temp("results/geneD/stats/{prefix}.{density_cat}.rdmSNP.{chr}.stats"),
        rdm_fai = "results/geneD/stats/{prefix}.{density_cat}.rdmSNP.{chr}.fai",
    conda:
        "../envs/gene_density.yml"
    params: 
      nsnps=100000
    shell:
        """
        bcftools view --header-only {input.vcf} > {output.unsorted_vcf}
        bcftools view --no-header {input.vcf} | \
            awk '{{ printf("%f\\t%s\\n",rand(),$0);}}' | ( sort -t $'\\t'  -T . -k1,1g || true) | \
            head -n {params.nsnps} | cut -f 2-  >> {output.unsorted_vcf}
            bcftools sort -o {output.rdm_vcf_gz}  {output.unsorted_vcf}  ## ensure output is sorted by position
        tabix -p vcf {output.rdm_vcf_gz}
        bgzip -cd {output.rdm_vcf_gz} > {output.rdm_vcf} # for next rules

        bcftools index -s {output.rdm_vcf_gz} > {output.rdm_stats}
        python3 {input.script} -i {input.stats} -r {output.rdm_stats} -f {input.fai} -o {output.rdm_fai} --method snp
        """