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
        script = "/shared/projects/plant_lewontin_paradox/vcf_processing/scripts/snakemake_prepareNe/workflow/modules/gene_density/scripts/get_gene_density.py"
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
        script = "/shared/projects/plant_lewontin_paradox/vcf_processing/scripts/snakemake_prepareNe/workflow/modules/gene_density/scripts/geneD_to_bed.py"
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


rule split_vcf_by_geneD:
    input:
        vcf = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz",
        vcf_idx = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz.tbi",
        bed = "results/geneD/{prefix}.{density_cat}.{chr}.bed"
    output:
        splitted_vcf = "results/geneD/vcf/{prefix}.{density_cat}.{chr}.vcf.gz",
        splitted_vcf_idx = "results/geneD/vcf/{prefix}.{density_cat}.{chr}.vcf.gz.tbi"
    conda: 
        "../envs/gene_density.yml"
    shell:
        """
        bcftools view -R {input.bed} {input.vcf} -o {output.splitted_vcf}
        tabix -p vcf {output.splitted_vcf}
        """

rule rdm_sample_vcf:
    " Random SNPs sampling "
    input: 
        vcf = "results/geneD/vcf/{prefix}.{density_cat}.{chr}.vcf.gz",
        vcf_idx = "results/geneD/vcf/{prefix}.{density_cat}.{chr}.vcf.gz.tbi",
        fai = "results/stats/rec/{prefix}.{density_cat}.{chr}.fai",
    output:
        unsorted_vcf = temp("results/geneD/vcf/{prefix}.{density_cat}.rdmSNP.unsorted.{chr}.vcf.gz"),
        rdm_vcf = temp("results/geneD/vcf/{prefix}.{density_cat}.rdmSNP.{chr}.vcf"),
        rdm_vcf_gz = "results/geneD/vcf/{prefix}.{density_cat}.rdmSNP.{chr}.vcf.gz",
        rdm_vcf_idx = "results/geneD/vcf/{prefix}.{density_cat}.rdmSNP.{chr}.vcf.gz.tbi",
        rdm_fai = "results/geneD/stats/{prefix}.{density_cat}.rdmSNP.{chr}.fai",
    conda:
        "../envs/gene_density.yml"
    params: 
      nsnps=320000 # ?
    shell:
        """
        bcftools view --header-only {input.vcf} > {output.unsorted_vcf}
        bcftools view --no-header {input.vcf} | \
            awk '{{ printf("%f\\t%s\\n",rand(),$0);}}' | ( sort -t $'\\t'  -T . -k1,1g || true) | \
            head -n {params.nsnps} | cut -f 2-  >> {output.unsorted_vcf}
            bcftools sort -o {output.rdm_vcf_gz}  {output.unsorted_vcf}  ## ensure output is sorted by position
        tabix -p vcf {output.rdm_vcf_gz}
        bgzip -cd {output.rdm_vcf_gz} > {output.rdm_vcf} # for next rules

        # resize (total length - lost snps) #TODO: PAS UNE SOUSTRACTION MAIS UN PRODUIT EN x ici !
        total_snps=$(bcftools index -s {input.vcf} | cut -f3)
        lost_snps=$((total_snps - {params.nsnps}))
        awk '{{ $2 = $2 - lost_snps; OFS="\\t"; print }}' {input.fai} > {output.rdm_fai}
        """