rule raw_pi:
    input:
        vcf_snps = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz",
        vcf_snps_idx = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi",
    output:
        sites_pi = "results/genetic_div/{prefix}.SNPS.{chr}.sites.pi"
        wdw_pi = "results/genetic_div/{prefix}.SNPS.{chr}.windowed.pi"
    conda:
        "../envs/genetic_div.yml"
    log:
        "logs/genetic_div/{prefix}.{chr}.txt"
    shell:
        """
        vcftools --gzvcf {input.vcf} --site-pi --out "results/genetic_div/{wildcards.prefix}.SNPS.{wildcards.chr}"
        vcftools --gzvcf {input.vcf} --window-pi 50000 --out "results/genetic_div/{wildcards.prefix}.SNPS.{wildcards.chr}"
        """

rule filtered_pi:
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.NA.{chr}.vcf.gz",
        vcf_idx = "results/vcf/snps/{prefix}.SNPS.NA.{chr}.vcf.gz.tbi",
    output:
        sites_pi = "results/genetic_div/{prefix}.SNPS.NA.{chr}.sites.pi"
        wdw_pi = "results/genetic_div/{prefix}.SNPS.NA.{chr}.windowed.pi"
    conda:
        "../envs/genetic_div.yml"
    log:
        "logs/genetic_div/{prefix}.{chr}.txt"
    shell:
        """
        vcftools --gzvcf {input.vcf} --site-pi --out "results/genetic_div/{wildcards.prefix}.SNPS.NA.{wildcards.chr}"
        vcftools --gzvcf {input.vcf} --window-pi 50000 --out "results/genetic_div/{wildcards.prefix}.SNPS.NA.{wildcards.chr}"
        """

    # Pi by gene
rule get_genes:
    input:
        gff = "results/genome/{prefix}.{chr}.gff",
        scripts = workflow.source_path("../scripts/get_genes.py")
    output:
        bed = "results/genome/{prefix}.{chr}.genes.bed"
    conda:
        "../envs/genetic_div.yml"
    log:
    shell:
        """
        python3 {input.script} -i {input.gff} -o {output.bed}
        """

rule pi_by_genes:
    input:
        pi = "results/genetic_div/{prefix}.SNPS.{chr}.sites.pi"
        bed = "results/geneD/{prefix}.{chr}.genes.bed"
        script = workflow.source_path("../scripts/pi_by_genes.py")
    output:
        pi_genes = 
    conda:
        "../envs/genetic_div.yml"
    log:
    shell:
        """
        python3 {input.script} --pi \
        -i {input.pi} \
        -s {input.bed} \
        -o {output.pi_genes}
        """

rule split_genome:
    "divide reference genome by chromosomes"
    input:
        genome = config["genome_path"]
    output:
        genome_by_chr = "results/genome/{prefix}.{chr}.fna"
    log:

    shell:
        """
        awk '/^>{wildcards.chr}/ {print; found=1; next} /^>/ {found=0} found' {input.genome} > {output.genome_by_chr}
        """

rule degeneracy:
    input:
        gff = "results/genome/{prefix}.{chr}.gff",
        genome = "results/genome/{prefix}.{chr}.fna",
        script = config["degenotate_path"]
    output:
        degeneracy_dir = dir("results/genetic_div/degeneracy/{prefix}.{chr}"),
        degeneracy_04 = "results/genetic_div/degeneracy/{prefix}.{chr}.degeneracy-04-sites.bed"
    conda:
        "../envs/genetic_div.yml"
    log:
    shell:
        """
        python3 {input.script} \
		-a {input.gff} \
		-g {input.genome} \
		-o {output.degeneracy_dir} -d " "

        awk '$5 == 4 || $5 == 0' {output.degeneracy_dir}/degeneracy-all-sites.bed > {output.degeneracy_04}
        """

rule get_degen_summary:
    input:
        transcript_count = "results/genetic_div/degeneracy/{prefix}.{chr}.transcript-counts.tsv",
        bed = "results/genome/{prefix}.{chr}.genes.bed"
    output:
        yes = temp("results/genetic_div/degeneracy/{prefix}.{chr}.transcript-counts.yes.tsv"),
        tmp = temp("results/genetic_div/degeneracy/{prefix}.{chr}.tmp"),
        degen_summary = "results/genetic_div/degeneracy/{prefix}.{chr}.degeneracy_summary.tsv"
    log:
        "logs/genetic_div/{prefix}.{chr}.log"
    shell:
        """
        grep "yes" {input.transcript_count} > {output.yes}
        awk 'NR==FNR{a[$2]; next} ($4 in a)' {output.yes} {input.bed} > {output.tmp} \
		&& paste tmp <(cut -f3-9 {output.yes}) > {output.degen_summary}
        """

rule p0p4:
    "Measure p0/p4"
    input:
        pi = "results/genetic_div/{prefix}.SNPS.{chr}.sites.pi",
        degeneracy_04 = "results/genetic_div/{prefix}/degeneracy-04-sites.bed",
        degeneracy_summary = "results/genetic_div/degeneracy/{prefix}.{chr}.degeneracy_summary.tsv",
        script = workflow.source_path("../scripts/pi_by_genes.py")
    output:
        p0p4 = "results/genetic_div/{prefix}.SNPS.{chr}.p0p4"
    conda:
        "../envs/genetic_div.yml"
    log:
    shell:
        """
        python3 {input.script} --p0p4 \
            -i {input.pi} \
            -d {input.degeneracy_04} \
            -s {input.degeneracy.summary} \
            -o {p0p4}
        """
