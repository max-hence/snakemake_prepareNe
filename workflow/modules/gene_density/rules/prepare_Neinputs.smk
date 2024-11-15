# Prepare inputs files for stairway_plot and smc++

rule sfs_projection:
    """
        Run easySFS to do SFS projection and find the best sample/snps ratio
    """
        input:
            vcf = "results/vcf/geneD/{prefix}.{density_cat}.rdmSNP.{chr}.vcf",
            pop_path = "results/{prefix}.pop",
            easySFS = config["easySFS_path"],
        output:
            preview = "results/sfs/snps/{prefix}.{density_cat}.rdmSNP.preview.{chr}.txt",
        conda:
            "../envs/easySFS.yml"
        shell:
            """
                python3 {input.easySFS} -i {input.vcf} -p {input.pop_path} \
                --preview -v -a > {output.preview}
            """

rule get_best_params:
    """
        merge all results from easySFS run by chr and find best params based on two methods
        - one maximizes nbr of snps
        - one maximizes a likelihood value (including n_snps and n_indiv as params)
    """
    input:
        previews = get_previews("results/sfs/snps/{prefix}.{density_cat}.rdmSNP.preview.{chr}.txt"),
        get_sfs_param = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/get_sfs_param.py"

    output:
        best_sample = "results/sfs/geneD/{prefix}.geneD.best_params.txt",
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        previews=$(echo {input.previews} | tr ' ' ',')
        python3 {input.get_sfs_param} -i $previews -m "ml" -o {output.best_sample_ml}
        """


rule trim_bed:
    """
        Remove regions in callability.bed where less than <best_sample> indiv have been well called
        Resize chr length based on trimmed bed
    """
    input:
        best_sample = "results/sfs/snps/{prefix}.geneD.best_params.txt",
        vcf = "results/vcf/geneD/{prefix}.{density_cat}.rdmSNP.{chr}.vcf.gz",
        vcf_idx = "results/vcf/geneD/{prefix}.{density_cat}.rdmSNP.{chr}.vcf.gz.tbi",
        fai = "results/stats/geneD/{prefix}.{density_cat}.rdmSNP.{chr}.fai",
        raw_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/rescale_genlen.py"
    output:
        trimmed_bed = "bed/snps/max/{prefix}.{density_cat}.geneD.{chr}.callable.bed",
        rescaled_fai = "results/stats/geneD/{prefix}.{density_cat}.rdmSNP.trimmed.{chr}.fai"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 ) / 2 ))
        bed_length=$(tail -1 {input.raw_bed} | cut -f3)
        awk -v n=$sampling_size '$4 >= n' {input.raw_bed} > {output.trimmed_bed}
        
        python3 {input.script} -i {output.trimmed_bed} -f {input.fai} \
        -o {output.rescaled_fai}
        """

    #####################
    ### Stairway Plot ###
    #####################

rule sfs:
    """
        Run easySFS with previously estimated sample size
    """
    input:
        vcf = "results/vcf/geneD/{prefix}.{density_cat}.rdmSNP.{chr}.vcf",
        pop_path = "results/{prefix}.pop",
        fai = "results/stats/geneD/{prefix}.{density_cat}.rdmSNP.trimmed.{chr}.fai", # for chr length
        easySFS = config["easySFS_path"]
    output:
        sfs_dir = temp(directory("results/sfs/geneD/{density_cat}/{prefix}.{density_cat}.{chr}")),
        final_sfs = "results/sfs/geneD/{density_cat}/{prefix}.{density_cat}.{chr}.sfs"
    conda:
        "../envs/easySFS.yml"
    shell:
        """
            sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 )))
            length=$(cut -f2 {input.fai})
            python3 {input.easySFS} -i {input.vcf} -p {input.pop_path} \
            --dtype 'int' \
            --proj $sampling_size \
            --ploidy 2 \
            --total-length $length \
            -v -a -f -o {output.sfs_dir}

            n_seq=$(sed -n '2 p' {output.sfs_dir}/fastsimcoal2/*MSFS.obs | \
                awk '{{print $2}}')
            echo $(sed -n '3 p' {output.sfs_dir}/fastsimcoal2/*MSFS.obs | \
                awk -v n_seq="$n_seq" '{{ for(i=2; i<=n_seq/2+1; i++) printf $i" "}}') \
                > {output.final_sfs}
        """


rule plot_sfs_small:
    """
        Plot SFS 
    """
    input:
        sfs = "results/sfs/geneD/{density_cat}/{prefix}.{density_cat}.{chr}.sfs",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/plot_sfs.py"
    output:
        sfs_plot = "results/sfs/geneD/{density_cat}/{prefix}.{density_cat}.{chr}.png"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
            python3 {input.script} -i {input.sfs} -o {output.sfs_plot}
        """

rule prepare_strway_plot_small:
    """
        Writes inputs for StairwayPlot
    """
    input:
        sfs = "results/sfs/geneD/{density_cat}/{prefix}.{density_cat}.{chr}.sfs",
        fai = "results/stats/geneD/{prefix}.{density_cat}.rdmSNP.trimmed.{chr}.fai", # for chr length
        stairway_plot_dir = config["stairway_plot_dir"],
    output:
        blueprint = "results/ne_inference/strway_plt/geneD/{density_cat}/{prefix}.{density_cat}.{chr}.blueprint"
    shell:
        """
            n_seq=$(( $(wc -w < {input.sfs}) * 2))
            total_sites=$(cut -f2 {input.fai})
            echo "# Settings for {wildcards.prefix}; chr : {wildcards.chr}; simple filtering on biallelic snps
popid: {wildcards.prefix}
nseq: $n_seq
L: $total_sites
whether_folded: true
SFS: $(cat {input.sfs})
smallest_size_of_SFS_bin_used_for_estimation: 1
largest_size_of_SFS_bin_used_for_estimation: $((n_seq/2 - 1))
pct_training: 0.67
nrand: $(((n_seq-2)/4))	$(((n_seq-2)/2))	$(((n_seq-2)*3/4))	$((n_seq-2))
project_dir: $(pwd)/results/stairway_plot/snps/small/{wildcards.prefix}.SNPS.small.{wildcards.chr}/
stairway_plot_dir: {input.stairway_plot_dir}
ninput: 200
#random_seed: 6
mu: {mu}
year_per_generation: 1
plot_title: {wildcards.prefix}.SNPS.small.{wildcards.chr}
xrange: 0.1,10000
yrange: 0,0
xspacing: 2
yspacing: 2
fontsize: 12" > {output.blueprint}
        """