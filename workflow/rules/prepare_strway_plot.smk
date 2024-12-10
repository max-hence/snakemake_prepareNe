# Prepare inputs for Stairway Plot
# Per chr analysis

include: "common.smk"
mu = config["mutation_rate"]
generation = config["generation_time"]


    #####################################
    ### For filter on Bi-allelic SNPs ###
    #####################################

    ### Resampled vcf (10 indiv)
rule sfs_small:
    """
        Run easySFS with previously calculated sample size
    """
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf",
        pop_path = "results/{prefix}.pop",
        fai = "results/stats/snps/small/{prefix}.SNPS.small.{chr}.fai", # for chr length
        easySFS = config["easySFS_path"]
    output:
        sfs_dir = directory("results/sfs/snps/small/{prefix}.small.{chr}"),
        final_sfs = "results/sfs/snps/small/{prefix}.small.{chr}.sfs"
    conda:
        "../envs/easySFS.yml"
    shell:
        """
            sampling_size=10
            python3 {input.easySFS} -i {input.vcf} -p {input.pop_path} \
            --dtype 'int' \
            --proj $sampling_size \
            --ploidy 2 \
            -v -a -f -o {output.sfs_dir}

            n_seq=$(sed -n '2 p' {output.sfs_dir}/fastsimcoal2/*MSFS.obs | \
                awk '{{print $2}}')
            echo $(sed -n '3 p' {output.sfs_dir}/fastsimcoal2/*MSFS.obs | \
                awk -v n_seq="$n_seq" '{{ for(i=2; i<=n_seq/2+1; i++) printf $i" "}}') \
                > {output.final_sfs}
        """

    ### Resampled vcf (all indiv)
rule sfs_strict:
    """
        Run easySFS with previously calculated sample size
    """
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf",
        pop_path = "results/{prefix}.pop",
        fai = "results/stats/snps/strict/{prefix}.SNPS.strict.{chr}.fai", # for chr length
        easySFS = config["easySFS_path"],
    output:
        sfs_dir = temp(directory("results/sfs/snps/strict/{prefix}.strict.{chr}")),
        final_sfs = "results/sfs/snps/strict/{prefix}.strict.{chr}.sfs"
    conda:
        "../envs/easySFS.yml"
    shell:
        """
            sampling_size=$(( $(wc -l < {input.pop_path}) * 2 ))
            python3 {input.easySFS} -i {input.vcf} -p {input.pop_path} \
            --dtype 'int' \
            --proj $sampling_size \
            --ploidy 2 \
            -v -a -f -o {output.sfs_dir}

            n_seq=$(sed -n '2 p' {output.sfs_dir}/fastsimcoal2/*MSFS.obs | \
                awk '{{print $2}}')
            echo $(sed -n '3 p' {output.sfs_dir}/fastsimcoal2/*MSFS.obs | \
                awk -v n_seq="$n_seq" '{{ for(i=2; i<=n_seq/2+1; i++) printf $i" "}}') \
                > {output.final_sfs}
        """

    ### Resampled vcf (max snps)
rule sfs_max_ml:
    """
        Run easySFS with previously calculated sample size
    """
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf",
        best_sample = "results/sfs/snps/{prefix}.SNPS.best_sample_{sfs_params_method}.txt",
        pop_path = "results/{prefix}.pop",
        fai = "results/stats/snps/{sfs_params_method}/{prefix}.SNPS.{sfs_params_method}.{chr}.fai", # for chr length
        easySFS = config["easySFS_path"]
    output:
        sfs_dir = temp(directory("results/sfs/snps/{sfs_params_method}/{prefix}.{sfs_params_method}.{chr}")),
        final_sfs = "results/sfs/snps/{sfs_params_method}/{prefix}.{sfs_params_method}.{chr}.sfs"
    conda:
        "../envs/easySFS.yml"
    shell:
        """
            sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 )))
            python3 {input.easySFS} -i {input.vcf} -p {input.pop_path} \
            --dtype 'int' \
            --proj $sampling_size \
            --ploidy 2 \
            -v -a -f -o {output.sfs_dir}

            n_seq=$(sed -n '2 p' {output.sfs_dir}/fastsimcoal2/*MSFS.obs | \
                awk '{{print $2}}')
            echo $(sed -n '3 p' {output.sfs_dir}/fastsimcoal2/*MSFS.obs | \
                awk -v n_seq="$n_seq" '{{ for(i=2; i<=n_seq/2+1; i++) printf $i" "}}') \
                > {output.final_sfs}
        """

rule plot_sfs:
    """
        Plot SFS
    """
    input:
        sfs = "results/sfs/snps/{subsample}/{prefix}.{subsample}.{chr}.sfs",
        script = workflow.source_path("../scripts/plot_sfs.py")
    output:
        sfs_plot = "results/sfs/snps/{subsample}/{prefix}.{subsample}.{chr}.png"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
            python3 {input.script} -i {input.sfs} -o {output.sfs_plot}
        """


rule prepare_strway_plot:
    """
        Writes inputs for StairwayPlot
    """
    input:
        sfs = "results/sfs/snps/{subsample}/{prefix}.{subsample}.{chr}.sfs",
        fai = "results/stats/snps/{subsample}/{prefix}.SNPS.{subsample}.{chr}.fai", # for chr length
        stairway_plot_dir = config["stairway_plot_dir"],
    output:
        blueprint = "results/ne_inference/strway_plt/{subsample}/{prefix}.SNPS.{subsample}.{chr}.blueprint"
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
project_dir: $(pwd)/results/stairway_plot/{wildcards.subsample}/{wildcards.prefix}.SNPS.{wildcards.subsample}.{wildcards.chr}/
stairway_plot_dir: {input.stairway_plot_dir}
ninput: 200
#random_seed: 6
mu: {mu}
year_per_generation: {generation}
plot_title: {wildcards.prefix}.SNPS.{wildcards.subsample}.{wildcards.chr}
xrange: 0.1,10000
yrange: 0,0
xspacing: 2
yspacing: 2
fontsize: 12" > {output.blueprint}
        """

    #################################
    ### For filter on Callability ###
    #################################

rule sfs_na:
    """
        Run easySFS with previously calculated sample size
    """
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz",
        best_sample = "results/sfs/snps_na/{prefix}.SNPS.NA.best_sample.txt",
        pop_path = "results/{prefix}.pop",
        fai = "results/stats/snps_na/{prefix}.SNPS.NA.{chr}.fai", # for chr length
        easySFS = config["easySFS_path"],
    output:
        sfs_dir = directory("results/sfs/snps_na/{prefix}.{chr}"),
        final_sfs = "results/sfs/snps_na/{prefix}.SNPS.NA.{chr}.sfs"
    conda:
        "../envs/easySFS.yml"
    shell:
        """
            sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 )))
            python3 {input.easySFS} -i {input.vcf} -p {input.pop_path} \
            --dtype 'int' \
            --proj $sampling_size \
            --ploidy 2 \
            -v -a -f -o {output.sfs_dir}

            n_seq=$(sed -n '2 p' {output.sfs_dir}/fastsimcoal2/*MSFS.obs | \
                awk '{{print $2}}')
            echo $(sed -n '3 p' {output.sfs_dir}/fastsimcoal2/*MSFS.obs | \
                awk -v n_seq="$n_seq" '{{ for(i=2; i<=n_seq/2+1; i++) printf $i" "}}') \
                > {output.final_sfs}
        """

rule plot_sfs_na:
    """
        Plot SFS
    """
    input:
        sfs = "results/sfs/snps_na/{prefix}.SNPS.NA.{chr}.sfs",
        script = workflow.source_path("../scripts/plot_sfs.py")
    output:
        sfs_plot = "results/sfs/snps_na/{prefix}.{chr}.png"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
            python3 {input.script} -i {input.sfs} -o {output.sfs_plot}
        """

rule prepare_strway_plot_na:
    input:
        sfs = "results/sfs/snps_na/{prefix}.SNPS.NA.{chr}.sfs",
        fai = "results/stats/snps_na/{prefix}.SNPS.NA.{chr}.fai", # for chr length
        stairway_plot_dir = config["stairway_plot_dir"]
    output:
        blueprint = "results/ne_inference/strway_plt/snps_na/{prefix}.SNPS.NA.{chr}.blueprint"
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
project_dir: $(pwd)/results/stairway_plot/snps_na/{wildcards.prefix}.SNPS.NA.{wildcards.chr}/
stairway_plot_dir: {input.stairway_plot_dir}
ninput: 200
#random_seed: 6
mu: {mu}
year_per_generation: {generation}
plot_title: {wildcards.prefix}.SNPS.NA.{wildcards.chr}
xrange: 0.1,10000
yrange: 0,0
xspacing: 2
yspacing: 2
fontsize: 12" > {output.blueprint}
        """