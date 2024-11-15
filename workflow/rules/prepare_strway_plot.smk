# Prepare inputs for Stairway Plot
# Per chr analysis

include: "common.smk"
mu = config["mutation_rate"]


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
        sfs_dir = temp(directory("results/sfs/snps/small/{prefix}.small.{chr}")),
        final_sfs = "results/sfs/snps/small/{prefix}.small.{chr}.sfs"
    conda:
        "../envs/easySFS.yml"
    shell:
        """
            sampling_size=20
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
        sfs = "results/sfs/snps/small/{prefix}.small.{chr}.sfs",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/plot_sfs.py"
    output:
        sfs_plot = "results/sfs/snps/small/{prefix}.small.{chr}.png"
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
        sfs = "results/sfs/snps/small/{prefix}.small.{chr}.sfs",
        fai = "results/stats/snps/small/{prefix}.SNPS.small.{chr}.fai", # for chr length
        stairway_plot_dir = config["stairway_plot_dir"],
    output:
        blueprint = "results/ne_inference/strway_plt/snps/small/{prefix}.SNPS.small.{chr}.blueprint"
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

    ### Resampled vcf (max snps)
rule sfs_max:
    """
        Run easySFS with previously calculated sample size
    """
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf",
        best_sample = "results/sfs/snps/{prefix}.SNPS.best_sample_max.txt",
        pop_path = "results/{prefix}.pop",
        fai = "results/stats/snps/max/{prefix}.SNPS.max.{chr}.fai", # for chr length
        easySFS = config["easySFS_path"]
    output:
        sfs_dir = temp(directory("results/sfs/snps/max/{prefix}.max.{chr}")),
        final_sfs = "results/sfs/snps/max/{prefix}.max.{chr}.sfs"
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

rule plot_sfs_max:
    """
        Plot SFS
    """
    input:
        sfs = "results/sfs/snps/max/{prefix}.max.{chr}.sfs",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/plot_sfs.py"
    output:
        sfs_plot = "results/sfs/snps/max/{prefix}.max.{chr}.png"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
            python3 {input.script} -i {input.sfs} -o {output.sfs_plot}
        """

rule prepare_strway_plot_max:
    """
        Writes inputs for StairwayPlot
    """
    input:
        sfs = "results/sfs/snps/max/{prefix}.max.{chr}.sfs",
        fai = "results/stats/snps/max/{prefix}.SNPS.max.{chr}.fai", # for chr length
        stairway_plot_dir = config["stairway_plot_dir"],
    output:
        blueprint = "results/ne_inference/strway_plt/snps/max/{prefix}.SNPS.max.{chr}.blueprint"
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
project_dir: $(pwd)/results/stairway_plot/snps/max/{wildcards.prefix}.SNPS.max.{wildcards.chr}/
stairway_plot_dir: {input.stairway_plot_dir}
ninput: 200
#random_seed: 6
mu: {mu}
year_per_generation: 1
plot_title: {wildcards.prefix}.SNPS.max.{wildcards.chr}
xrange: 0.1,10000
yrange: 0,0
xspacing: 2
yspacing: 2
fontsize: 12" > {output.blueprint}
        """

    ### Resampled vcf (ml method)
rule sfs_ml:
    """
        Run easySFS with previously calculated sample size
    """
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf",
        best_sample = "results/sfs/snps/{prefix}.SNPS.best_sample_ml.txt",
        pop_path = "results/{prefix}.pop",
        fai = "results/stats/snps/ml/{prefix}.SNPS.ml.{chr}.fai", # for chr length
        easySFS = config["easySFS_path"]
    output:
        sfs_dir = temp(directory("results/sfs/snps/ml/{prefix}.ml.{chr}")),
        final_sfs = "results/sfs/snps/ml/{prefix}.ml.{chr}.sfs"
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


rule plot_sfs_ml:
    """
        Plot SFS 
    """
    input:
        sfs = "results/sfs/snps/ml/{prefix}.ml.{chr}.sfs",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/plot_sfs.py"
    output:
        sfs_plot = "results/sfs/snps/ml/{prefix}.ml.{chr}.png"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
            python3 {input.script} -i {input.sfs} -o {output.sfs_plot}
        """

rule prepare_strway_plot_ml:
    """
        Writes inputs for StairwayPlot
    """
    input:
        sfs = "results/sfs/snps/ml/{prefix}.ml.{chr}.sfs",
        fai = "results/stats/snps/ml/{prefix}.SNPS.ml.{chr}.fai", # for chr length
        stairway_plot_dir = config["stairway_plot_dir"],
    output:
        blueprint = "results/ne_inference/strway_plt/snps/ml/{prefix}.SNPS.ml.{chr}.blueprint"
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
project_dir: $(pwd)/results/stairway_plot/snps/ml/{wildcards.prefix}.SNPS.ml.{wildcards.chr}/
stairway_plot_dir: {input.stairway_plot_dir}
ninput: 200
#random_seed: 6
mu: {mu}
year_per_generation: 1
plot_title: {wildcards.prefix}.SNPS.ml.{wildcards.chr}
xrange: 0.1,10000
yrange: 0,0
xspacing: 2
yspacing: 2
fontsize: 12" > {output.blueprint}
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

rule plot_sfs_strict:
    """
        Plot SFS
    """
    input:
        sfs = "results/sfs/snps/strict/{prefix}.strict.{chr}.sfs",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/plot_sfs.py"
    output:
        sfs_plot = "results/sfs/snps/strict/{prefix}.strict.{chr}.png"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
            python3 {input.script} -i {input.sfs} -o {output.sfs_plot}
        """

rule prepare_strway_plot_strict:
    """
        Writes inputs for StairwayPlot
    """
    input:
        sfs = "results/sfs/snps/strict/{prefix}.strict.{chr}.sfs",
        fai = "results/stats/snps/strict/{prefix}.SNPS.strict.{chr}.fai", # for chr length
        stairway_plot_dir = config["stairway_plot_dir"],
    output:
        blueprint = "results/ne_inference/strway_plt/snps/strict/{prefix}.SNPS.strict.{chr}.blueprint"
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
project_dir: $(pwd)/results/stairway_plot/snps/strict/{wildcards.prefix}.SNPS.strict.{wildcards.chr}/
stairway_plot_dir: {input.stairway_plot_dir}
ninput: 200
#random_seed: 6
mu: {mu}
year_per_generation: 1
plot_title: {wildcards.prefix}.SNPS.strict.{wildcards.chr}
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

rule plot_sfs_na:
    """
        Plot SFS
    """
    input:
        sfs = "results/sfs/snps_na/{prefix}.SNPS.NA.{chr}.sfs",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/plot_sfs.py"
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
project_dir: $(pwd)/results/stairway_plot/SNPS.NA/{wildcards.prefix}.SNPS.NA.{wildcards.chr}/
stairway_plot_dir: {input.stairway_plot_dir}
ninput: 200
#random_seed: 6
mu: {mu}
year_per_generation: 1
plot_title: {wildcards.prefix}.SNPS.NA.{wildcards.chr}
xrange: 0.1,10000
yrange: 0,0
xspacing: 2
yspacing: 2
fontsize: 12" > {output.blueprint}
        """