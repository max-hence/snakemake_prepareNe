# Prepare inputs files for stairway_plot and smc++ for vcf splitted by gene density

include: "common.smk"
mu = config["mutation_rate"]
generation = config["generation_time"]

rule sfs_projection:
    """
        Run easySFS to do SFS projection and find the best sample/snps ratio
    """
        input:
            vcf = "results/geneD/vcf/{prefix}.{density}.{chr}.vcf",
            pop_path = "results/{prefix}.pop",
            easySFS = config["easySFS_path"],
        output:
            preview = "results/geneD/sfs/{prefix}.{density}.preview.{chr}.txt",
        conda:
            "../envs/gene_density.yml"
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
        previews = get_previews,
        get_sfs_param = workflow.source_path("../scripts/get_sfs_param.py")
    output:
        best_sample = "results/geneD/sfs/{prefix}.best_params.txt",
    conda:
        "../envs/gene_density.yml"
    shell:
        """
        python3 {input.get_sfs_param} -i {input.previews} -m "ml" -o {output.best_sample}
        """

rule intersect_beds:
    """
    Trim callable bed to keep only regions inside the given rec threshold (merging btw 2 beds)
    """
    input:
        bed = "results/geneD/bed/{prefix}.{density}.{chr}.bed",
        callable_bed = "results/raw/bed/{prefix}.raw.{chr}.callable.bed",
        script = workflow.source_path("../scripts/merge_bed.py")
    output:
        intersect_bed = "results/geneD/bed/{prefix}.{density}.intersect.{chr}.bed"
    conda:
        "../envs/gene_density.yml"
    log:
        "logs/rec/{prefix}.{density}.{chr}.log"
    shell:
        """
        python3 {input.script} -i {input.bed} -b {input.callable_bed} -o {output.intersect_bed}
        """

rule trim_bed:
    """
        Remove regions in callability.bed where less than <best_sample> indiv have been well called
        Resize chr length based on trimmed bed
    """
    input:
        best_sample = "results/geneD/sfs/{prefix}.best_params.txt",
        fai = "results/geneD/stats/{prefix}.{density}.{chr}.fai",
        intersect_bed = "results/geneD/bed/{prefix}.{density}.intersect.{chr}.bed",
        script = workflow.source_path("../scripts/rescale_genlen.py")
    output:
        trimmed_bed = "results/geneD/bed/{prefix}.{density}.{chr}.callable.bed",
        rescaled_fai = "results/geneD/stats/{prefix}.{density}.resized.{chr}.fai"
    conda:
        "../envs/gene_density.yml"
    shell:
        """
        sampling_size=$(( $(tail -1 {input.best_sample} | cut -f1 ) / 2 ))
        awk -v n=$sampling_size '$4 >= n' {input.intersect_bed} > {output.trimmed_bed}
        
        python3 {input.script} -i {output.trimmed_bed} -f {input.fai} \
        -o {output.rescaled_fai} --method bed
        """

    #####################
    ### Stairway Plot ###
    #####################

rule sfs:
    """
        Run easySFS with previously estimated sample size
    """
    input:
        vcf = "results/geneD/vcf/{prefix}.{density}.{chr}.vcf",
        pop_path = "results/{prefix}.pop",
        fai = "results/geneD/stats/{prefix}.{density}.resized.{chr}.fai", # for chr length
        best_sample = "results/geneD/sfs/{prefix}.best_params.txt",
        easySFS = config["easySFS_path"]
    output:
        sfs_dir = temp(directory("results/geneD/sfs/{density}/{prefix}.{density}.{chr}")),
        final_sfs = "results/geneD/sfs/{density}/{prefix}.{density}.{chr}.sfs"
    conda:
        "../envs/gene_density.yml"
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


rule plot_sfs:
    """
        Plot SFS 
    """
    input:
        sfs = "results/geneD/sfs/{density}/{prefix}.{density}.{chr}.sfs",
        script = workflow.source_path("../scripts/plot_sfs.py")
    output:
        sfs_plot = "results/geneD/sfs/{density}/{prefix}.{density}.{chr}.png"
    conda:
        "../envs/gene_density.yml"
    shell:
        """
            python3 {input.script} -i {input.sfs} -o {output.sfs_plot}
        """

rule prepare_strway_plot:
    """
        Writes inputs for StairwayPlot
    """
    input:
        sfs = "results/geneD/sfs/{density}/{prefix}.{density}.{chr}.sfs",
        fai = "results/geneD/stats/{prefix}.{density}.resized.{chr}.fai", # for chr length
        stairway_plot_dir = config["stairway_plot_dir"],
    output:
        blueprint = "results/geneD/ne_inference/strway_plt/{density}/{prefix}.{density}.{chr}.blueprint"
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
project_dir: $(pwd)/results/geneD/stairway_plot/{wildcards.prefix}.{wildcards.density}.{wildcards.chr}/
stairway_plot_dir: {input.stairway_plot_dir}
ninput: 100
#random_seed: 6
mu: {mu}
year_per_generation: {generation}
plot_title: {wildcards.prefix}.{wildcards.density}.{wildcards.chr}
xrange: 0.1,10000
yrange: 0,0
xspacing: 2
yspacing: 2
fontsize: 12" > {output.blueprint}
        """

    ##############
    ### SMC ++ ###
    ##############

rule flip_bed:
    """
        Reverse bed file to show positions that will be excluded for SMC++ input
    """
    input:
        bed = "results/geneD/bed/{prefix}.{density}.{chr}.callable.bed",
        raw_bed = "results/raw/bed/{prefix}.raw.{chr}.callable.bed",
        script = workflow.source_path("../scripts/reverse_bed.py")
    output:
        mask = temp("results/geneD/bed/{prefix}.{density}.callable.flipped.{chr}.bed"),
        mask_gz = "results/geneD/bed/{prefix}.{density}.callable.flipped.{chr}.bed.gz",
        mask_idx = "results/geneD/bed/{prefix}.{density}.callable.flipped.{chr}.bed.gz.tbi"
    conda:
        "../envs/gene_density.yml"
    shell:
        """
        last_pos=$(tail -1 {input.raw_bed} | cut -f3)
        python3 {input.script} -i {input.bed} -l $last_pos -o {output.mask}

        bgzip -c {output.mask} > {output.mask_gz}
        tabix -p bed {output.mask_gz}
        """

rule prepare_smcpp:
    """
        Prepare .smc file per chr with the callability.bed as mask
    """
    input:
        vcf = "results/geneD/vcf/{prefix}.{density}.{chr}.vcf.gz",
        vcf_idx = "results/geneD/vcf/{prefix}.{density}.{chr}.vcf.gz.tbi",
        pop_path = "results/{prefix}.pop",
        fai = "results/geneD/stats/{prefix}.{density}.resized.{chr}.fai",
        mask = "results/geneD/bed/{prefix}.{density}.callable.flipped.{chr}.bed.gz",
        mask_idx = "results/geneD/bed/{prefix}.{density}.callable.flipped.{chr}.bed.gz.tbi",
        smcpp = config["smcpp"]
    output:
        smcpp_input = "results/geneD/ne_inference/smcpp/{density}/{prefix}.{density}.{chr}.smc.gz"
    conda:
        "../envs/gene_density.yml"
    shell:
        """
        length=$(cut -f2 {input.fai})
        pop_input=$(awk '{{a[$2] = a[$2] (a[$2] ? "," : "") $1}} END {{ for (key in a) {{print key ":" a[key]}} }}' {input.pop_path}) 

        singularity exec {input.smcpp} smc++ vcf2smc \
            {input.vcf} \
            {output.smcpp_input} \
            {wildcards.chr} \
            $pop_input \
            --length $length \
            --mask {input.mask} \
            -v
        """