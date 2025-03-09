# Prepare inputs for Stairway Plot
# Per chr analysis

include: "common.smk"
mu = config["mutation_rate"]
generation = config["generation_time"]

rule sfs:
    """
        easySFS with previously calculated sample size
    """
    input:
        vcf = "results/snps/vcf/{prefix}.SNPS.{chr}.vcf",
        best_sample = "results/snps/sfs/{prefix}.SNPS.best_sample.txt",
        pop_path = "results/{prefix}.pop",
        easySFS = config["easySFS_path"]
    output:
        sfs_dir = temp(directory("results/snps/sfs/{prefix}.SNPS.{chr}")),
        final_sfs = "results/snps/sfs/{prefix}.{chr}.sfs"
    conda:
        "../envs/easySFS.yml"
    log:
        "logs/{prefix}.{chr}.log"
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

rule prepare_strway_plot:
    """
        Writes inputs for StairwayPlot
    """
    input:
        sfs = "results/snps/sfs/{prefix}.{chr}.sfs",
        fai = "results/snps/stats/{prefix}.SNPS.resized.{chr}.fai", # for chr length
        stairway_plot_dir = config["stairway_plot_dir"]
    output:
        blueprint = "results/snps/ne_inference/strway_plt/{prefix}.SNPS.{chr}.blueprint"
    log:
        "logs/{prefix}.{chr}.log"
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
project_dir: $(pwd)/results/ne_inference/stairway_plot/{wildcards.prefix}.SNPS.{wildcards.chr}/
stairway_plot_dir: {input.stairway_plot_dir}
ninput: 100
#random_seed: 6
mu: {mu}
year_per_generation: {generation}
plot_title: {wildcards.prefix}.SNPS.{wildcards.chr}
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
        Reverse bed file to show positions that will be excluded
    """
    input:
        bed = "results/snps/bed/{prefix}.SNPS.{chr}.callable.bed",
        raw_bed = "results/raw/bed/{prefix}.raw.{chr}.callable.bed",
        script = workflow.source_path("../scripts/reverse_bed.py")
    output:
        mask = temp("results/snps/bed/{prefix}.SNPS.callable.flipped.{chr}.bed"),
        mask_gz = "results/snps/bed/{prefix}.SNPS.callable.flipped.{chr}.bed.gz",
        mask_idx = "results/snps/bed/{prefix}.SNPS.callable.flipped.{chr}.bed.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.{chr}.log"
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
        vcf = "results/snps/vcf/{prefix}.SNPS.{chr}.vcf.gz",
        vcf_idx = "results/snps/vcf/{prefix}.SNPS.{chr}.vcf.gz.tbi",
        pop_path = "results/{prefix}.pop",
        fai = "results/snps/stats/{prefix}.SNPS.resized.{chr}.fai",
        mask = "results/snps/bed/{prefix}.SNPS.callable.flipped.{chr}.bed.gz",
        mask_idx = "results/snps/bed/{prefix}.SNPS.callable.flipped.{chr}.bed.gz.tbi",
        smcpp = config["smcpp"]
    output:
        smcpp_input = "results/snps/ne_inference/smcpp/{prefix}.SNPS.{chr}.smc.gz"
    conda:
        "../envs/vcf_processing.yml"
    log:
        "logs/{prefix}.{chr}.log"
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