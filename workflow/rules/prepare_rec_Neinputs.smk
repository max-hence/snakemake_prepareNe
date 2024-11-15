# Prepare inputs files for stairway_plot and smc++

    #####################
    ### Stairway Plot ###
    #####################

rule sfs_rec:
    """
        Run easySFS with best param values
    """
    input:
        vcf = "results/vcf/rec/{prefix}.{rec_quant}.320kSNP.{chr}.vcf",
        pop_path = "results/{prefix}.pop",
        best_sample = "results/sfs/rec/{prefix}.{rec_quant}.best_sample.txt",
        fai = "results/stats/rec/ml/{prefix}.{rec_quant}.320kSNP.{chr}.fai", # for chr length
        easySFS = config["easySFS_path"]
    output:
        sfs_dir = temp(directory("results/sfs/rec/{rec_quant}/{prefix}.{rec_quant}.{chr}")),
        final_sfs = "results/sfs/rec/{rec_quant}/{prefix}.{rec_quant}.{chr}.sfs"
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

rule plot_sfs_rec:
    """
        Plot SFS
    """
    input:
        sfs = "results/sfs/rec/{rec_quant}/{prefix}.{rec_quant}.{chr}.sfs",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/plot_sfs.py"
    output:
        sfs_plot = "results/sfs/rec/{rec_quant}/{prefix}.{rec_quant}.{chr}.png"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
            python3 {input.script} -i {input.sfs} -o {output.sfs_plot}
        """

rule prepare_strway_plot_rec:
    """
        Writes inputs for StairwayPlot
    """
    input:
        sfs = "results/sfs/rec/{rec_quant}/{prefix}.{rec_quant}.{chr}.sfs",
        fai = "results/stats/rec/ml/{prefix}.{rec_quant}.320kSNP.{chr}.fai", # for chr length
        stairway_plot_dir = config["stairway_plot_dir"],
    output:
        blueprint = "results/ne_inference/strway_plt/rec/{rec_quant}/{prefix}.{rec_quant}.{chr}.blueprint"
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
project_dir: $(pwd)/results/stairway_plot/rec/{wildcards.rec_quant}/{wildcards.prefix}.{wildcards.rec_quant}.{wildcards.chr}/
stairway_plot_dir: {input.stairway_plot_dir}
ninput: 200
#random_seed: 6
mu: {mu}
year_per_generation: 1
plot_title: {wildcards.prefix}.{wildcards.rec_quant}.{wildcards.chr}
xrange: 0.1,10000
yrange: 0,0
xspacing: 2
yspacing: 2
fontsize: 12" > {output.blueprint}
        """


    #############
    ### SMC++ ###
    #############


rule flip_bed_rec:
    """
        Reverse bed file to show positions that will be excluded for SMC++ input
    """
    input:
        bed = "results/bed/rec/ml/{prefix}.{rec_quant}.{chr}.callable.bed",
        raw_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/reverse_bed.py"
    output:
        mask = temp("results/bed/rec/{prefix}.{rec_quant}.callable.flipped.{chr}.bed"),
        mask_gz = "results/bed/rec/{prefix}.{rec_quant}.callable.flipped.{chr}.bed.gz",
        mask_idx = "results/bed/rec/{prefix}.{rec_quant}.callable.flipped.{chr}.bed.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        last_pos=$(tail -1 {input.raw_bed} | cut -f3)
        python3 {input.script} -i {input.bed} -l $last_pos -o {output.mask}

        bgzip -c {output.mask} > {output.mask_gz}
        tabix -p bed {output.mask_gz}
        """

rule prepare_smcpp_rec:
    """
        Prepare .smc file per chr with the callability.bed as mask
    """
    input:
        vcf = "results/vcf/rec/{prefix}.{rec_quant}.320kSNP.{chr}.vcf.gz",
        vcf_idx = "results/vcf/rec/{prefix}.{rec_quant}.320kSNP.{chr}.vcf.gz.tbi",
        pop_path = "results/{prefix}.pop",
        fai = "results/stats/rec/ml/{prefix}.{rec_quant}.320kSNP.{chr}.fai",
        mask = "results/bed/rec/{prefix}.{rec_quant}.callable.flipped.{chr}.bed.gz",
        mask_idx = "results/bed/rec/{prefix}.{rec_quant}.callable.flipped.{chr}.bed.gz.tbi",
        smcpp = config["smcpp"]
    output:
        smcpp_input = "results/ne_inference/smcpp/rec/{rec_quant}/{prefix}.{rec_quant}.{chr}.smc.gz"
    conda:
        "../envs/vcf_processing.yml"
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