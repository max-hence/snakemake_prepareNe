
mu = config["mutation_rate"]
generation = config["generation_time"]
    #################################
    ### For filter on Callability ###
    #################################

rule sfs_small:
    """
        Run easySFS with 10 samples
    """
    input:
        vcf = "results/callability/vcf/{prefix}.SNPS.NA.{chr}.vcf",
        pop_path = "results/{prefix}.pop",
        easySFS = config["easySFS_path"]
    output:
        sfs_dir = temp(directory("results/callability/sfs/small/{prefix}.small.{chr}")),
        final_sfs = "results/callability/sfs/small/{prefix}.small.{chr}.sfs"
    conda:
        "../envs/callability.yml"
    log:
        "logs/{prefix}.{chr}.log"
    shell:
        """
            sampling_size=20
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
        Run easySFS with all samples
    """
    input:
        vcf = "results/callability/vcf/{prefix}.SNPS.NA.{chr}.vcf",
        pop_path = "results/{prefix}.pop",
        easySFS = config["easySFS_path"],
    output:
        sfs_dir = temp(directory("results/callability/sfs/strict/{prefix}.strict.{chr}")),
        final_sfs = "results/callability/sfs/strict/{prefix}.strict.{chr}.sfs"
    conda:
        "../envs/callability.yml"
    log:
        "logs/{prefix}.{chr}.log"
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
rule sfs_ml:
    """
        Run easySFS with previously calculated sample size
    """
    input:
        vcf = "results/callability/vcf/{prefix}.SNPS.NA.{chr}.vcf",
        best_sample = "results/callability/sfs/{prefix}.SNPS.NA.best_sample.txt",
        pop_path = "results/{prefix}.pop",
        easySFS = config["easySFS_path"]
    output:
        sfs_dir = temp(directory("results/callability/sfs/{prefix}.ml.{chr}")),
        final_sfs = "results/callability/sfs/ml/{prefix}.ml.{chr}.sfs"
    conda:
        "../envs/callability.yml"
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
    input:
        sfs = "results/callability/sfs/{call_filter}/{prefix}.{call_filter}.{chr}.sfs",
        fai = "results/callability/stats/{prefix}.SNPS.NA.{call_filter}.{chr}.fai", # for chr length
        stairway_plot_dir = config["stairway_plot_dir"]
    output:
        blueprint = "results/callability/ne_inference/strway_plt/{call_filter}/{prefix}.SNPS.NA.{call_filter}.{chr}.blueprint"
    log:
        "logs/{prefix}.{call_filter}.{chr}.log"
    shell:
        """
        n_seq=$(( $(wc -w < {input.sfs}) * 2))
        total_sites=$(cut -f2 {input.fai})
        echo $n_seq >&2
        echo $total_sites >&2
        echo $n_seq
        echo $total_sites
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
project_dir: $(pwd)/results/callability/ne_inference/stairway_plt/{wildcards.prefix}.SNPS.NA.{wildcards.call_filter}.{wildcards.chr}/
stairway_plot_dir: {input.stairway_plot_dir}
ninput: 100
#random_seed: 6
mu: {mu}
year_per_generation: {generation}
plot_title: {wildcards.prefix}.SNPS.NA.{wildcards.call_filter}.{wildcards.chr}
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
        Reverse bed file to show position that will be excluded
    """
    input:
        bed = "results/callability/bed/{call_filter}/{prefix}.SNPS.NA.{call_filter}.{chr}.callable.bed",
        raw_bed = "results/raw/bed/{prefix}.raw.{chr}.callable.bed",
        script = workflow.source_path("../scripts/reverse_bed.py")
    output:
        mask = temp("tmp/{prefix}.SNPS.NA.{chr}.{call_filter}.callable.flipped.bed"),
        mask_gz = "results/callability/bed/{call_filter}/{prefix}.SNPS.NA.{call_filter}.{chr}.callable.flipped.bed.gz",
        mask_idx = "results/callability/bed/{call_filter}/{prefix}.SNPS.NA.{call_filter}.{chr}.callable.flipped.bed.gz.tbi"
    conda:
        "../envs/callability.yml"
    log:
        "logs/{prefix}.{call_filter}.{chr}.log"
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
        vcf = "results/callability/vcf/{prefix}.SNPS.NA.{chr}.vcf.gz",
        vcf_idx = "results/callability/vcf/{prefix}.SNPS.NA.{chr}.vcf.gz.tbi",
        pop_path = "results/{prefix}.pop",
        fai = "results/callability/stats/{prefix}.SNPS.NA.{call_filter}.{chr}.fai",
        mask = "results/callability/bed/{call_filter}/{prefix}.SNPS.NA.{call_filter}.{chr}.callable.flipped.bed.gz",
        mask_idx = "results/callability/bed/{call_filter}/{prefix}.SNPS.NA.{call_filter}.{chr}.callable.flipped.bed.gz.tbi",
        smcpp = config["smcpp"]
    output:
        smcpp_input = "results/callability/ne_inference/smcpp/{call_filter}/{prefix}.SNPS.NA.{call_filter}.{chr}.smc.gz"
    conda:
        "../envs/callability.yml"
    log:
        "logs/{prefix}.{call_filter}.{chr}.log"
    shell:
        """
        length=$(cut -f2 {input.fai})
        pop_input=$(awk '{{a[$2] = a[$2] (a[$2] ? "," : "") $1}} END {{ for (key in a) {{print key ":" a[key]}} }}' {input.pop_path}) 

        if [[ $(zcat {input.mask} | wc -c) -eq 0 ]]; then
            singularity exec {input.smcpp} smc++ vcf2smc \
                {input.vcf} \
                {output.smcpp_input} \
                {wildcards.chr} \
                $pop_input \
                --length $length \
                -v
        else
            singularity exec {input.smcpp} smc++ vcf2smc \
                {input.vcf} \
                {output.smcpp_input} \
                {wildcards.chr} \
                $pop_input \
                --length $length \
                --mask {input.mask} \
                -v
        fi
        """