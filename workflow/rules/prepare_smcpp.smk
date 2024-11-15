# Correct genotypes based on callability bed file 
# Run easy_SFS to determine the best indiv number to take to maximize SNPs and minimize NA
# Then resample the bed file

include: "common.smk"

    #####################################
    ### For filter on Bi-allelic SNPs ###
    #####################################

    ### Full vcf

rule prepare_smcpp_full:
    """
        Prepare .smc file per chr with the callability.bed as mask
    """
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz",
        vcf_idx = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi",
        pop_path = "results/{prefix}.pop",
        fai = "results/stats/snps/full/{prefix}.SNPS.full.{chr}.fai",
        smcpp = config["smcpp"]
    output:
        smcpp_input = "results/ne_inference/smcpp/snps/full/{prefix}.SNPS.full.{chr}.smc.gz"
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
            -v
        """

    ### Resampled vcf (10 indiv)

rule flip_bed_small:
    """
        Reverse bed file to show positions that will be excluded
    """
    input:
        bed = "results/bed/snps/small/{prefix}.SNPS.small.{chr}.callable.bed",
        raw_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/reverse_bed.py"
    output:
        mask = temp("results/bed/{prefix}.SNPS.small.callable.flipped.{chr}.bed"),
        mask_gz = "results/bed/snps/small/{prefix}.SNPS.small.callable.flipped.{chr}.bed.gz",
        mask_idx = "results/bed/snps/small/{prefix}.SNPS.small.callable.flipped.{chr}.bed.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        last_pos=$(tail -1 {input.raw_bed} | cut -f3)
        python3 {input.script} -i {input.bed} -l $last_pos -o {output.mask}

        bgzip -c {output.mask} > {output.mask_gz}
        tabix -p bed {output.mask_gz}
        """

rule prepare_smcpp_small:
    """
        Prepare .smc file per chr with the callability.bed as mask
    """
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz",
        vcf_idx = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi",
        pop_path = "results/{prefix}.pop",
        fai = "results/stats/snps/full/{prefix}.SNPS.full.{chr}.fai",
        mask = "results/bed/snps/small/{prefix}.SNPS.small.callable.flipped.{chr}.bed.gz",
        mask_idx = "results/bed/snps/small/{prefix}.SNPS.small.callable.flipped.{chr}.bed.gz.tbi",
        smcpp = config["smcpp"]
    output:
        smcpp_input = "results/ne_inference/smcpp/snps/small/{prefix}.SNPS.small.{chr}.smc.gz"
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

    ### Resampled vcf (max snps)

rule flip_bed_max:
    """
        Reverse bed file to show positions that will be excluded
    """
    input:
        bed = "results/bed/snps/max/{prefix}.SNPS.max.{chr}.callable.bed",
        raw_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/reverse_bed.py"
    output:
        mask = temp("results/bed/{prefix}.SNPS.max.callable.flipped.{chr}.bed"),
        mask_gz = "results/bed/{prefix}.SNPS.max.callable.flipped.{chr}.bed.gz",
        mask_idx = "results/bed/{prefix}.SNPS.max.callable.flipped.{chr}.bed.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        last_pos=$(tail -1 {input.raw_bed} | cut -f3)
        python3 {input.script} -i {input.bed} -l $last_pos -o {output.mask}

        bgzip -c {output.mask} > {output.mask_gz}
        tabix -p bed {output.mask_gz}
        """

rule prepare_smcpp_max:
    """
        Prepare .smc file per chr with the callability.bed as mask
    """
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz",
        vcf_idx = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi",
        pop_path = "results/{prefix}.pop",
        fai = "results/stats/snps/max/{prefix}.SNPS.max.{chr}.fai",
        mask = "results/bed/{prefix}.SNPS.max.callable.flipped.{chr}.bed.gz",
        mask_idx = "results/bed/{prefix}.SNPS.max.callable.flipped.{chr}.bed.gz.tbi",
        smcpp = config["smcpp"]
    output:
        smcpp_input = "results/ne_inference/smcpp/snps/max/{prefix}.SNPS.max.{chr}.smc.gz"
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


    ### Resampled vcf (ml method)

rule flip_bed_ml:
    """
        Reverse bed file to show positions that will be excluded
    """
    input:
        bed = "results/bed/snps/ml/{prefix}.SNPS.ml.{chr}.callable.bed",
        raw_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/reverse_bed.py"
    output:
        mask = temp("results/bed/{prefix}.SNPS.ml.callable.flipped.{chr}.bed"),
        mask_gz = "results/bed/{prefix}.SNPS.ml.callable.flipped.{chr}.bed.gz",
        mask_idx = "results/bed/{prefix}.SNPS.ml.callable.flipped.{chr}.bed.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        last_pos=$(tail -1 {input.raw_bed} | cut -f3)
        python3 {input.script} -i {input.bed} -l $last_pos -o {output.mask}

        bgzip -c {output.mask} > {output.mask_gz}
        tabix -p bed {output.mask_gz}
        """

rule prepare_smcpp_ml:
    """
        Prepare .smc file per chr with the callability.bed as mask
    """
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz",
        vcf_idx = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi",
        pop_path = "results/{prefix}.pop",
        fai = "results/stats/snps/ml/{prefix}.SNPS.ml.{chr}.fai",
        mask = "results/bed/{prefix}.SNPS.ml.callable.flipped.{chr}.bed.gz",
        mask_idx = "results/bed/{prefix}.SNPS.ml.callable.flipped.{chr}.bed.gz.tbi",
        smcpp = config["smcpp"]
    output:
        smcpp_input = "results/ne_inference/smcpp/snps/ml/{prefix}.SNPS.ml.{chr}.smc.gz"
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


rule flip_bed_strict:
    """
        Reverse bed file to show positions that will be excluded
    """
    input:
        bed = "results/bed/snps/strict/{prefix}.SNPS.strict.{chr}.callable.bed",
        raw_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/reverse_bed.py"
    output:
        mask = temp("results/bed/{prefix}.SNPS.strict.callable.flipped.{chr}.bed"),
        mask_gz = "results/bed/{prefix}.SNPS.strict.callable.flipped.{chr}.bed.gz",
        mask_idx = "results/bed/{prefix}.SNPS.strict.callable.flipped.{chr}.bed.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """
        last_pos=$(tail -1 {input.raw_bed} | cut -f3)
        python3 {input.script} -i {input.bed} -l $last_pos -o {output.mask}

        bgzip -c {output.mask} > {output.mask_gz}
        tabix -p bed {output.mask_gz}
        """

rule prepare_smcpp_strict:
    """
        Prepare .smc file per chr with the callability.bed as mask
    """
    input:
        vcf = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz",
        vcf_idx = "results/vcf/snps/{prefix}.SNPS.{chr}.vcf.gz.tbi",
        pop_path = "results/{prefix}.pop",
        fai = "results/stats/snps/strict/{prefix}.SNPS.strict.{chr}.fai",
        mask = "results/bed/{prefix}.SNPS.strict.callable.flipped.{chr}.bed.gz",
        mask_idx = "results/bed/{prefix}.SNPS.strict.callable.flipped.{chr}.bed.gz.tbi",
        smcpp = config["smcpp"]
    output:
        smcpp_input = "results/ne_inference/smcpp/snps/strict/{prefix}.SNPS.strict.{chr}.smc.gz"
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

    #################################
    ### For filter on Callability ###
    #################################

rule flip_bed_na:
    """
        Reverse bed file to show position that will be excluded
    """
    input:
        bed = "results/bed/snps_na/{prefix}.SNPS.NA.{chr}.callable.bed",
        raw_bed = "results/bed/raw/{prefix}.raw.{chr}.callable.bed",
        script = "/groups/plantlp/vcf_processing/scripts/snakemake_prepareNe/workflow/scripts/reverse_bed.py"
    output:
        mask = temp("tmp/{prefix}.SNPS.NA.{chr}.callable.flipped.bed"),
        mask_gz = "results/bed/{prefix}.SNPS.NA.{chr}.callable.flipped.bed.gz",
        mask_idx = "results/bed/{prefix}.SNPS.NA.{chr}.callable.flipped.bed.gz.tbi"
    conda:
        "../envs/vcf_processing.yml"
    shell:
        """ 
        last_pos=$(tail -1 {input.raw_bed} | cut -f3)
        python3 {input.script} -i {input.bed} -l $last_pos -o {output.mask}

        bgzip -c {output.mask} > {output.mask_gz}
        tabix -p bed {output.mask_gz}
        """

rule prepare_smcpp_na:
    """
        Prepare .smc file per chr with the callability.bed as mask
    """
    input:
        vcf = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz",
        vcf_idx = "results/vcf/snps_na/{prefix}.SNPS.NA.{chr}.vcf.gz.tbi",
        pop_path = "results/{prefix}.pop",
        fai = "results/stats/snps_na/{prefix}.SNPS.NA.{chr}.fai",
        mask = "results/bed/{prefix}.SNPS.NA.{chr}.callable.flipped.bed.gz",
        mask_idx = "results/bed/{prefix}.SNPS.NA.{chr}.callable.flipped.bed.gz.tbi",
        smcpp = config["smcpp"]
    output:
        smcpp_input = "results/ne_inference/smcpp/snps_na/{prefix}.SNPS.NA.{chr}.smc.gz"
    conda:
        "../envs/vcf_processing.yml"
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