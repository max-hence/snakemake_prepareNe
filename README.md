#### PrepareNE ####

Pipeline to prepare inputs for Ne inference. From raw vcf (snpArcher outputs) to stairway plot and smc++ inputs

Executes the following functions:

    - Filters variants to keep only bi-allelic SNPS
    - Rescales chr length to take into account removing of other variants
    - Corrects genotypes based on callability by sample
    - Subsamples based on amount of missing data
    - Rescales chr length to take into account loss of regions which have too many missing data
    - Generates SFS on every chr with different treatments
    - Prepare StairwayPlot and SMC++ inputs