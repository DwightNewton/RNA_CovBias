# RNA_CovBias
Pipeline for pulling Picard QC metrics of RNA-sequencing data, and generating inferred measures of RNA quality via transcript coverage bias.

Estimating RNA quality in pooled-cell samples is difficult given low starting material. RNAseq data itself can provide an inferred measure of RNA quality by assessing 5' or 3' bias in coverage. As RNA degradation primarily occurs 3'-> 5', a greater 5' bias would indicate greater degradation. This pipeline uses the Picard toolkit (Broad Institute) to generate normalized coverage distributions across highly expressed transcripts of all RNAseq libraries in an experiment. These distributions are then compared (Kolmogorov–Smirnov test) across contrasts of interest (primarily experimental condition or cell-type) to determine if there are any confounding differences in quality.

Required packages:
* GATK/Picard toolkit (CAMH SCC has pre-installed modules)
* SAMtools
* Base R, ggplot2, reshape2, multtest, gridExtra

Scripts below were written for use with Sibille Lab SCT-RNAseq PITT Tetrad dataset.

Steps:

1. Convert existing GTF file (used for alignment) to *refFlat* format and re-order for compatibility.
   * Run `<gtftoGenePred>` function by UCSC to generate raw *refFlat* file.
   * Use `<refFlat convert.R>`.
   * Adapt input and output filenames as appropriate.

2. Sort .bam files using SAMtools if not already sorted.
   * `<./samtools_sort_script.sh > samsort.cmdlist>` to create joblist.
   * Run in parallel `<samsort_parallelization.sh>` - **6+ hours is a safe wall-time generally**

3. Run Picard in parallel for all libraries using SLURM.
   * Create job list using `<picard_cmdgen.sh>` and run with `<picard_parallelization.sh>`.
   * Change job # and working directories as appropriate.
   * Direct output to new directory.

5. Picard output is non-tabular - some wrangling required.
   * Run `<picard output merge.R>`.

6. Run analysis script - `<picard results and plots.R>`
   * Produces plots comparing coverage distributions across contrasts of interest.
   * Excel files generated contain *post-hoc* Bonferroni-corrected t-tests and effect sizes. Only distributions which are significantly different by KS-test should be extracted.
