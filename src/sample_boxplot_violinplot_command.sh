#!/bin/bash
#
# Example for generating box plots of the adjusted fold-enrichment for eQTLs and sQTLs across 49 GTEx tissues (v8) and Height GWAS from the UK Biobank.
#
# Author: Andrew Hamel
# Date: 24 January 2020
#
# Written in the Segre lab, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
#
# Required libraries: ggplot2, cowplot, readr, dplyr in R
#
# Details of all options can be found in src/READM.md
#

qtlenrich_file="../examples/Output_QTLEnrich_example_49_GTEx_tissues_Feb12_2020/results/QTLEnrich_best_eqtl_50_irnt.gwas.imputed_v3.both_sexes.GRCh38.tsv.gz_example_full_tissues_Feb12_2020.txt"
qtlenrich_file2="../examples/Output_QTLEnrich_example_49_GTEx_tissues_Feb12_2020/results/QTLEnrich_best_sqtl_50_irnt.gwas.imputed_v3.both_sexes.GRCh38.tsv.gz_example_full_tissues_Feb12_2020.txt"
label="eQTL"
label2="sQTL"
gtex_tissue_names_file="../data/gtex_tissue_colors_v8.txt"
remove_subset_tissues="../examples/sample_tissues_remove.txt" # subset of tissues to remove of full list of tissues tested, if any
exp_label="Height_49_tiss"
significant_cutoff=0.00051 # Bonferroni	correction for number of tissues tested and two QTL types at alpha=0.05

# load R package

module load R/3.6.1

./forest_boxplot_module.R --qtlenrich_file1 $qtlenrich_file \
                              --label_file1 $label \
                              --qtlenrich_file2 $qtlenrich_file2 \
                              --label_file2 $label2 \
                              --significant_option \
                              --significant_tissues_cutoff $significant_cutoff \
                              --gtex_tissue_names $gtex_tissue_names_file \
                              --remove_subset_tissues $remove_subset_tissues \
                              --exp_label $exp_label \
                              --box_plot \
                              --paired_plot
