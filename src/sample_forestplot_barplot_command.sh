#!/bin/bash
#
# Python script for generating forest plots and barplots of adjusted fold-enrichment, estimated number of traits associations amongt QTLs, or empirical Pi1
# from QTLEnrich output files
#
# Author: Andrew Hamel
# Date: 24 January 2020
#
# Written in the Segre lab, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
#
# Required libraries: ggplo2, cowplot, readr, and dplry in R
#
# Details of all options can be found in src/README.md

qtlenrich_file="../examples/Output_QTLEnrich_example_49_GTEx_tissues_Feb12_2020/results/QTLEnrich_best_eqtl_50_irnt.gwas.imputed_v3.both_sexes.GRCh38.tsv.gz_example_full_tissues_Feb12_2020.txt"
qtlenrich_file2="../examples/Output_QTLEnrich_example_49_GTEx_tissues_Feb12_2020/results/QTLEnrich_best_sqtl_50_irnt.gwas.imputed_v3.both_sexes.GRCh38.tsv.gz_example_full_tissues_Feb12_2020.txt"
label1="eQTL"
label2="sQTL"
gtex_tissue_names_file="../data/gtex_tissue_colors_v8.txt"
exp_label="example"
significant_cutoff=0.00051 # Bonferroni correction for testing 49 tissues and two QTL types at alpha=0.05

# load R package 

./forest_boxplot_module.R --qtlenrich_file1 $qtlenrich_file \
                              --qtlenrich_file2 $qtlenrich_file2 \
                              --label_file1 $label1 \
                              --label_file2 $label2 \
                              --significant_option \
                              --significant_tissues_cutoff $significant_cutoff \
                              --identify_number_significant_tissues 49 \
                              --gtex_tissue_names $gtex_tissue_names_file \
                              --exp_label $exp_label \
                              --forest_plot \
                              --bar_plot \
                              --paired_plot \
                              --side_by_side

