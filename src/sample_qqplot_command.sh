#!/bin/bash
#
# Python script for generating Q-Q plots from QTLEnrich output files
#
# Author: Andrew Hamel
# Date: 24 January 2020
#
# Written in the Segre lab, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
#
# Required libraries: Python3, scipy, numpy, pandas, and matplotlib
#
#
# Example script for generating a Q-Q plot for Height GWAS p-values of a subest of Whole Blood eQTLs from GTEx v8.
# This code plots 2 Q-Q plots, one with GWAS p-values of all variants tested in the GWAS as background, and another
# with GWAS p-values of 100 permuted confounder-matched null variant sets as background.
#
# Details of all options can be found in src/READM.md
#

gwas_pvalues="../data/FULL_GWAS_PVAL_FILENAME" # GWAS p-value of the full GWAS, required column: 'gwas_p_value'
tissue_pvalues="../examples/Output_QTLEnrich_example_Mar02_2020/output_files/QTLEnrich_significant_gwas_p_values_best_eqtl_50_irnt.gwas.imputed_v3.both_sexes.GRCh38_Whole_Blood_Mar02_2020.txt" # Required columns: 'variant', 'gwas_p_value'
matrix="../examples/Output_QTLEnrich_example_Mar02_2020/output_files/QTLEnrich_matched_null_gwas_p_values_best_eqtl_50_irnt.gwas.imputed_v3.both_sexes.GRCh38_Whole_Blood_Mar02_2020.txt" # GWAS p-values of 1000 permuted confounder-matched null variant sets
num_permutations=100 # number permutatoins to include from matrix file for confounder-matched null background plotting
tissue="Whole_Blood"
qtl="best_eQTL"
title="Height GWAS"

./qqplot_qtlenrich.py --gwas_file $gwas_pvalues \
                      --qtl_file $tissue_pvalues \
		      --pvalue_null_matrix $matrix \
                      --gwas_background \
                      --matched_null_background \
		      --num_permutations \
                      --confidence_interval \
                      --show_title \
                      --qqplot_title $title \
                      --qtl_type $qtl \
                      --tissue_label $tissue

