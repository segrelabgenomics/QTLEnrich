#!/bin/bash
#
# This is a sample command for runninng QTLEnrich
#
# Author: Andrew Hamel
# Date: 24 January 2020
#
# Written in the Segre lab, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
#
# Required libraries listed in src/README.md
#
# Details of all options can be found in src/README.md
#

file='../data/GWASFILE'
QTL_Directory='../examples/'
File_Extension='.v8.egenes.txt.gz' # suffix of eGenes file in GTEx v8
QTL_Type='best_eqtl'
Confounders_Table='../data/GTEx_v8_Confounders_Table_49Tiss_Global.txt.gz'
Null_Table='../data/GTEx_v8_Null_Table_49Tiss_Global.txt.gz'
gencode_file='../data/GENCODE_26_df.txt'
exp_label='example'

#purge and load modules as necessary

./QTLEnrichV2.py --gwas_file $file \
                 --qtl_directory $QTL_Directory \
                 --file_name $File_Extension \
                 --qtl_type $QTL_Type \
                 --confounders_table $Confounders_Table \
                 --null_table $Null_Table \
                 --gencode_file $gencode_file \
                 --exp_label $exp_label \
                 --GeneEnrich_input
