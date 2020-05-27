#!/bin/bash

# This is a sample command for runninng QTLEnrich on a UK Biobank Height GWAS and whole blood eQTLs from GTEx v8
# Input GWAS and eQTL files should be saved in examples/ directory

# Author: Andrew Hamel
# Date: 24 January 2020
#
# Written in the Segre lab, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
#
# Required libraries listed in src/README.md
#
# Details of all options can be found in src/README.md
#

GWAS_file="../examples/50_irnt.gwas.imputed_v3.both_sexes.GRCh38.tsv.gz" # genome build GRCh38
QTL_Directory='../examples/'
File_Extension='.v8.egenes.txt.gz'
QTL_Type='best_eqtl'
Confounders_Table='../data/GTEx_v8_Confounders_Table_49Tiss_Global.txt.gz' # Distance to TSS, MAF, number of LD proxy variants per variant
Null_Table='../data/GTEx_V8_Null_Table_49Tiss_Global.txt.gz'  # table that specifies non-significant (null) eQTL or sQTL variants in each tissue (based on genes expressed in given tissue) or globally
gencode_file='../data/GENCODE_26_df.txt'
exp_label='test_run'

#purge and load modules as necessary
#module purge
#module load <>

./QTLEnrichV2.py --gwas_file $GWAS_file \
                 --qtl_directory $QTL_Directory \
                 --file_name $File_Extension \
                 --qtl_type $QTL_Type \
                 --confounders_table $Confounders_Table \
                 --null_table $Null_Table \
                 --gencode_file $gencode_file \
                 --exp_label $exp_label \
                 --keep_null_variants
