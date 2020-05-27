#!/bin/bash

# Script used to generate confounders table for GTEx v8
# Author: Andrew Hamel
# Date: 30 Oct. 2019
#
# Written in the Segre lab at Massachusetts Eye and Ear, Department of Ophthalmology
# and Ocular Genomics Institute
# Harvard Medical School, Boston, MA
#

QTL_Directory=""
File_Extension=".allpairs.txt.gz"
Variants_List="../../data/GTEx_V8_unique_variants.txt"
LD_Proxy="../../data/GTEx_V8_ld_proxy_variants.txt"
GENCODE='../../data/GENCODE_26_df.txt'
Output_File="../../data/GTEx_V8_confounders_table.txt.gz"

#purge and load modules as necessary
#module purge
#module <>

./Generate_Confounders_Table.py --QTL_Directory $Q \
  --File_Extension $File_Extension \
  --Variants_List $Variants_List \
  --LD_Proxy $LD_Proxy \
  --Output_File $Output_File \
  --GENCODE_File $GENCODE  


