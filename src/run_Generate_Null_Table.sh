#!/bin/bash

# Script used to generate null table for GTEx v8
# Author: Andrew Hamel
# Date: 30 Oct. 2019
#
# Written in the Segre lab at Massachusetts Eye and Ear, Department of Ophthalmology
# and Ocular Genomics Institute
# Harvard Medical School, Boston, MA
#



QTL_Directory=""
File_Extension=".allpairs.txt.gz"
Variants_List="../../data/GTEx_V8_null_variants.txt"
significant_variants="../../data/significant_variants.txt"
Output_File="../../data/GTEx_V8_Null_Table.txt.gz"

#purge and load modules as necessary
#module purge
#module <>

./Generate_Null_Table.py --QTL_Directory $QTL_Directory \
  --File_Extension $File_Extension \
  --Null_Variants_List $Variants_List \
  --QTL_Significant_Variants_List $significant_variants \
  --Output_File $Output_File

