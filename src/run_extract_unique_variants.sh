#!/bin/bash

#extract_unique_variants.py is a python script to extract a list of unique variants from multiple tissues
#
#This script requires pandas.
#
# Author: Andrew Hamel
# Written in the lab of Ayellet Segre, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
# Date: 22 May 2020
#
##### Inputs #####
#
# directory: directory containing qtl variants for all relevant tissues
# file_name: the suffix of all QTL files that will be read in the directory
# output_file: Name of outputted unique variant file, e.g. unique_variants.txt
#
#Note: files within the directory should contain the header 'variant_id' and are to be tab-delimited, such as in GTEx

directory=""
file_name=""
output_file=""

#purge and load modules as needed
#module purge
#module load <>

./extract_unique_variants.py --directory $directory \
                             --file_name $file_name \
                             --output_file $output_file


