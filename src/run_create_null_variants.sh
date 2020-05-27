#!/bin/bash

# Shell script to that runs create_null_variants.py
# create_null_variants.py was written in python3 and requires pandas

# Author: Andrew Hamel
# Written in the lab of Ayellet Segre, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
# Date: 22 May 2020

##### Inputs #####
# all_variants_file: File containing list of all variants.
# significant_variants_file: File containing list of significant variants.
# output_file: name of outputted null variant file, e.g. null_variants.txt

# Note: all_variants_file and significant_variants_file consist of a table with one column
# with header "variant" and one variant id per line
#

all_variants_file=""
significant_variants_file=""
output_file=""

#purge and load modules as needed
#module purge
#module load <>

./create_null_variants.py --all_variants_file $all_variants_file \
                          --significant_variants_file $significant_variants_file \
                          --output_file $output_file

