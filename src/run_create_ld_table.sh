#!/bin/bash

#create_ld_table.py is a python script to compute a set of table of variants and the number of associated ld proxy variants
#
#This script requires pandas.
#
# Author: Andrew Hamel
# Written in the lab of Ayellet Segre, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
# Date: 22 May 2020
#
##### Inputs #####
#
# PLINK Output Directory: directory containing outputs from plink command used to compute ld proxy variants
# Output file: name of outputted ld variant table, e.g. ld_proxy_table.txt
#

plink_directory=""
output_file=""

#purge and load modules as needed
#module purge
#module load <>

./create_ld_table.py --plink_output_directory $plink_directory \
                          --output_file $output_file

