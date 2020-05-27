#!/bin/bash

# Shell script to concatenate outputs of QTLEnrich into one file

# Author: Andrew Hamel
# Written in the lab of Ayellet Segre, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
# Date: 12 February 2020

##### Inputs #####
# directory: directory where QTLEnrich output files were generated
# concatenated_filename: designated name of concatenated QTLEnrich set, e.g. /path/to/file/concatenated_set.txt

# Note: depending on the presence of other files in directory,
# user may need to change *.txt to concatenate only relevant files
# Currently assumes directory consists of relevant files

directory=$1
concatenated_filename=$2

#change to relevant directory
#concatenate files
cd $directory
awk 'FNR==1 && NR!=1{next;}{print}' *.txt > $concatenated_filename


