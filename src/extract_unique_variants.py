#!/usr/bin/env python3

"""
extract_unique_variants.py is a python script to extract a list of unique variants from multiple tissues

This script requires pandas.

# Author: Andrew Hamel
# Written in the lab of Ayellet Segre, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
# Date: 22 May 2020

##### Inputs #####

1. directory: directory containing qtl variants for all relevant tissues
2. file_name: the suffix of all QTL files that will be read in the directory
3. output_file: Name of outputted unique variant file, e.g. unique_variants.txt

Note: files within the directory should contain the header 'variant_id' and are to be tab-delimited, such as in GTEx 

"""

import sys
import os
import os.path
import argparse
import pandas as pd

def extract_tissue_files(directory,file_name):
    """
    return list of files in directory
    """
    return([file for file in os.listdir(directory) if file.endswith(file_name)]

def extract_unique_variants(directory,tissue_files):
    """
    Generate a single columnn file of unique variants for a given directory
    """
    #format directory
    if directory[-1] != "/":
        directory = directory + "/"

    unique_variant_df = pd.DataFrame()

    for tissue in tissue_files:
        tissue_df = pd.read_csv(directory+tissue,sep="\t")
        tissue_df = tissue_df[['variant_id']]
        tissue_df = tissue_df.drop_duplicates()

        unique_variant_df = unique_variant_df.append(tissue_df)
        unique_variant_df = unique_variant_df.drop_duplicates()

    #rename variant_id to variant if using for future scripts in QTLEnrich
    unique_variant_df = unique_variant_df.rename(columns={"variant_id":"variant"})
    return(unique_variant_df)

def parse_args():
    """
    arguments to pass
    """
    parser = argparse.ArgumentParser(description="Determine unique variants.")

    parser.add_argument("--directory",type=str,help="directory containing qtl variants for all relevant tissues",required=True)
    parser.add_argument("--file_name",type=str,help="the suffix of all QTL files that will be read in the directory",required=True)
    parser.add_argument("--output_file",type=str,help="Name of outputted unique variant file, e.g. unique_variants.txt",required=True)

    return(parser.parse_args())


if __name__ == "__main__":

    #retrieve arguments
    args = parse_args()

    #retrieve list of files
    tissue_files = extract_tissue_files(args.directory,args.file_name)

    #obtain unique variants from tissue files
    unique_variants = extract_unique_variants(directory,tissue_files)

    #write to file
    unique_variants.to_csv(args.output_file,header=True,sep="\t",index=None)
