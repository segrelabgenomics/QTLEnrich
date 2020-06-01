#!/usr/bin/env python3

"""
create_null_variants.py is a python script to compute a set of null variants from all variants and significant variant

This script requires pandas.

# Author: Andrew Hamel
# Written in the lab of Ayellet Segre, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
# Date: 22 May 2020

##### Inputs #####

1. all_variants_file: File containing list of all variants.
2. significant_variants_file: File containing list of significant variants. 
3. output_file: name of outputted null variant file, e.g. null_variants.txt

Note: all_variants_file and significant_variants_file consist of a table with one column
with header "variant" and one variant id per line

"""


import sys
import os
import os.path
import argparse
import pandas as pd
import numpy as np

def drop_duplicates(df):
    """
    drops duplicates in variant column
    """
    df = df.drop_duplicates(["variant"])
    return(df)

def generate_null_variants(all_variants_file,significant_variants_file):
    """
    set difference between all variants and significant variants
    """
    all_variants_df = pd.read_csv(all_variants_file)
    all_variants_df = drop_duplicates(all_variants_df)   

    significant_variants_df = pd.read_csv(significant_variants_file)
    significant_variants_df = drop_duplicates(significant_variants_df) 

    null_variants = all_variants_df.loc[~all_variants_df["variant"].isin(significant_variants_df["variant"])].copy()
    return(null_variants)

def parse_args():
    """
    arguments to pass
    """
    parser = argparse.ArgumentParser(description="Determine null variants.")

    parser.add_argument("--all_variants_file",type=str,help="File containing list of all variants.",required=True)
    parser.add_argument("--significant_variants_file",type=str,help="File containing list of significant variants.",required=True)
    parser.add_argument("--output_file",type=str,help="Name of outputted null vairant fie, e.g. null_variants.txt",required=True)

    return(parser.parse_args())

if __name__ == "__main__":

    #retrieve arguments
    args = parse_args()

    #compute null variants in table format
    null_df = generate_null_variants(args.all_variants_file,args.significant_variants_file)

    #write to file
    null_df.to_csv(args.output_file,header=True,index=None,sep='\t')

