#!/usr/bin/env python3

"""
create_ld_table.py is a python script to compute a set of table of variants and the number of associated ld proxy variants

This script requires pandas.

# Author: Andrew Hamel
# Written in the lab of Ayellet Segre, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
# Date: 22 May 2020

##### Inputs #####

1. PLINK Output Directory: directory containing outputs from plink command used to compute ld proxy variants
2. Output file: name of outputted ld variant table, e.g. ld_proxy_table.txt

"""

import os
import os.path
import pandas as pd

def retrieve_plink_output_files(plink_output_directory):
    """
    return list of files that end in '.ld'
    """
    return([file for file in os.listdir(plink_output_directory) if file.endswith(".ld")])

def create_ld_proxy_variants_table(plink_output_directory,ld_files):
    """
    Generates a table of variant with number of LD Proxy variants
    """
    if plink_output_directory[-1] != "/":
        plink_output_directory = plink_output_directory + "/"

    ld_variant_table = pd.DataFrame()

    for file in ld_files:
        df = pd.read_csv(plink_output_directory+file,sep="\s+")

        #create value counts
        value_counts_plink = df["SNP_A"].value_counts()
        df_value_counts = pd.DataFrame(value_counts_plink).reset_index(level = 0)
        df_value_counts.columns = ["variant","num_LD_variants"]

        #append to one dataframe
        ld_variant_table = ld_variant_table.append(df_value_counts)

    return(ld_variant_table)

def parse_args():
    """
    arguments to pass
    """
    parser = argparse.ArgumentParser(description="Compute LD variant table.")

    parser.add_argument("--plink_output_directory",type=str,help="directory containing outputs from plink command used to compute ld proxy variants",required=True)
    parser.add_argument("--output_file",type=str,help="Name of outputted ld variant table, e.g. ld_proxy_table.txt",required=True)

    return(parser.parse_args())

if __name__ == "__main__":

    #retrieve arguments
    args = parse_args()

    #retrieve plink output files
    ld_files = retrieve_plink_output_files(args.plink_output_directory)

    #compute ld variant table
    ld_variant_table = create_ld_proxy_variants_table(args.plink_output_directory,ld_files)

    #write to file
    ld_variant_table.to_csv(args.output_file,header=True,index=None,sep='\t')



