#!/usr/bin/env python3
"""
This python script parses a sample GWAS from the UK Biobank computed by the Neale lab (http://www.nealelab.is/uk-biobank) and runs liftOver from UCSC so that it is compatible with GTEx release v8 genome build GRCh38.

Author: Andrew Hamel
Date: 10 February 2020

"""

# ADD HERE command to load liftOver module

import os
import os.path
import re
import argparse
import tempfile
import subprocess
import pandas as pd

def parse_gwas(gwas_file):
    """
    parses GWAS:
     1. replaces colon (:) with underscore (_)
     2. splits variant header to extract chr, pos, non_effect, and effect_allele
     3. changes chr and pos to appropriate datatypes
     4. p-value heaader to gwas_p_value
     5. id header is used to merge back to older gwas
    """
    gwas = pd.read_csv(gwas_file,sep="\t")

    gwas["variant"] = gwas["variant"].str.replace(":","_")
    gwas["chr"],gwas["pos"],gwas["non_effect_allele"],gwas["effect_allele"] = gwas["variant"].str.split("_",3).str
    gwas["chr"] = "chr"+gwas["chr"]
    gwas["pos"] = gwas["pos"].astype(int)
    gwas["id"] = gwas["chr"]+"_"+gwas["pos"].astype(str)
    gwas = gwas.rename(columns={"pval":"gwas_p_value"})
    
    return(gwas)

def prepare_gwas_liftover(gwas):
    """
    prepares gwas so that it can be run on liftOver:
     1. extracts chr,pos headers
     2. start header is pos-1
    """
    gwas_liftover = gwas[["chr","pos","id"]].copy()
    gwas_liftover["start"] = gwas_liftover["pos"]-1
    gwas_liftover = gwas_liftover[["chr","start","pos","id"]]
    return(gwas_liftover)    

def write_temp_file(gwas_liftover):
    """
    write gwas_liftover df to temp file
    """
    temp = tempfile.NamedTemporaryFile()
    gwas_liftover.to_csv(temp.name,sep="\t",header=None,index=None)
    return(temp)

def run_liftover(temp,chain_file):
    """
    run liftover
    """
    liftover_command = ["liftOver",temp.name,chain_file,"mapped","unmapped"]
    subprocess.call(liftover_command)

def reparse_gwas(gwas):
    """
    obtain new chromosome/position pairs
    """
    #add headers to mapped file
    new_gwas = pd.read_csv("mapped",header=None,sep="\t")
    new_gwas.columns = ["new_chr","start","new_pos","id"]

    #merge with old gwas, remove old headers, rename new headers
    new_gwas = pd.merge(new_gwas,gwas,on=["id"])
    new_gwas = new_gwas.drop(["chr","start","pos","id"],axis=1)
    new_gwas = new_gwas.rename(columns={"new_chr":"chr","new_pos":"pos"})

    return(new_gwas)

def parse_args():
    """
    arguments to run script
    """
    parser = argparse.ArgumentParser(description="liftOver GWAS")

    parser.add_argument("--gwas_file",type=str,help="gwas trait file of interest.",required=True)
    parser.add_argument("--chain_file",type=str,help="chain file used in liftOver.",required=True)
    parser.add_argument("--output_file",type=str,help="name of output file.",required=True)

    return parser.parse_args()

if __name__ == "__main__":

    args = parse_args()

    print("reformatting gwas file")
    gwas = parse_gwas(args.gwas_file)

    print("formatting gwas for liftOver")
    gwas_liftover = prepare_gwas_liftover(gwas)

    print("creating temporary file...")
    temp = write_temp_file(gwas_liftover)

    print("liftover...")
    run_liftover(temp,args.chain_file)

    print("reparsing gwas and writing to new file")
    new_gwas = reparse_gwas(gwas)
    new_gwas.to_csv(args.output_file,sep="\t",index=None,compression="gzip")


