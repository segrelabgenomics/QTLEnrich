#!/usr/bin/env python3

"""
Author: Andrew Hamel
Date: 27 March 2019

Python3 script to parse .allpairs.txt.gz files in GTEx to create confounders table to run in QTLEnrich.

Requires at least pandas release 0.24

Arguments:
    1. Directory
    2. file_extension (files in the directory to be parsed e.g. .allpairs.txt.gz)
    3. Variant_List (unique list of all variants)
    4. LD_Proxy (file containing variants and the number of LD proxy variants)
    5. Output_File (name of file created with path. If path not included, defaults to current working directory.)

Returns:
    1. Confounders table

Confounder Table format:
    Headers
        A. variant (variant_id in the form chrN_pos_ref_alt_b38)
        B. LD_proxy (number of LD proxy variants for a given variant)
        C. MAF_Tissue (For each tissue-variant pair, minor allele frequency is indicated)
        D. TSS_Tissue (For each tissue-variant pair, nearest TSS distance is chosen)

"""

import os
import os.path
import argparse
import re
import pandas as pd
import numpy as np

def Format_GENCODE_File(GENCODE):
    """
    Prepare GENCODE File to merge with retina eQTL file
    """
    GENCODE_df = pd.read_csv(GENCODE,sep='\t')

    #select genes
    GENCODE_df_gene = GENCODE_df.loc[GENCODE_df['feature'] == 'gene'].copy()

    #extract relevant columns
    GENCODE_df_gene = GENCODE_df_gene[['gene_id','start','end','strand']]
    GENCODE_df_gene = GENCODE_df_gene.rename(columns = {'start':'gene_start','end':'gene_end'})

    return GENCODE_df_gene

def Merge_GENCODE_QTL(variant_df,gencode):
    """
    Merge GENCODE file and retina file on ensembl gene_id, without decimal suffix
    """
    print('merging')
    Merged_File = pd.merge(variant_df,gencode,on=['gene_id'])
    return Merged_File

def Compute_TSS_Distance(variant_df,gencode):
    """
    compute tss distance

    if strand is positive:
        tss distance = start of variant -- start of gene
    if strand is negative:
        tss distance = (start of variant -- end of gene)*(-1)
    """
    GENCODE_QTL = Merge_GENCODE_QTL(variant_df,gencode)
    print('computing TSS_Distance')

    #compute TSS Distance for positive and negative strands separately
    #append back together
    GENCODE_QTL_pos_strand = GENCODE_QTL.loc[GENCODE_QTL['strand'] == '+'].copy()
    GENCODE_QTL_neg_strand = GENCODE_QTL.loc[GENCODE_QTL['strand'] == '-'].copy()

    GENCODE_QTL_pos_strand['TSS_Distance'] = GENCODE_QTL_pos_strand['pos'] - GENCODE_QTL_pos_strand['gene_start']
    GENCODE_QTL_neg_strand['TSS_Distance'] = (GENCODE_QTL_neg_strand['pos'] - GENCODE_QTL_neg_strand['gene_end'])*(-1)

    GENCODE_QTL_Full = GENCODE_QTL_pos_strand.append(GENCODE_QTL_neg_strand)
    GENCODE_QTL_Full = GENCODE_QTL_Full.drop(['gene_start','gene_end','strand'],axis=1)

    return GENCODE_QTL_Full

def Parse_QTL_File(QTL_File,QTL_Directory,Tissue,gencode):
    """
    Parses QTLs

    1. Removes rows with MAF=0
    2. If strand < 0: multiply TSS*-1
    3. Sorts by absolute value of TSS distance
    3. Removes duplicates, selecting rows with smallest TSS distance
    4. Renames columns in form: variant,MAF_Tissue,TSS_Tissue    
    """
    QTL = pd.read_csv(QTL_Directory+QTL_File,sep='\t',usecols=['variant_id','gene_id','maf'])

    #remove rows with maf=0
    QTL = QTL.loc[QTL['maf'] != 0.0]

    #extract ensemblID
    QTL['pos'] = QTL['variant_id'].str.split('_',2).str[1].astype(int)

    print(QTL.head())

    #computing tss distance
    QTL_computed_tss_distance = Compute_TSS_Distance(QTL,gencode)

    del QTL

    #create column of absolute value of tss_distance
    #sort by abs value, ascending
    QTL_computed_tss_distance['Abs_tss_distance'] = QTL_computed_tss_distance['TSS_Distance'].abs()
    QTL_computed_tss_distance = QTL_computed_tss_distance.sort_values(by=['Abs_tss_distance'],ascending=True)

    #drop duplicates
    #keep smallest tss_distance
    QTL_computed_tss_distance = QTL_computed_tss_distance.drop_duplicates(subset=['variant_id'],keep='first')
    print('length after removing duplicates:',len(QTL_computed_tss_distance))
    #rename columns
    #drop Abs_Tss column
    print(QTL_computed_tss_distance.head())
    QTL_computed_tss_distance = QTL_computed_tss_distance.rename(columns={'variant_id':'variant','TSS_Distance':'TSS_'+Tissue,'maf':'MAF_'+Tissue})
    QTL_computed_tss_distance = QTL_computed_tss_distance.drop(['Abs_tss_distance','gene_id','pos'],axis=1)
    QTL_computed_tss_distance['MAF_'+Tissue] = QTL_computed_tss_distance['MAF_'+Tissue].round(decimals=5)

    return QTL_computed_tss_distance

def parse_args():
    """
    Arguments passed to generate confounders table
    """
    parser = argparse.ArgumentParser(description="Generate confounders table to be used on QTLEnrich")
    parser.add_argument('-D','--QTL_Directory',type=str,help='Directory containing QTL files to parse',required=True)
    parser.add_argument('-F','--File_Extension',type=str,help='File to extension to parse files in QTL directory',required=True)
    parser.add_argument('-V','--Variants_List',type=str,help='Unique list of all variants',required=True)
    parser.add_argument('-L','--LD_Proxy',type=str,help='File containing number of LD proxy variants',required=True)
    parser.add_argument('-G','--GENCODE_File',type=str,help='GENCODE file to extract strand orientation',required=True)
    parser.add_argument('-O','--Output_File',type=str,help='Name of generated file',required=True)

    return parser.parse_args()

if __name__ == '__main__':

    args = parse_args()

    print('parsing gencode')
    gencode = Format_GENCODE_File(args.GENCODE_File)

    #select files to parse
    QTL_Files = [x for x in os.listdir(args.QTL_Directory) if x.endswith(args.File_Extension)]

    #read unique variants file
    #denoted as confounders table
    Confounders_Table_df = pd.read_csv(args.Variants_List,sep='\t') 
    Confounders_Table_df = Confounders_Table_df.drop_duplicate()

    #Merge with LD Proxy file
    LD_Proxy_df = pd.read_csv(args.LD_Proxy,sep='\t')
    Confounders_Table_df = pd.merge(Confounders_Table_df,LD_Proxy_df,on=['variant'],how='outer')

    print('length of confounders table before merging:',len(Confounders_Table_df))

    #parse files
    for Count,QTL_File in enumerate(QTL_Files,1):
        print(Count,QTL_File)

        #Extract tissue name
        tissue=re.sub(args.File_Extension,'',QTL_File)

        #parse file
        #merge with preexisting table
        df = Parse_QTL_File(QTL_File,args.QTL_Directory,tissue,gencode)
        Confounders_Table_df = pd.merge(Confounders_Table_df,df,on=['variant'],how='outer')
        print('length after merging:',len(Confounders_Table_df))

    #write to file
    Confounders_Table_df.to_csv(args.Output_File,sep='\t',header=True,index=None,compression='gzip')
        
        
