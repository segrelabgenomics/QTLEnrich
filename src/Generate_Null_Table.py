#!/usr/bin/env python3

"""
Author: Andrew Hamel
Date: 28 March 2019

Python3 script to generate null variants in the correct format to be used in QTLEnrich.

Requires at least pandas release 0.24

Arguments:
    1. Directory
    2. file_extension (files in the directory to be parsed e.g. .allpairs.txt.gz)
    3. Null_Variant_List (global null)
    4. sQTL significant variants list (here the null variants intersected with significant sQTL variants are in global null)
    5. Output_File (name of file created with path. If path not included, defaults to current working directory.)

Returns:
    1. Null table

Null Table format:
    Headers
        A. variant (variant_id in the form chrN_pos_ref_alt_b38)
        B. Tissue_tissueName (For each tissue, 1 indicates null variant, 0 indicates other, example column name 'Tissue_Whole_Blood')

"""

import os
import os.path
import re
import argparse
import pandas as pd
import numpy as np

def format_variants(significant_variants_file):
    """
    Set sQTL list or null variants list to numpy array
    """
    x = pd.read_csv(significant_variants_file,sep='\t')
    return (np.array(x['variant']))

def parse_QTLs(QTL_File,QTL_Directory,Null_Variants_List_Array,Tissue,qtl_significant_variants=""):
    """
    Parses .allpairs file 

    Finds intersection with allpairs and Null_Variants_List
    Removes significant variants

    Revert back to pandas dataframe
    Each variant will be labeled as 1
    """
    print('reading file')
    QTL = pd.read_csv(QTL_Directory+QTL_File,sep='\t',usecols=['variant_id'])
    QTL = QTL.rename(columns={'variant_id':'variant'})

    #remove duplicates
    print('removing duplicates')
    QTL = QTL.drop_duplicates()

    #find intersection with global null
    print('finding intersection')
    Null = QTL[QTL['variant'].isin(Null_Variants_List_Array)].copy()

    #remove variants found in sQTL list
    if qtl_significant_variants:
        print('removing variants found in sQTL list')
        Null = Null[~Null['variant'].isin(qtl_significant_variants)].copy()

    #back to dataframe and set those variants to 1
    Null['Tissue_'+Tissue] = 1

    return Null_df

def parse_args():
    """
    Arguments to generate null table
    """
    parser = argparse.ArgumentParser(description="Generate null table for QTLEnrich")
    parser.add_argument('-D','--QTL_Directory',type=str,help='Directory containing QTL files to parse',required=True)
    parser.add_argument('-F','--File_Extension',type=str,help='File to extension to parse files in QTL directory',required=True)
    parser.add_argument('--Null_Variants_List',type=str,help='List of null variants',required=True)
    parser.add_argument('--QTL_Significant_Variants_List',type=str,default="None",help='significant QTL variants')
    parser.add_argument('-O','--Output_File',type=str,help='Name of generated file',required=True)

    return parser.parse_args()

if __name__ == '__main__':

    args = parse_args()

    #extract relevant files
    print('selecting appropriate files.')
    QTL_Files = [x for x in os.listdir(args.QTL_Directory) if x.endswith(args.File_Extension)]

    #read in null variants list
    print('Reading global null.')
    Null_Table_df = pd.read_csv(args.Null_Variants_List)

    #Create np array for Null_Variants_List
    Null_Variants_List_Array = np.array(Null_Table_df['variant'])

    #Extract significant sQTL variants
    if args.QTL_Significant_Variants_List != "None":
        significant_variants = format_variants(args.QTL_Significant_Variants_List)
    else:
        significant_variants = ""

    for i in QTL_Files:

        #extract tissue name
        tissue_name=re.sub(args.File_Extension,'',i)

        #Determine null variants for each tissue
        Null_df = parse_QTLs(i,args.QTL_Directory,Null_Variants_List_Array,tissue_name,qtl_significant_variants=significant_variants)

        #merge with Null_Table_df
        Null_Table_df = pd.merge(Null_Table_df,Null_df,on=['variant'],how='outer')

        #set na to 0
        Null_Table_df = Null_Table_df.fillna(0)

        #set floats to integers
        Null_Table_df['Tissue_'+tissue_name] = Null_Table_df['Tissue_'+tissue_name].astype('Int64')

        print(Null_Table_df.head())
 
    Null_Table_df.to_csv(args.Output_File,header=True,index=None,sep='\t',compression='gzip')


