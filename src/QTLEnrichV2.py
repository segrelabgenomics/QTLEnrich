#!/usr/bin/env python3

"""
QTLEnrichV2 was written under the direction and supervision of Ayellet Segre
at Massachusetts Eye and Ear, Department of Ophthalmology and Ocular Genomics Institute

Author: Andrew Hamel

Code written in Python version 3.6

QTLEnrichV2 assesses enrichment of gwas variant associations amongst qtls in
a given tissue correcting for potential confounding factors 
(MAF, distance to TSS, local LD).
"""

import sys
import os
import os.path
import re
import argparse
import logging
import datetime
import pandas as pd 
import numpy as np
import os_sys_qtlenrich
import qtl_gwas_parsing 
import rand_null_samp_permute_matching_confounders
import compute_fold_enrichments
import output_table
import GeneEnrich_inputs

if __name__ == "__main__":

    args = os_sys_qtlenrich.parse_args()

    #check specific arguments
    os_sys_qtlenrich.raise_error_gencode_file(args)

    #retrieve dates
    date,month_date = os_sys_qtlenrich.retrieve_current_date()

    #parse exp label
    exp_label = os_sys_qtlenrich.parse_exp_label(args.exp_label)

    #Extract gwas Trait name
    gwas_trait = qtl_gwas_parsing.extract_gwas_trait_name(args.gwas_file,args.trait_name)

    #Create output directory
    output_directory = os_sys_qtlenrich.create_output_directory(args,month_date,exp_label)

    #create pvalues directory
    output_files_directory = os_sys_qtlenrich.create_output_files_directory(output_directory)

    #extract arguments
    arguments,parameters = os_sys_qtlenrich.extract_qtlenrich_arguments(args)

    #create and write to log file
    log_file = os_sys_qtlenrich.write_initial_logfile(output_directory,date,gwas_trait,args.qtl_type,arguments,parameters,exp_label)

    #Select appropriate files
    qtl_files = qtl_gwas_parsing.extract_qtl_files(args.qtl_directory,args.file_name)

    if args.subset_tissues:
        qtl_files = os_sys_qtlenrich.subset_qtl_files(qtl_files,args.file_name,args.tissue_names)

    #interaction qtl dict
    interaction_qtl_dict = qtl_gwas_parsing.interaction_qtl_dictionary(args.qtl_type,args.file_name,qtl_files)

    #parse gwas
    print("preparing gwas")
    with open(log_file,"a+") as f:
        f.write("Preparing gwas\n")
    gwas = qtl_gwas_parsing.prepare_gwas_file(args.gwas_file,args.genome_build)

    print("Preparing confounders and null table...\nThis may take a few minutes\n")
    with open(log_file,"a+") as f:
        f.write("Preparing confounders and null table...\nThis may take a few minutes\n")

    #parse confounders and null tables
    confounders_table = qtl_gwas_parsing.read_confounders_table(args.confounders_table)
    null_table = qtl_gwas_parsing.prepare_null_table(args.null_table,gwas)

    print("parsing qtls\n")
    with open(log_file,"a+") as f:
        f.write("parsing qtls\n")

    gencode = qtl_gwas_parsing.declare_gencode(args)

    significant_qtl_dict,significant_qtl_gwas_dict,geneenrich_input_gwas_dict = qtl_gwas_parsing.parse_qtls(output_files_directory,month_date,gwas_trait,qtl_files,args.qtl_directory,args.qtl_type,args.file_name,gwas,args.qtl_q_value,args.independent_ranking,args.compute_tss_distance,gencode,args.subset_genes,args.null_option,keep_gene_id_suffix)

    #Generate GeneEnrich Inputs if necessary
    if args.GeneEnrich_input:
        if args.qtl_type == "dapg_independent" or args.qtl_type == "conditional_independent":

            #GeneEnrich input files for independent eQTLs requires additional files
            os_sys_qtlenrich.check_independent_GeneEnrich_parameters(args)

            GeneEnrich_inputs.create_GeneEnrich_inputs_independent_qtls(args.eGenes_directory,args.eGene_file_name,gwas,gwas_trait,output_directory,date,geneenrich_input_gwas_dict,args.qtl_type,p=args.gwas_p_value,q=args.qtl_q_value,independent_ranking=args.independent_ranking,gencode=gencode,subset_genes=args.subset_genes)
        else:
            GeneEnrich_inputs.create_GeneEnrich_inputs_best_eqtl_sqtl(gwas_trait,output_directory,date,geneenrich_input_gwas_dict,args.qtl_type,exp_label,p=args.gwas_p_value,q=args.qtl_q_value)

    del geneenrich_input_gwas_dict, gwas

    print("computing observed fold-enrichment")
    with open(log_file,"a+") as f:
        f.write("computing observed fold-enrichment\n")

    #Compute qtl Statistics
    trait_dict,original_length_dict,observed_fold,obs_dict,best_eqtl_variants_dict = compute_fold_enrichments.compute_qtl_statistics(gwas_trait,significant_qtl_dict,significant_qtl_gwas_dict,args.gwas_p_value)

    #perform random sampling of matched null variants
    enrichment_p_value_dict,adjusted_fold_enrichment_dict,upper_bound_confidence_interval_dict,lower_bound_confidence_interval_dict,pi_1_dict,true_trait_dict = rand_null_samp_permute_matching_confounders.sample_match_null_confounders(log_file,month_date,output_files_directory,args.qtl_type,gwas_trait,confounders_table,null_table,significant_qtl_gwas_dict,observed_fold,args.gwas_p_value,null_option=args.null_option,num_permutations=args.lower_bound_permutations,upper_bound_permutations=args.upper_bound_permutations,lambda_factor=args.lambda_factor,independent_ranking=args.independent_ranking,compute_tss_distance=args.compute_tss_distance,keep_null_variants=args.keep_null_variants,keep_pvalue_matrix=args.keep_pvalue_matrix,interaction_qtl_dict=interaction_qtl_dict,num_quantiles=args.num_quantiles)

    #create Output table
    output_table = output_table.create_output_table(trait_dict,original_length_dict,best_eqtl_variants_dict,observed_fold,obs_dict,enrichment_p_value_dict, adjusted_fold_enrichment_dict,upper_bound_confidence_interval_dict,lower_bound_confidence_interval_dict,pi_1_dict,true_trait_dict)

    #write output tables
    results_directory = os_sys_qtlenrich.create_results_directory(output_directory)
    output_filename = os_sys_qtlenrich.create_output_filename(results_directory,month_date,gwas_trait,args.qtl_type,args.independent_ranking,exp_label)

    with open(log_file,"a+") as f:
        f.write("Creating QTLEnrich output table. QTLEnrich is complete.\n")

    output_table.to_csv(output_filename,index=None,sep="\t")

