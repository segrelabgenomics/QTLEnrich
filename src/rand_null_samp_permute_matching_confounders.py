
"""
this script consists of functions related to the permutation scheme of QTLEnrich
"""

import os
import os.path
import pandas as pd
import numpy as np
import sys
import time
import re
from functools import partial
from multiprocessing import Pool
import collections
from decimal import Decimal
from os_sys_qtlenrich import make_directory

def print_pvalue_matrix(log_file,output_files_directory,date,return_matrix,qtl_type,trait,tissue,independent_ranking):
    """
    returns pvalue matrix
    """
    tissue_removed = re.sub("Tissue_","",tissue)
    with open(log_file,"a+") as f:
        f.write("returning length of matrix: "+str(len(return_matrix)))
    return_matrix.index.name = "significant_variant_id"
    if qtl_type == "dapg_independent" or qtl_type == "conditional_independent":
        return_matrix.to_csv(output_files_directory+"QTLEnrich_matched_null_gwas_p_values_"+qtl_type+"_"+str(independent_ranking)+"_"+trait+"_"+tissue_removed+"_"+date+".tsv", index=True, header=True,sep="\t")
    else:
        return_matrix.to_csv(output_files_directory+"QTLEnrich_matched_null_gwas_p_values_"+qtl_type+"_"+trait+"_"+tissue_removed+"_"+date+".tsv", index=True, header=True,sep="\t")

def print_sampled_null_variants(log_file,output_files_directory,date,output_names_df,qtl_type,trait,tissue,independent_ranking):
    """
    returns sampled null variants
    """
    tissue_removed = re.sub("Tissue_","",tissue)

    with open(log_file,"a+") as f:
        f.write("returning length of null matrix: "+str(len(output_names_df)) + "\n")

    output_names_df.index.name = "significant_variant_id"

    if qtl_type == "dapg_independent" or qtl_type == "conditional_independent":
        output_names_df.to_csv(output_files_directory+"QTLEnrich_"+"null_variant_id_matrix_"+qtl_type+"_"+str(independent_ranking)+"_"+trait+"_"+tissue_removed+"_"+date+".tsv", index=True, header=True,sep="\t")
    else:
        output_names_df.to_csv(output_files_directory+"QTLEnrich_"+"null_variant_id_matrix_"+qtl_type+"_"+trait+"_"+tissue_removed+"_"+date+".tsv", index=True, header=True,sep="\t")

def compute_upper_bound_confidence_interval(adjusted_fold_enrichment,median_fold_enrichment,first_fold_enrichments):
    """
    compute upper bound confidence interval
    adjusted_fold_enrichment + adjusted_fold_enrichment*(upper 95 percentile)
    """
    u95 = (np.percentile(first_fold_enrichments,q=97.5) - median_fold_enrichment)/median_fold_enrichment
    upper_95 = adjusted_fold_enrichment + (adjusted_fold_enrichment*u95)
    return(round(upper_95,4))

def compute_lower_bound_confidence_interval(adjusted_fold_enrichment,median_fold_enrichment,first_fold_enrichments):
    """
    compute lower bound confidence interval
    adjusted fold-enrichment - (adjusted fold-enrichment* lower 95 percentile)
    """
    l95 = (median_fold_enrichment - np.percentile(first_fold_enrichments,q=2.5))/median_fold_enrichment
    lower_95 = adjusted_fold_enrichment - (adjusted_fold_enrichment*l95)
    return(round(lower_95,4))

def calculate_adjusted_fold_enrichment(num_permutations,fold_enrichments,observed_fold_enrichment):
    """
    Adjusted fold-enrichment: Observed fold enrichment / median fold enrichment of first 1000
    """
    if num_permutations <= 1000:
        use_folds = num_permutations
    else:
        use_folds = 1000
    first_fold_enrichments = fold_enrichments[:use_folds]
    fold_enrichments_median = np.median(first_fold_enrichments)

    adjusted_fold_enrichment = round(observed_fold_enrichment/fold_enrichments_median,4)
    return(first_fold_enrichments,adjusted_fold_enrichment,fold_enrichments_median)

def calculate_enrichment_pvalue(fold_enrichments,observed_fold_enrichment):
    """
    calculate the enrichment p-value 

    Enrichment p-value is denoted as the ratio of number of fold enrichments greater than the observed fold enrichment
    """
    enrichment_p_value = float(len([x for x in fold_enrichments if x >= observed_fold_enrichment]))/float(len(fold_enrichments))
    #to avoid an enrichment p-value of 0 after sampling
    if enrichment_p_value == 0.0:
        enrichment_p_value = 0.99/float(len(fold_enrichments))
        return("{:.2e}".format(Decimal(enrichment_p_value)))
    else:
        return("{:.2e}".format(Decimal(enrichment_p_value)))

def compute_true_trait(pi1,qtl_set):
    """
    # true trait associations = pi1*(# of qtl variants)
    """
    return(round(pi1*len(qtl_set)))

def compute_pi1(null_p_value_matrix,qtl_set,lambda_factor=0.8):
    """
    Inputs:
        1. Matrix containing p-values of sampled null variants
        2. QTLs for a given tissue with associated GWAS p-values
        3. Lambda factor cutoff: default 0.8

    pi_0 = proportion(p_gwas > lambda) / proportion(p_null > lambda)
    p_gwas: GWAS p-values of the significant QTLs
    p_null: GWAS p-values of matched SNPs.
    proportion(p_gwas > lambda) = (# of QTLs with a GWAS p-value above lambda)/(total number of QTLs)
    proportion(p_null > lambda) = (median number of matched null variants with a GWAS p-value above lambda)/(total number of QTLs)

    pi1 = (1-pi0)

    estimated # true trait associations amongst QTLs = pi1*(number of significant variants with a GWAS p-value)
    """
    #compute proportion(p_GWAS > lambda)
    #compute proportion of null GWAS of matched SNPs
    if lambda_factor == 0.0:
        proportion_gwas = (len(qtl_set.loc[qtl_set["gwas_p_value"] > lambda_factor]))/(len(qtl_set))
        null_proportion = np.median((null_p_value_matrix > lambda_factor).sum())/(len(null_p_value_matrix))
    else:
        proportion_gwas = (len(qtl_set.loc[qtl_set["gwas_p_value"] > lambda_factor]))/(len(qtl_set))/lambda_factor
        null_proportion = np.median((null_p_value_matrix > lambda_factor).sum())/(len(null_p_value_matrix))/lambda_factor

    #compute pi0
    if null_proportion != 0.0:
        pi0 = float(proportion_gwas)/float(null_proportion)     
        if pi0 > 1:
            pi0 = 1
            pi1 = 0
        else:
            pi1 = round((1-pi0),4) 
    else:
        #pi1 = np.nan
        pi1 = 0
    #compute estimated # true trait associations
    true_trait = compute_true_trait(pi1,qtl_set)
    return(pi1,true_trait)

def prepare_global_null_table(sub_table,null_table):
    """
    Inputs:
    1. sub_table: table consisting of variant and number of ld proxy variants (LD_Proxy) for each tissue and their corresponding maf and tss_distance
    2. null_table: table consisting of variants indicating if null per tissue (1 is null, 0 is not null)
    """
    null_df = null_table[["variant","Global","gwas_p_value"]].copy()
    null_df = null_df.loc[null_df["Global"] == 1]
    null_table = pd.merge(null_df,sub_table,on=["variant"])
    null_table["LD_Proxy"] = null_table["LD_Proxy"].fillna(0)
    null_table = null_table.dropna()

    return(null_table)

def prepare_null_table_sampling(sub_table,null_table,tissue):
    """
    Inputs:
        1. sub_table: table consisting of variant and number of ld proxy variants (LD_Proxy) for each tissue and their corresponding maf and tss_distance
        2. null_table: table consisting of variants indicating if null per tissue (1 is null, 0 is not null)
        3. Tissue: relevant tissue

    Selects null variants for a given tissue.  
    Removes extra columns.

    Returns null_table
    """
    null_df = null_table[["variant",tissue,"gwas_p_value"]].copy()
    null_df[tissue] = null_df[tissue].astype(np.uint8)
    null_df = null_df.loc[null_df[tissue] == 1]
    null_table = pd.merge(null_df,sub_table,on=["variant"])
    null_table["LD_Proxy"] = null_table["LD_Proxy"].fillna(0)
    null_table = null_table.dropna()
    return(null_table)

def prepare_significant_table_sampling(sub_table,significant_table,tissue,MAF,TSS,compute_tss_distance):
    """
    Inputs:
        1. sub_table: table consisting of variant and number of ld proxy variants (LD_Proxy) for each tissue and their corresponding maf and tss_distance
        2. significant_table: table consisting of variants indicating if null per tissue (1 is null, 0 is not null)
        3. Tissue: relevant tissue

    Selects significant variants for a given tissue.
    Removes extra columns.
    """
    if not compute_tss_distance:
        significant_table = pd.merge(sub_table,significant_table,on=["variant"],how="right")
        significant_table["LD_Proxy"] = significant_table["LD_Proxy"].fillna(0)
    else:
        sub_table = sub_table[["variant",MAF,"LD_Proxy"]]
        significant_table = pd.merge(sub_table,significant_table,on=["variant"],how="right")
        significant_table["LD_Proxy"] = significant_table["LD_Proxy"].fillna(0)
    return(significant_table)

def get_bins(df,col_name,num_quantiles=10.):
    """
    Find the percent of data in each quantile in order to appropriately bin data by bin number
    """
    num_quantiles = float(10.)
    percent_data_per_quantile = 1./num_quantiles
    return(df[col_name].rank(pct = True).apply(lambda x: x/percent_data_per_quantile).apply(np.ceil).astype(np.uint8))

def random_sampling_for_variant(variant, options_tables, variant_to_bins, num_permutations, use_num, p_value_cutoff):
    """
    Inner function to permute_null_variants_single
    equivalent to permutation for 1 variant
    """
    #Get nulls that are in the same bin as our variant
    permuted_null_variant_table = options_tables.get_group(variant_to_bins[variant])
    #If number of permutations > number of nulls, sample with replacement
    del options_tables
    del variant_to_bins
    if num_permutations >= len(permuted_null_variant_table):
        replace = True
    else:
        replace = False
    #Select rows for permutation, then add as columns to our temporary dfs
    picked_indices = permuted_null_variant_table.sample(n=num_permutations, replace=replace)
    del permuted_null_variant_table #dump memory
    temp_names_df = picked_indices["variant"].values[:use_num]
    temp_p_values_df = picked_indices["gwas_p_value"].values
    low_p_value_indices = np.where(temp_p_values_df <= p_value_cutoff)[0]
    temp_p_values_df = temp_p_values_df[:use_num]
    return(temp_names_df, temp_p_values_df, low_p_value_indices)

def permute_null_variants(log_file,date,output_files_directory,qtl_type,tissue,trait,observed_fold_enrichment, null_table, significant_table, options_tables, num_permutations, lower_threshold, upper_bound_permutations, p_value_cutoff,independent_ranking,keep_null_variants):
    """
    Internal function to rand_sample_null_match_confound to allow recursive calling
    """
    ranking=independent_ranking

    #Verify this makes sense
    if upper_bound_permutations < num_permutations:
        num_permutations = upper_bound_permutations
        with open(log_file,"a+") as f:
            f.write("You called an upper boundary of permutations smaller than the requested number of permutations."
              " Setting number of permutations to upper boundary.\n")
        print("You called an upper boundary of permutations smaller than the requested number of permutations."
              " Setting number of permutations to upper boundary.")

    start_time = time.time()
    with open(log_file,"a+") as f:
        f.write("Running with "+str(num_permutations)+" permutations for "+str(len(significant_table))+ " QTLs with an observed fold enrichment of "+str(observed_fold_enrichment) +"...\n")
    print("Running with "+str(num_permutations)+" permutations for "+str(len(significant_table))+ " QTLs with an observed fold enrichment of "+str(observed_fold_enrichment) +"...")

    #Loop through significant QTLs. Run permutation, generate output matrices
    columns  = [("P"+str(x)) for x in range(1,num_permutations+1)]

    temp_names_df = pd.DataFrame([]) #to speed up filling
    temp_p_values_df = pd.DataFrame([]) #to speed up filling

    #make a dict of tables to speed up permutation
    variant_to_bins = dict(zip(significant_table["variant"], significant_table["bin_id"])) #Dict of variant => bin
    #Remove the variants that did not fall into any bin
    variant_list = list(set(significant_table["variant"]))
    variant_list = [x for x in variant_list if "NA" not in variant_to_bins[x]]

    if len(variant_list) != len(significant_table["variant"]):
        difference = abs(len(variant_list) - len(significant_table["variant"]))
        with open(log_file,"a+") as f:
            f.write("Numbers do not align."+str(difference) + " missing variants\n")
    else:
        with open(log_file,"a+") as f:
            f.write("Numbers align. No missing variants.\n")

    stop_idx = len(variant_list)
    if num_permutations <= 1000:
       use_num = num_permutations
    else:
       use_num = 1000
    temp_names_df, temp_p_values_df, low_p_value_indices = map(np.array, zip(*[random_sampling_for_variant(variant, options_tables, variant_to_bins, num_permutations, use_num, p_value_cutoff) for variant in variant_list]))
    low_p_value_indices = np.array([item for sublist in low_p_value_indices for item in sublist]).astype(np.float32)

    temp_p_values_df = pd.DataFrame(data=temp_p_values_df).T #TODO: Remove this T and just write to output_p_values_df
    temp_p_values_df.columns = variant_list
    temp_names_df = pd.DataFrame(data=temp_names_df).T #TODO: Remove this T and just write to output_names_df
    temp_names_df.columns = variant_list

    output_names_df = temp_names_df.T #TODO: redundant transpose
    output_names_df.columns = columns[:use_num]
    output_p_values_df = temp_p_values_df.T #TODO: redundant transpose
    output_p_values_df.columns = columns[:use_num]

    #Find fold enrichment for each permutation    
    fold_timer = time.time()
    unique, fold_counter = np.unique(low_p_value_indices, return_counts=True)
#    fold_enrichments = np.array([(x/np.float32(len(output_p_values_df)))/p_value_cutoff for x in fold_counter])
    fold_enrichments = fold_counter/np.float32(len(variant_list))/p_value_cutoff
    with open(log_file,"a+") as f:
        f.write("--- %s seconds ---\n" % (np.round(time.time() - start_time)))
    print("--- %s seconds ---" % (time.time() - start_time))

    #Number of permuted fold enrichments larger than observed fold enrichment
    permutations_below_threshold = len([x for x in fold_enrichments if x >= observed_fold_enrichment])
    #If not enough null variants below threshold, repeat with num_permutations x 10
    if permutations_below_threshold < lower_threshold:
        with open(log_file,"a+") as f:
            f.write("lower threshold of permutations >= observed_fold_enrichment not reached.\nMultiplying number of permutations by 10.\n")
        print("lower threshold of permutations >= observed_fold_enrichment not reached.\nMultiplying number of permutations by 10.")
        num_permutations = num_permutations * 10
        #verify num_permutations has not exceeded maximum
        if num_permutations <= upper_bound_permutations:
            return(permute_null_variants(log_file,date,output_files_directory,qtl_type,tissue,trait,observed_fold_enrichment=observed_fold_enrichment, significant_table = significant_table, options_tables=options_tables,
                                                                 null_table = null_table, num_permutations = num_permutations,
                                                                 lower_threshold = lower_threshold, upper_bound_permutations = upper_bound_permutations,
                                                                 p_value_cutoff = p_value_cutoff,independent_ranking=ranking,keep_null_variants=keep_null_variants))
        else:
            with open(log_file,"a+") as f:
                f.write("upper boundary of number of permutations reached or exceeded. Returning p-value-matrix.\n")
            print("upper boundary of number of permutations reached or exceeded. Returning p-value-matrix.")
            # PRINT TABLE WE WANT PRINTED
            if keep_null_variants:
                print_sampled_null_variants(log_file,output_files_directory,date,output_names_df,qtl_type,trait,tissue,independent_ranking=ranking)  
            return(output_p_values_df, fold_enrichments)
    else:
        # PRINT TABLE WE WANT PRINTED
        if keep_null_variants:
            print_sampled_null_variants(log_file,output_files_directory,date,output_names_df,qtl_type,trait,tissue,independent_ranking=ranking)
        return(output_p_values_df, fold_enrichments)
    return(output_p_values_df, fold_enrichments)

def retrieve_tissue_columns(qtl_type,qtl_table,null_table,observed_fold_enrichments,null_option="tissue_dependent",interaction_qtl_dict=""):
    """
    obtain tissue columns from list of tissues or qtls being sampled
    """
    if null_option == "global":
        tissue_columns = [x for x in list(qtl_table.keys())]
        return(tissue_columns)
    elif qtl_type == "ieQTL" or qtl_type == "isQTL":
        tissue_columns = [x for x in interaction_qtl_dict.keys()]
        return(tissue_columns)
    else:
        tissue_columns = [x for x in list(null_table.columns.values) if x.startswith("Tissue_")]
        tissue_columns = [x for x in tissue_columns if x in observed_fold_enrichments]
        return(tissue_columns)

def identify_maf_tss_bin_headers(tissue,qtl_type,null_option="tissue_dependent",interaction_qtl_dict=""):
    """
    MAF = MAF+tissue, e.g. MAF_Nerve_Tibial
    TSS = TSS+tissue, e.g. TSS_Nerve_Tibial

    for global, 
    MAF = MAF_Global
    TSS = TSS_Global

    for interaction QTLs, there is a further expansion
    """
    if null_option == "global":
        MAF = "MAF_Global"
        TSS = "TSS_Global"
        return(tissue,MAF,TSS)
    elif (qtl_type == "ieQTL" or qtl_type == "isQTL") and null_option != "global":
        tissue_name = interaction_qtl_dict[tissue]
        MAF = "MAF_" + tissue_name
        TSS = "TSS_" + tissue_name
        return("Tissue_" + tissue_name,MAF,TSS)
    else:
        tissue_name = "_"+"_".join(tissue.split("_")[1:])
        MAF = "MAF" + tissue_name
        TSS = "TSS" + tissue_name
        return(tissue_name,MAF,TSS)

def reconfigure_null_table(null_table,sub_table,tissue,tissue_name,null_option="tissue_dependent"):
    """
    Selecting appropriate columns for null table
    """
    if null_option == "global":
        null_table_filtered = prepare_global_null_table(sub_table,null_table)
        return(null_table_filtered)
    else:
        null_table_filtered = prepare_null_table_sampling(sub_table,null_table,tissue)
        return(null_table_filtered)

def sample_match_null_confounders(log_file,date,output_files_directory,qtl_type,trait,confounders_table,null_table,qtl_table,observed_fold_enrichments, p_value_cutoff = np.float32(0.05),null_option="tissue_dependent",num_permutations=1000, upper_bound_permutations=100000, lower_threshold=5,num_quantiles=10,lambda_factor=0.5,independent_ranking=1,compute_tss_distance=False,interaction_qtl_dict="",keep_null_variants=False,keep_pvalue_matrix=False):
    """
    Inputs:
    in_table: table containing variant_ID, number of LD variants, minor allele frequency of variant, 
    transcription start site of variant (distance to TSS of eGene), 0 or 1 column for whether variant
    is a significant eQTL (to determine nulls), and GWAS p-value.
    
    lower_threshold: minimum number of permutations which should have a median p-value lower than cutoff. 
    If threshold is not met, permutation is repeated with size x 10

    fold_enrichment: value indicating fold-enrichment of GWAS p-values for significant QTLs
    num_quantiles: how many bins should each confounding factor be split into?
    num_permutations: sample size taken from null variants for each significant QTL
    upper_bound: max number of permutations (in case need increasing)
    """
    p_value_cutoff = np.float32(p_value_cutoff)
    ranking=independent_ranking

    print("Preparing data for permutations... This may take a few minutes.\n")
    with open(log_file,"a+") as f:
        f.write("Preparing data for permutations... This may take a few minutes.")
    start_time = time.time()

    #select tissue names
    tissue_columns = retrieve_tissue_columns(qtl_type,qtl_table,null_table,observed_fold_enrichments,null_option=null_option,interaction_qtl_dict=interaction_qtl_dict) 

    #declare empty dictionaries
    enrichment_p_value_dict = {}
    adjusted_fold_enrichment_dict = {}
    upper_bound_confidence_interval_dict = {}
    lower_bound_confidence_interval_dict = {}
    pi1_dict = {}
    true_trait_dict = {}

    for tissue in tissue_columns:
        with open(log_file,"a+") as f:
            f.write(tissue)
        print(tissue)
      
        tissue_name,MAF,TSS = identify_maf_tss_bin_headers(tissue,qtl_type,null_option=null_option,interaction_qtl_dict=interaction_qtl_dict)

        #intermediate table
        sub_table = confounders_table[["variant",MAF,TSS,"LD_Proxy"]].copy()

        #Touch up significant and null tables
        null_table_temp = reconfigure_null_table(null_table,sub_table,tissue,tissue_name,null_option=null_option)

        significant_table = prepare_significant_table_sampling(sub_table,qtl_table[tissue],tissue_name,MAF,TSS,compute_tss_distance) 

        del sub_table

        observed_fold_enrichment = observed_fold_enrichments[tissue]

        #Get bin values null table
        null_table_temp["MAF_bin"] = get_bins(null_table_temp, MAF,num_quantiles)
        null_table_temp["TSS_bin"] = get_bins(null_table_temp, TSS,num_quantiles,)
        null_table_temp["LD_Proxy_bin"] = get_bins(null_table_temp, "LD_Proxy",num_quantiles)
        null_table_temp["bin_id"] = null_table_temp["MAF_bin"].astype(str) + ":" + null_table_temp["TSS_bin"].astype(str) + ":" + null_table_temp["LD_Proxy_bin"].astype(str)

        #Find all bin options, then create a dict of bins
        options = list(set(null_table_temp["bin_id"]))
        options_tables = null_table_temp.groupby("bin_id")

        options = list((np.array(range(num_quantiles))+1)) #Get bin options
        #Generate bin_values for signif_table
        significant_table["MAF_bin"] = "NA"
        significant_table["TSS_bin"] = "NA"
        significant_table["LD_Proxy_bin"] = "NA"

        #Populate bins for signif_table
        for option in options:
            MAF_bin = null_table_temp.loc[null_table_temp["MAF_bin"] == option]
            TSS_bin = null_table_temp.loc[null_table_temp["TSS_bin"] == option]
            LD_Proxy_bin = null_table_temp.loc[null_table_temp["LD_Proxy_bin"] == option]

            str_option = str(option)

            if str_option == "1":
                significant_table.loc[(0 <= significant_table[MAF]) & (significant_table[MAF] <= MAF_bin[MAF].max()), ["MAF_bin"]] = option
                significant_table.loc[(significant_table[TSS].min() <= significant_table[TSS]) & (significant_table[TSS] <= TSS_bin[TSS].max()), ["TSS_bin"]] = option
                significant_table.loc[(0 <= significant_table["LD_Proxy"]) & (significant_table["LD_Proxy"] <= LD_Proxy_bin["LD_Proxy"].max()), ["LD_Proxy_bin"]] = option
            elif str_option == str(range(10)[-1]+1):
                significant_table.loc[(MAF_bin[MAF].min() <= significant_table[MAF]) & (significant_table[MAF] <= 1), ["MAF_bin"]] = option
                significant_table.loc[(TSS_bin[TSS].min() <= significant_table[TSS]) & (significant_table[TSS] <= significant_table[TSS].max()), ["TSS_bin"]] = option
                significant_table.loc[(LD_Proxy_bin["LD_Proxy"].min() <= significant_table["LD_Proxy"]) & (significant_table["LD_Proxy"] <= significant_table["LD_Proxy"].max()), ["LD_Proxy_bin"]] = option
            else:
                significant_table.loc[(MAF_bin[MAF].min() <= significant_table[MAF]) & (significant_table[MAF] <= MAF_bin[MAF].max()), ["MAF_bin"]] = option
                significant_table.loc[(TSS_bin[TSS].min() <= significant_table[TSS]) & (significant_table[TSS] <= TSS_bin[TSS].max()), ["TSS_bin"]] = option
                significant_table.loc[(LD_Proxy_bin["LD_Proxy"].min() <= significant_table["LD_Proxy"]) & (significant_table["LD_Proxy"] <= LD_Proxy_bin["LD_Proxy"].max()), ["LD_Proxy_bin"]] = option

        significant_table["bin_id"] = significant_table["MAF_bin"].astype(str) + ":" + significant_table["TSS_bin"].astype(str) + ":" + significant_table["LD_Proxy_bin"].astype(str)
        observed_fold_enrichment = observed_fold_enrichments[tissue]

        #permute null variants
        #return p_value matrix and fold_enrichments
        return_matrix, fold_enrichments= permute_null_variants(log_file,date,output_files_directory,qtl_type,tissue,trait,observed_fold_enrichment=observed_fold_enrichment,
                                                               significant_table = significant_table,options_tables=options_tables,null_table = null_table_temp,
                                                               num_permutations = num_permutations,lower_threshold = lower_threshold, upper_bound_permutations = upper_bound_permutations,
                                                               p_value_cutoff = p_value_cutoff,independent_ranking=ranking,keep_null_variants=keep_null_variants)

        #print p-value matrices
        if keep_pvalue_matrix:
            print_pvalue_matrix(log_file,output_files_directory,date,return_matrix,qtl_type,trait,tissue,independent_ranking=ranking)

        with open(log_file,"a+") as f:
            f.write("\n computing enrichment p-value, adjusted fold-enrichment, confidence intervals, and pi1")
        print("\n computing enrichment p-value, adjusted fold-enrichment, confidence intervals, and pi1")


        #compute enrichment p-value,adjusted fold-enrichment,confidence intervals,pi1,true traits
        #Enrichment p-value
        enrichment_p_value_dict[tissue] = calculate_enrichment_pvalue(fold_enrichments,observed_fold_enrichment)
        
        #Adjusted fold enrichment: Observed fold enrichment / median fold enrichment of first 1000
        first_fold_enrichments,adjusted_fold_enrichment_dict[tissue],fold_enrichments_median = calculate_adjusted_fold_enrichment(num_permutations,fold_enrichments,observed_fold_enrichment)

        #compute confidence intervals    
        upper_bound_confidence_interval_dict[tissue] = compute_upper_bound_confidence_interval(adjusted_fold_enrichment_dict[tissue],fold_enrichments_median,first_fold_enrichments)
        lower_bound_confidence_interval_dict[tissue] = compute_lower_bound_confidence_interval(adjusted_fold_enrichment_dict[tissue],fold_enrichments_median,first_fold_enrichments)

        #fill in return dicts        
        pi1_dict[tissue],true_trait_dict[tissue] = compute_pi1(return_matrix,significant_table,lambda_factor)

    return(enrichment_p_value_dict,adjusted_fold_enrichment_dict,upper_bound_confidence_interval_dict,lower_bound_confidence_interval_dict,pi1_dict,true_trait_dict)
    
