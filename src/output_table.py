
"""
This script consists of functions related to computing the QTLEnrich results output table
"""

import pandas as pd
import numpy as np

def estimated_true_trait(output):
    """
    Compute estimated num true trait associations

    equation: ((adjusted fold-enrichment - 1)/adjusted fold-enrichment * OBS)
    """
    output["Estimated_num_trait_associations_with_GWAS_p<0.05"] = ((output["Adjusted_Fold_Enrichment"]-1)/output["Adjusted_Fold_Enrichment"])*(output["Num_OBS_QTLs_with_GWAS_p<0.05"])
    output["Estimated_num_trait_associations_with_GWAS_p<0.05"] = output["Estimated_num_trait_associations_with_GWAS_p<0.05"].astype(int)
    return(output)

def expected_trait(output):
    """
    Compute expected num trait associations based on adjusted fold-enrichment

    equation: OBS * (1-((adjusted fold-enrichment - 1)/adjusted fold-enrichment))
    """
    output["Num_EXP_QTLs_with_GWAS_p<0.05"] = (output["Num_OBS_QTLs_with_GWAS_p<0.05"])*(1-((output["Adjusted_Fold_Enrichment"]-1)/output["Adjusted_Fold_Enrichment"]))
    output["Num_EXP_QTLs_with_GWAS_p<0.05"] = output["Num_EXP_QTLs_with_GWAS_p<0.05"].astype(int)
    return(output)

def dict_to_df(dict,col_name):
    """
    Converts dictionary into a dataframe
    """
    return(pd.DataFrame(dict.items(),columns = ["Tissue_(QTL)",col_name]))

def trait_dict_to_df(trait_dict):
    """
    Special case of converting dictionary to dataframe
    """
    return(pd.DataFrame(trait_dict).melt().rename(columns = {"variable":"Trait_(GWAS)","value":"Tissue_(QTL)"}))

def sort_output_table(output):
    """
    Take QTLEnrich output table sort by enrichment p-value and adjusted fold-enrichment.

    Enrichment p-value is sorted in ascending order.
    Adjusted fold-enrichment is sorted in descending order.
    """
    output["intermediate_columns"] = output["Enrichment_P_value"].apply(lambda x: np.float32(x))
    output = output.sort_values(["intermediate_columns","Adjusted_Fold_Enrichment"],ascending=[True,False])
    output = output.drop(["intermediate_columns"],axis=1)
    return(output)

def merge_tables(df_list):
    """
    Merges list of dataframes into one dataframe
    """
    output = df_list[0]
    for df_ in df_list[1:]:
        output = output.merge(df_, on="Tissue_(QTL)")
    output["Tissue_(QTL)"] = output["Tissue_(QTL)"].str.replace("Tissue_","")
    return(output)

def create_output_table(trait_dict,original_length_dict,best_QTL_variants_dict,observed_fold_enrichment,obs_dict,enrichment_p_value_dict, adjusted_fold_enrichment_dict,upper_bound_confidence_interval_dict,lower_bound_confidence_interval_dict,pi_1_dict,true_trait_dict):

    #make dictionary inputs into dataframes
    trait_dict_df = trait_dict_to_df(trait_dict)
    original_length_dict_df = dict_to_df(original_length_dict,col_name="Num_QTL_variants")
    best_QTL_variants_dict_df = dict_to_df(best_QTL_variants_dict,col_name="Num_QTL_variants_in_GWAS")
    obs_dict_df = dict_to_df(obs_dict,col_name="Num_OBS_QTLs_with_GWAS_p<0.05")
    observed_fold_enrichment_df = dict_to_df(observed_fold_enrichment,col_name="Fold-Enrichment_(OBS/EXP)")
    adjusted_fold_enrichment_dict_df = dict_to_df(adjusted_fold_enrichment_dict,col_name="Adjusted_Fold_Enrichment")
    lower_bound_confidence_interval_dict_df = dict_to_df(lower_bound_confidence_interval_dict,col_name="Lower_bound_95%_CI")
    upper_bound_confidence_interval_dict_df = dict_to_df(upper_bound_confidence_interval_dict,col_name="Upper_bound_95%_CI")
    enrichment_p_value_dict_df = dict_to_df(enrichment_p_value_dict,col_name="Enrichment_P_value")
    pi1_dict_df = dict_to_df(pi_1_dict,col_name="Empirical_Pi1")
    true_trait_dict_df = dict_to_df(true_trait_dict,col_name="Estimated_total_num_trait_associations")

    df_list = [trait_dict_df,original_length_dict_df,best_QTL_variants_dict_df,
               obs_dict_df,observed_fold_enrichment_df,adjusted_fold_enrichment_dict_df,
               lower_bound_confidence_interval_dict_df,upper_bound_confidence_interval_dict_df,enrichment_p_value_dict_df,pi1_dict_df,true_trait_dict_df]

    output_table = merge_tables(df_list)

    output_table = estimated_true_trait(output_table)
    output_table = expected_trait(output_table)

    #reorder columns
    output_table = output_table[["Trait_(GWAS)","Tissue_(QTL)","Num_QTL_variants","Num_QTL_variants_in_GWAS",
                                 "Num_OBS_QTLs_with_GWAS_p<0.05","Num_EXP_QTLs_with_GWAS_p<0.05","Estimated_num_trait_associations_with_GWAS_p<0.05",
                                 "Empirical_Pi1","Estimated_total_num_trait_associations","Fold-Enrichment_(OBS/EXP)",
                                 "Adjusted_Fold_Enrichment","Lower_bound_95%_CI","Upper_bound_95%_CI","Enrichment_P_value"]].copy()

    output_table = sort_output_table(output_table) 

    return(output_table)

