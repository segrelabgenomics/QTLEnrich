
"""
This script computes fold enrichments for each tissue

"""


from __future__ import division
import pandas as pd

def compute_obs(qtl_set,p_value=0.05):
    """
    observed number of qtls below a gwas cutoff
    """
    x = int(len(qtl_set[qtl_set['gwas_p_value'] <= p_value]))
    obs=x
    return(obs,x)

def compute_exp(qtl_set,p_value=0.05):
    """
    expected number of qtls below a gwas cutoff
    exp = (# qtls)*p
    """
    x = int(round(len(qtl_set)*p_value))
    exp=x
    return(exp,x)

def compute_fold(obs,exp):
    """
    Compute observed fold enrichment
    Fold = obs/exp
    """
    if exp != 0:
        return(round(obs/exp,4))
    return(0)
    
def compute_qtl_statistics(gwas_trait,significant_qtl_dict,significant_qtl_gwas_dict,p_value=0.05):

    """
    Computes following statistics: 
        1.trait-tissue pair
        2.length of qtls
        3.length of qtls with gwas cutoff
        4.obs
        5.exp
        6.fold-enrichment (obs/exp)
    """
    #set empty dictionaries for each statistic
    trait_dict = {}
    best_eqtl_variants_dict = {}
    observed_fold = {}
    obs_dict = {}
    exp_dict = {}
    original_length_dict = {}

    #extract keys from dictionary
    dict_keys = [x for x in significant_qtl_dict.keys()] 

    gwas_trait_list = [gwas_trait]

    for i in gwas_trait_list:
        trait_dict[i]=dict_keys

    #loop through tissues and calculate statistics
    for i in dict_keys:
        original_length_dict[i] = len(significant_qtl_dict[i])
        best_eqtl_variants_dict[i] = len(significant_qtl_gwas_dict[i])

        obs,obs_dict[i] = compute_obs(significant_qtl_gwas_dict[i],p_value)
        exp,exp_dict[i] = compute_exp(significant_qtl_gwas_dict[i],p_value)

        observed_fold[i] = compute_fold(obs,exp)

    return(trait_dict,original_length_dict,observed_fold,obs_dict,best_eqtl_variants_dict)






