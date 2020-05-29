import sys
import os
import os.path
import datetime
import re
import pandas as pd
from os_sys_qtlenrich import make_directory
import qtl_gwas_parsing

def create_GeneEnrich_directory(output_directory):
    """
    creates directory to place GeneEnrich input files
    """
    GeneEnrich_directory_name = str(output_directory)+"GeneEnrich_Input_Files"
    GeneEnrich_directory = make_directory(GeneEnrich_directory_name)
    return(GeneEnrich_directory+"/")

def parse_significant_independent_variants(qtl_set,qtl_type,p=0.05,independent_ranking=1):
    """
    significant independent eQTLs: eQTL rank per eGene and below GWAS p-value
    """
    if qtl_type == "dapg_independent":
        significant_qtl_set = qtl_set.loc[qtl_set["ranking_by_beta"] == independent_ranking].copy()
    else:
        significant_qtl_set = qtl_set.loc[qtl_set["rank"] == independent_ranking].copy()
    significant_qtl_set = significant_qtl_set.loc[significant_qtl_set["gwas_p_value"] <= p]
    significant_qtl_set = significant_qtl_set.drop_duplicates("variant_id")
    significant_qtl_set = significant_qtl_set[["variant_id","gene_id"]]
    significant_qtl_set.columns = ["eVariant","eGene"]
    return(significant_qtl_set)

def parse_egenes_set(tissue_name,eGenes_directory,eGene_file_name,gwas,q=0.05,p=0.05,gencode="",subset_genes=False):
    """
    parses egenes for null genes from independent qtls
    """
    qtl_set = pd.read_csv(eGenes_directory+tissue_name+eGene_file_name,sep="\t")

    if "qval" in qtl_set:
        qtl_set = qtl_set[["gene_id","variant_id","qval"]]
        qtl_set = qtl_set.loc[qtl_set["qval"] > q]
    else:
        qtl_set = qtl_set[["gene_id","variant_id"]]

    #split chromosome and position
    #in case to match on chrosome and position instead of variant
    qtl_set_rename = qtl_gwas_parsing.rename_variant(qtl_set)
    qtl_set_split = qtl_gwas_parsing.split_chromosome_position(qtl_set_rename)

    qtl_set_split["gene_id"] = qtl_set_split["gene_id"].str.extract(r".*(ENSG\d+)")

    #subset on protein coding and lincRNAs
    if subset_genes:
        qtl_set_split = filter_protein_coding_lincRNAs(gencode,qtl_set_split)

    #merge with gwas
    if "variant" in gwas:
        qtl_gwas = pd.merge(qtl_set_split,gwas, on = ["variant"])
    else:
        qtl_gwas = pd.merge(qtl_set_split,gwas, on = ["chr","pos"])
    qtl_gwas = qtl_gwas.drop(["chr","pos","Allele"],axis=1)
    qtl_gwas = qtl_gwas.drop_duplicates(["variant"],keep="first")

    #p-value greater than 0.05
    qtl_gwas = qtl_gwas.loc[qtl_gwas["gwas_p_value"] > p]
    qtl_gwas = qtl_gwas[["gene_id"]]
    return(qtl_gwas)

def parse_non_significant_egenes(qtl_list,p,independent_ranking,qtl_type):
    """
    non significant independent genes: genes at a given rank with gwas p-value above cutoff
    """
    if qtl_type == "dapg_independent":
        independent_non_significant_egenes = qtl_list.loc[qtl_list["ranking_by_beta"] == independent_ranking].copy()
    else:
        independent_non_significant_egenes = qtl_list.loc[qtl_list["rank"] == independent_ranking].copy()    

    independent_non_significant_egenes = independent_non_significant_egenes.loc[independent_non_significant_egenes["gwas_p_value"] > p]
    independent_non_significant_egenes = independent_non_significant_egenes[["gene_id"]]
    return(independent_non_significant_egenes)


def parse_null_independent_variants(qtls,tissue_name,eGenes_directory,eGene_file_name,gwas,qtl_type,p,q,independent_ranking,gencode,subset_genes):
    """
    parses null genes for independent eQTLs

    non-significant genes (q > cutoff) whose best QTL has GWAS p-value above a given cutoff and significant eGenes (q < cutoff) whose    independent qtl has gwas p-value above cutoff 
    """
    non_significant_egenes = parse_egenes_set(tissue_name,eGenes_directory,eGene_file_name,gwas,q=q,p=p,gencode=gencode,subset_genes=subset_genes)   

    independent_non_significant_egenes = parse_non_significant_egenes(qtls,p,independent_ranking,qtl_type) 

    null_qtl_set = non_significant_egenes.append(independent_non_significant_egenes)
    null_qtl_set = null_qtl_set.drop_duplicates("gene_id",keep="first")
    null_qtl_set.column = ["eGene"]
    return(null_qtl_set)

def parse_significant_variants(qtl_set,q=0.05,p=0.05):
    """
    significant genes for GeneEnrich Inputs for best eQTL or best sQTL
    """
    if "qval" in qtl_set:
        qtl_set = qtl_set.loc[(qtl_set["qval"]<=q) & (qtl_set["gwas_p_value"]<=p)].copy()

    qtl_set = qtl_set.loc[qtl_set["gwas_p_value"]<=p]

    qtl_set = qtl_set[["variant","gene_id"]]
    qtl_set.columns = ["eVariant","eGene"]
    return(qtl_set)

def parse_null_variants(qtl_set,significant_qtl_set,q=0.05,p=0.05):
    """
    null variants for GeneEnrich Inputs for best eQTL or best sQTL
    defined as: genes with q-value above significant cutoff whose best QTL has GWAS p-value above cutoff
    """
    if "qval" in qtl_set:
        qtl_set = qtl_set.loc[qtl_set["qval"] > q].copy()

    qtl_set = qtl_set.loc[qtl_set["gwas_p_value"]>p].copy()
    qtl_set = qtl_set[["gene_id"]].copy()
    qtl_set.columns = ["eGene"]
    qtl_set = qtl_set.loc[~qtl_set["eGene"].isin(significant_qtl_set["eGene"])]
    qtl_set = qtl_set.drop_duplicates(["eGene"],keep="first")
    return(qtl_set)

def create_GeneEnrich_inputs_best_eqtl_sqtl(trait,output_directory,date,qtl_set_gwas_dict,qtl_type,exp_label,q=0.05,p=0.05):
    """
    prints input files for GeneEnrich   
    """
    #get directory
    GeneEnrich_directory = create_GeneEnrich_directory(output_directory)

    #get tissue names
    #keys of dictionary
    dict_keys = [x for x in qtl_set_gwas_dict.keys()]

    for tissue in dict_keys:
        significant_qtl_set = parse_significant_variants(qtl_set_gwas_dict[tissue],q,p)
        null_qtl_set = parse_null_variants(qtl_set_gwas_dict[tissue],significant_qtl_set,q,p)

        significant_qtl_set.to_csv(GeneEnrich_directory+"GeneEnrich_input_"+qtl_type+"_"+tissue+"_"+trait+"_"+exp_label+"_significant.txt",index=None,header=True,sep="\t")
        null_qtl_set.to_csv(GeneEnrich_directory+"GeneEnrich_input_"+qtl_type+"_"+tissue+"_"+trait+"_"+exp_label+"_null.txt",index=None,header=True,sep="\t")

def create_GeneEnrich_inputs_independent_qtls(eGenes_directory,eGene_file_name,gwas,gwas_trait,output_directory,date,GeneEnrich_input_gwas_dict,qtl_type,p=0.05,q=0.05,independent_ranking=1,gencode="",subset_genes=False):
    """
    
    """
    #get directory
    GeneEnrich_directory = create_GeneEnrich_directory(output_directory)  

    #get tissue names
    #keys of dictionary
    dict_keys = [x for x in GeneEnrich_input_gwas_dict.keys()]

    for tissue in dict_keys:
        significant_qtl_set = parse_significant_independent_variants(GeneEnrich_input_gwas_dict[tissue],qtl_type,p,independent_ranking)

        tissue_name = re.sub("Tissue_","",tissue)
        null_qtl_set = parse_null_independent_variants(GeneEnrich_input_gwas_dict[tissue],tissue_name,eGenes_directory,eGene_file_name,gwas,qtl_type,p,q,independent_ranking,gencode=gencode,subset_genes=subset_genes)


        significant_qtl_set.to_csv(GeneEnrich_directory+"GeneEnrich_input_"+qtl_type+"_ranking_"+str(independent_ranking)+"_"+tissue+"_"+trait+"_significant.txt",index=none,header=true,sep="\t")
        null_qtl_set.to_csv(GeneEnrich_directory+"GeneEnrich_input_"+qtl_type+"_ranking_"+str(independent_ranking)+"_"+tissue+"_"+trait+"_null.txt",index=none,header=true,sep="\t")





