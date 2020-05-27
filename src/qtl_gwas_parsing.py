import os
import os.path
import re
import pandas as pd
import numpy as np
from os_sys_qtlenrich import make_directory

def extract_qtl_files(qtl_directory,qtl_file_name):
    """
    returns list of files from directory
    """
    return([file for file in os.listdir(qtl_directory) if file.endswith(qtl_file_name)])

def interaction_qtl_dictionary(qtl_type,qtl_file_name,qtl_files):
    """
    due to the nomenclature of interaction qtls, e.g. whole_blood.monocytes,
    an additional dictionary is defined
    """
    if qtl_type == "ieqtl" or qtl_type == "isqtl":
        parsed_qtl_files = [re.sub(qtl_file_name,"",y) for y in qtl_files]
        return({key:key.split(".",1)[0] for key in parsed_qtl_files})
    return("")

def declare_gencode(args):
    """ 
    format gencode file or return empty string
    """
    if args.subset_genes or args.compute_tss_distance:
        gencode = format_gencode(args.gencode_file)
        return(gencode)
    return("")

def extract_gwas_trait_name(gwas_file,gwas_trait):
    """
    gwas trait is name of file without suffix.
    e.g. type_2_diabetes.txt -> type_2_diabetes
    """
    if gwas_trait is None:
        trait_name = re.sub(".txt.*$","", os.path.basename(gwas_file))
        trait_name = re.sub(".csv.*$","",trait_name)
        trait_name = re.sub(".tsv.*$","",trait_name)
        return(trait_name)
    return(gwas_trait)

def generate_interaction_qtl_dictionary(qtl_type,qtl_files,file_name):
    """
    Interaction qtls (ieqtl or isqtl) are formatted differently, where they include the tissue subtype.
    e.g. instead of Whole_Blood -> Whole_Blood.keratinocytes
    """
    if (qtl_type == "ieqtl") or (qtl_type == "isqtl"):
        parsed_qtl_files = [re.sub(file_name,"",y) for y in qtl_files]
        return({key:key.split(".",1)[0] for key in parsed_qtl_files})
    return("")

def check_if_column_headers_present(gwas):
    """
    Ensure gwas contains appropriate headers
    """
    for i in gwas.columns:
        if i.lower() in ("gwas_p_value","p","p_value","p-value","pvalue","pval","p_val","p-val"):
            gwas = gwas.rename(columns = {i:"gwas_p_value"})

        if i.lower() in ("ref","reference","non_effect_allele","non-effect-allele"):
            gwas = gwas.rename(columns = {i:"non_effect_allele"})

        if i.lower() in ("alt","alternate","effect_allele","effect-allele"):
            gwas = gwas.rename(columns = {i:"effect_allele"})

        if i.lower() in ("chr","chromosome","chrom"):
            gwas = gwas.rename(columns = {i:"chr"})

        if i.lower() in ("pos","position"):
            gwas = gwas.rename(columns = {i:"pos"})
    return(gwas)

def create_gwas_variant_id(gwas,build_version):
    """
    if non_effect_allele and effect_allele are present in gwas, create variant id.
    matching with qtl variants will be done 
    """
    #Create gwas Variant ID if non_effect_allele and effect_allele are found in gwas
    if "non_effect_allele" in gwas and "effect_allele" in gwas:
        gwas["variant"] = gwas["chr"].astype(str) + "_" + gwas["pos"].astype(str) + "_" + gwas["non_effect_allele"].astype(str) + "_" + gwas["effect_allele"].astype(str) + "_" + str(build_version)
        gwas = gwas[["variant","gwas_p_value"]]
        return(gwas)

    #if no ref/alt allele
    else:
        gwas = gwas[["chr","pos","gwas_p_value"]]
        return(gwas)

def check_gwas_header_datatypes(gwas):
    """
    chromosome must be "chr"+N, e.g. chr4
    position must be int, e.g. 153920
    """
    if not all(gwas["chr"].str.startswith("chr")):
        gwas["chr"] = gwas["chr"].apply(lambda x: "chr"+str(x) if not x.startswith("chr") else str(x))
    gwas["pos"] = gwas["pos"].astype(int)
    return(gwas)

def prepare_gwas_file(gwas_file,build_version):
    """
    Clean gwas Dataset to account for differences in 
    notation of chromosomes and select relevant columns.
    Position will be treated as integers
    """
    gwas = pd.read_csv(gwas_file,sep="\t")
    #loop through column names in gwas dataset
    #change column names
    gwas = check_if_column_headers_present(gwas)

    #accomodate for extrememly low p-values
    gwas["gwas_p_value"] = gwas["gwas_p_value"].replace(0.0,1e-280)

    #Ensure chromosome column indentifies as string datatype
    gwas = check_gwas_header_datatypes(gwas)

    gwas = create_gwas_variant_id(gwas,build_version)
    return(gwas)

def prepare_null_table(null_table,gwas):
    """
    read null table
    split variant to chr,pos,allele
 
    merge with gwas
    """
    null_table_df = pd.read_csv(null_table,sep="\t")
    null_table_split = split_chromosome_position(null_table_df)

    if "variant" in gwas:
        null_table_gwas = pd.merge(null_table_split,gwas, on = ["variant"])
    else:
        null_table_gwas = pd.merge(null_table_split,gwas, on = ["chr","pos"])
    null_table_gwas = null_table_gwas.drop(["chr","pos","Allele"],axis=1).copy()
    null_table_gwas = null_table_gwas.drop_duplicates(["variant"],keep="first")

    del null_table_df,null_table_split

    return(null_table_gwas)

def read_confounders_table(confounders_table):
    """
    read in confounders table
    """
    confounders_table_df = pd.read_csv(confounders_table,sep="\t")
    return(confounders_table_df)

def rename_variant(df):
    """
    renames variant_id to variant
    """
    return(df.rename(columns={"variant_id":"variant"}))

def split_chromosome_position(df):
    """
    Takes variant column of dataframe and splits into chromosome, position, and allele
    """
    df["chr"],df["pos"],df["Allele"] = df["variant"].str.split("_",2).str
    df["pos"] = df["pos"].astype(int)
    return(df)

def format_gencode(gencode):
    """
    Prepare gencode file to merge with retina eqtl file
    """
    gencode_df = pd.read_csv(gencode,sep="\t")

    #select genes
    gencode_df_gene = gencode_df.loc[gencode_df["feature"] == "gene"].copy()

    #extract relevant columns
    gencode_df_gene = gencode_df_gene[["gene_id","start","end","strand","gene_type"]]
    gencode_df_gene["gene_id"] = gencode_df_gene["gene_id"].str.extract(r".*(ENSG\d+)")
    gencode_df_gene = gencode_df_gene.rename(columns = {"start":"gene_start","end":"gene_end"})

    return(gencode_df_gene)

def filter_protein_coding_lincRNAs(gencode,qtl_set):
    """
    filter qtls based on protein coding and linc RNA genes
    """
    gencode_subset = gencode.loc[(gencode["gene_type"] == "lincRNA") | (gencode["gene_type"] == "protein_coding")].copy()
    qtl_set = qtl_set.loc[qtl_set["gene_id"].isin(gencode_subset["gene_id"])].copy()
    del gencode_subset
    return(qtl_set)

def merge_gencode_qtl(qtl_set,gencode):
    """
    merge gencode file and retina file on ensembl gene_id, without decimal suffix
    """
    print("merging")
    merged_file = pd.merge(qtl_set,gencode,on=["gene_id"])
    return(merged_file)

def perform_tss_distance_computation(qtl_set,tissue,gencode,null_option):
    """
    compute tss distance

    if strand is positive:
        tss distance = start of variant -- start of gene
    if strand is negative:
        tss distance = (start of variant -- end of gene)*(-1)
    """
    gencode_qtl = merge_gencode_qtl(qtl_set,gencode)

    #compute TSS Distance for positive and negative strands separately
    #append back together
    gencode_qtl_pos_strand = gencode_qtl.loc[gencode_qtl["strand"] == "+"].copy()
    gencode_qtl_neg_strand = gencode_qtl.loc[gencode_qtl["strand"] == "-"].copy()

    if null_option == "tissue_dependent":
        gencode_qtl_pos_strand["TSS_"+tissue] = gencode_qtl_pos_strand["pos"] - gencode_qtl_pos_strand["gene_start"]
        gencode_qtl_neg_strand["TSS_"+tissue] = (gencode_qtl_neg_strand["pos"] - gencode_qtl_neg_strand["gene_end"])*(-1)
    else:
        gencode_qtl_pos_strand["TSS_Global"] = gencode_qtl_pos_strand["pos"] - gencode_qtl_pos_strand["gene_start"]
        gencode_qtl_neg_strand["TSS_Global"] = (gencode_qtl_neg_strand["pos"] - gencode_qtl_neg_strand["gene_end"])*(-1)

    gencode_qtl_full = gencode_qtl_pos_strand.append(gencode_qtl_neg_strand)
    gencode_qtl_full = gencode_qtl_full.drop(["gene_start","gene_end","strand"],axis=1)

    return(gencode_qtl_full)

def parse_interaction_qtls(qtl_set,tissue,compute_tss_distance=False,gencode="",subset_genes=False,null_option="tissue_dependent"):
    """
    parses qtls for ieqtls per eGene and isqtls
    If pval_adj_bh is not present, then qval column is not taken.
    merges with gwas trait of interest to obtain associated gwas p_value
    """
    if "pval_adj_bh" in qtl_set:
        qtl_set = qtl_set[["gene_id","variant_id","pval_adj_bh"]]
    else:
        qtl_set = qtl_set[["gene_id","variant_id"]]

    qtl_set_rename = rename_variant(qtl_set)
    qtl_set_split = split_chromosome_position(qtl_set_rename)

    #subset on protein coding and lincRNAs
    if subset_genes:
        qtl_set_split = filter_protein_coding_lincRNAs(gencode,qtl_set_split)

    if compute_tss_distance:
        print("computing tss distance. Checking for gencode file.")
        qtl_set_split = perform_tss_distance_computation(qtl_set_split,tissue,gencode,null_option)

    return(qtl_set_split,qtl_set_split)

def parse_best_splice_qtls(qtl_set,tissue,compute_tss_distance=False,gencode="",subset_genes=False,null_option="tissue_dependent"):
    """
    parses qtls for best eqtl per eGene and sqtls

    If qvalue is set to None, then qval column is not taken.

    merges with gwas trait of interest to obtain associated gwas p_value
    """
    if "qval" in qtl_set:
        qtl_set = qtl_set[["gene_id","variant_id","qval"]]
    else:
        qtl_set = qtl_set[["gene_id","variant_id"]]

    #split chromosome and position
    #in case to match on chrosome and position instead of variant
    qtl_set_rename = rename_variant(qtl_set)
    qtl_set_split = split_chromosome_position(qtl_set_rename)

    qtl_set_split["gene_id"] = qtl_set_split["gene_id"].str.extract(r".*(ENSG\d+.*)")

    #subset on protein coding and lincRNAs
    if subset_genes:
        qtl_set_split = filter_protein_coding_lincRNAs(gencode,qtl_set_split)

    if compute_tss_distance:
        print("computing tss distance. Checking for gencode file.")
        qtl_set_split = perform_tss_distance_computation(qtl_set_split,tissue,gencode,null_option)
    return(qtl_set_split,qtl_set_split)

def parse_independent_eqtls(qtl_file,tissue,qtl_type,compute_tss_distance=False,gencode="",subset_genes=False,null_option="tissue_dependent"):
    """
    parses independent eqtls.

    merges with gwas trait of interest to obtain associated gwas p_value.
    """
    if qtl_type == "dapg_independent":
        qtl_set = qtl_file[["gene_id","variant_id","ranking_by_beta"]]
    else:
        qtl_set = qtl_file[["gene_id","variant_id","rank"]]

    qtl_set_rename = rename_variant(qtl_set)
    qtl_set_split = split_chromosome_position(qtl_set_rename)

    qtl_set_split["gene_id"] = qtl_set_split["gene_id"].str.extract(r".*(ENSG\d+.*)")

    #subset on protein coding and lincRNAs
    if subset_genes:
        qtl_set_split = filter_protein_coding_lincRNAs(gencode,qtl_set_split)

    if compute_tss_distance:
        qtl_set_split = perform_tss_distance_computation(qtl_set_split,tissue,gencode,null_option)

    return(qtl_set_split,qtl_set_split)

def find_min_tss_distance(qtl,tissue,null_option):
    """
    find min tss distance
    """
    if null_option == "tissue_dependent":
        qtl["abs_tss_distance"] = qtl["TSS_"+tissue].abs()
    else:
        qtl["abs_tss_distance"] = qtl["TSS_Global"].abs()

    #for gene
    qtl = qtl.loc[qtl.groupby(["gene_id"])["abs_tss_distance"].idxmin()]
    #for variant
    qtl = qtl.loc[qtl.groupby(["variant"])["abs_tss_distance"].idxmin()]
    qtl = qtl.drop(["abs_tss_distance"],axis=1)
    return(qtl)

def extract_significant_interaction_qtls(qtl,tissue,q_value,compute_tss_distance):
    """
    Significant ieqtl or isqtl is defined by a pval_cutoff.
    If no qvalue is present, all variant gene pairs will be passed.
    """
    q_value = np.float32(q_value)
    if "pval_adj_bh" in qtl:
        qtl["pval_adj_bh"] = qtl["pval_adj_bh"].astype(float)
        qtl = qtl.loc[qtl["pval_adj_bh"]<=q_value].copy()

    if compute_tss_distance:
        qtl = find_min_tss_distance(qtl,tissue)

    return(qtl)

def extract_significant_best_splice_qtls(qtl,tissue,q_value,compute_tss_distance,null_option):
    """
    Significant best eqtl or sqtl is defined by a q_value cutoff.
    If no qvalue is present, all variant gene pairs will be passed.
    """
    q_value = np.float32(q_value)
    if "qval" in qtl:
        qtl["qval"] = qtl["qval"].astype(float)
        qtl = qtl.loc[qtl["qval"]<=q_value].copy()

    if compute_tss_distance:
        qtl = find_min_tss_distance(qtl,tissue,null_option)

    return(qtl)

def extract_significant_independent_qtls(qtl,tissue,independent_ranking,qtl_type,compute_tss_distance):
    """
    Significant independent eqtl is defined as eqtls at a given ranking.
    """
    if qtl_type == "dapg_independent":    
        significant_qtl = qtl[qtl["ranking_by_beta"] == independent_ranking].copy()
    else:
        significant_qtl = qtl[qtl["rank"] == independent_ranking].copy()

    if compute_tss_distance:
        qtl = find_min_tss_distance(qtl,tissue)

    return(significant_qtl)

def print_gwas_pvalues(date,qtl,tissue,qtl_type,trait,output_files_directory,independent_ranking=1):
    """
    prints table containing gwas p-values for significant QTLs
    """
    qtl_df = qtl[["variant","gwas_p_value"]]

    if qtl_type == "dapg_independent" or qtl_type == "conditional_independent":
        qtl_df.to_csv(output_files_directory+"QTLEnrich_significant_gwas_p_values_"+qtl_type+"_"+str(independent_ranking)+"_"+trait+"_"+tissue+"_"+date+".txt", index=False, header=True,sep="\t")
    else:
        qtl_df.to_csv(output_files_directory+"QTLEnrich_significant_gwas_p_values_"+qtl_type+"_"+trait+"_"+tissue+"_"+date+".txt", index=False, header=True,sep="\t")

def merge_with_gwas(date,gwas_file,qtl,tissue,qtl_type,trait,output_files_directory,independent_ranking):
    """
    Matches variants from qtls with gwas file.

    If ref and alt are found in gwas, matches on variant id.

    If not, matches on chromosome and position pair.
    """
    if "variant" in gwas_file:
        qtl_gwas = pd.merge(qtl,gwas_file, on = ["variant"])
    else:
        qtl_gwas = pd.merge(qtl,gwas_file, on = ["chr","pos"])
    qtl_gwas = qtl_gwas.drop(["chr","pos","Allele"],axis=1)
    qtl_gwas = qtl_gwas.drop_duplicates(["variant"],keep="first")

    #print pvalues for significant
    print_gwas_pvalues(date,qtl_gwas,tissue,qtl_type,trait,output_files_directory,independent_ranking)
    return(qtl_gwas)

def parse_qtls(output_files_directory,date,trait,qtl_files,qtl_directory,qtl_type,file_name,gwas,q_value=0.05,independent_ranking=1,compute_tss_distance=False,gencode="",subset_genes=False,null_option="tissue_dependent"):
    """
    parses qtl files before performing random sampling on null variants from matched confounder bins.

    Expanded for additional types of  qtls

    qtl types are:
        1. "best_eqtl": best eqtl per eGene
        2. "splice": splice qtls
        3. "dapg_independent","conditional_independent"

    "best_eqtl" or "splice"
    if "best" or "splice" is chosen: qval must be indicated (default q-value is 0.05)
    if user does not want to subset on qval, qval argument must be set to "False"

    "independent"
    if independent is chosen, ranking must be indicated (default ranking is 1)

    Removes duplicates on variant id.  If qval, selects lowest qval
    """
    #set empty dictionary
    qtl_dict={}
    geneenrich_input_dict = {}
    geneenrich_input_gwas_dict = {}
    significant_qtl_dict = {}
    significant_qtl_gwas_dict = {}

    #loop through each file
    #select relevant columns 
    for i in qtl_files:
        #collect tissue name
        tissue = re.sub(file_name,"",i)
        tissue_name = "Tissue_"+str(tissue)

        print(tissue)

        qtl_file = pd.read_csv(qtl_directory+i,sep="\t")
        #parse qtl by type
        #best/splice
        if qtl_type == "best_eqtl" or qtl_type == "best_sqtl":
            geneenrich_input_dict[tissue_name],qtl_dict[tissue_name] = parse_best_splice_qtls(qtl_file,tissue,compute_tss_distance,gencode,subset_genes,null_option)

            geneenrich_input_gwas_dict[tissue_name] = merge_with_gwas(date,gwas,geneenrich_input_dict[tissue_name],tissue,qtl_type,trait,output_files_directory,independent_ranking)

            significant_qtl_dict[tissue_name]=extract_significant_best_splice_qtls(qtl_dict[tissue_name],tissue,q_value,compute_tss_distance,null_option)

            significant_qtl_gwas_dict[tissue_name] = merge_with_gwas(date,gwas,significant_qtl_dict[tissue_name],tissue,qtl_type,trait,output_files_directory,independent_ranking)

        elif qtl_type == "ieqtl" or qtl_type == "isqtl":
            geneenrich_input_dict[tissue_name],qtl_dict[tissue_name] = parse_interaction_qtls(qtl_file,tissue,compute_tss_distance,gencode,subset_genes,null_option)

            geneenrich_input_gwas_dict[tissue_name] = merge_with_gwas(date,gwas,geneenrich_input_dict[tissue_name],tissue,qtl_type,trait,output_files_directory,independent_ranking)

            significant_qtl_dict[tissue_name]=extract_significant_interaction_qtls(qtl_dict[tissue_name],tissue,q_value,compute_tss_distance)

            significant_qtl_gwas_dict[tissue_name] = merge_with_gwas(date,gwas,significant_qtl_dict[tissue_name],tissue,qtl_type,trait,output_files_directory,independent_ranking)
        #independent
        elif qtl_type == "dapg_independent" or qtl_type == "conditional_independent":
            geneenrich_input_dict[tissue_name],qtl_dict[tissue_name] = parse_independent_eqtls(qtl_file,tissue,qtl_type,compute_tss_distance,gencode,subset_genes,null_option)
            geneenrich_input_gwas_dict[tissue_name] = merge_with_gwas(date,gwas,geneenrich_input_dict[tissue_name],tissue,qtl_type,trait,output_files_directory,independent_ranking)
            significant_qtl_dict[tissue_name]=extract_significant_independent_qtls(qtl_dict[tissue_name],tissue,independent_ranking,qtl_type,compute_tss_distance)
            significant_qtl_gwas_dict[tissue_name] = merge_with_gwas(date,gwas,significant_qtl_dict[tissue_name],tissue,qtl_type,trait,output_files_directory,independent_ranking)
        else:
            raise ValueError("qtl_type must be best_eqtl, best_sqtl, ieqtl,isqtl,dapg_independent, or conditional_independent,")

    #delete intermediate dictionaries to save memory
    del geneenrich_input_dict,qtl_dict

    return(significant_qtl_dict,significant_qtl_gwas_dict,geneenrich_input_gwas_dict)


