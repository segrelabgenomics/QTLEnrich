
"""
this script consists of the functions related to creating directories and argument parsing in QTLEnrich
"""

import sys
import os
import os.path
import re
import argparse
import logging
import datetime

def parse_args():
    """
    Arguments to pass when running eqtlEnrich.
    """
    parser = argparse.ArgumentParser(description="Run QTLEnrich on qtls with top ranked gwas P-Values.")

    #Mandatory Arguments
    parser.add_argument("-G","--gwas_file", type=str, help="gwas trait file of interest", required=True)
    parser.add_argument("-D","--qtl_directory", type=str, help="Directory containing qtl files.", required=True)
    parser.add_argument("-F","--file_name",type=str,help="Suffix for files to parse in qtl directory.",required=True)
    parser.add_argument("-Q","--qtl_type",type=str,help="Options include: best,splice,independent.",required=True)
    parser.add_argument("-C","--confounders_table", type=str, help = "Table containing variants and their corresponding confounding factors.", required=True)
    parser.add_argument("-N","--null_table",type=str,help="Table indicating which variants for a given tissues are null.",required=True)

    #Optional Arguments
    parser.add_argument("-s","--subset_genes",action="store_false",help="Indicate if user wishes to subset on protein coding and lincRNA genes.")
    parser.add_argument("-c","--compute_tss_distance",action="store_false",help="Indicate if user wishes qtlEnrich to compute distance.")
    parser.add_argument("-g","--gencode_file",type=str,help="gencode_file used for computing tss distance.")
    parser.add_argument("--keep_gene_id_suffix",action="store_true",help="Option to keep the decimal version suffix on gene ids.")
    parser.add_argument("--null_option",type=str,default = "tissue_dependent",help="Global or tissue_dependent")
    parser.add_argument("-B","--genome_build", type=str,default="b38",help = "Indicate genome build: b37 or b38.")
    parser.add_argument("--trait_name",type=str,help="Name of trait in analysis. Default takes GWAS file name.")
    parser.add_argument("-q","--qtl_q_value",default=0.05,help="qtl q-value cutoff")
    parser.add_argument("-p","--gwas_p_value",default=0.05,type=float,help="gwas p-value cutoff")
    parser.add_argument("-i","--independent_ranking",default=1,type=int,help="Independent eqtl ranking")
    parser.add_argument("-l","--lambda_factor",default=0.8,type=float,help="Lambda factor cutoff when calculating pi1")
    parser.add_argument("-L","--lower_bound_permutations",default=1000,type=int,help="Lower bound number of permutations")
    parser.add_argument("-U","--upper_bound_permutations",default=100000,type=int, help="Upper bound number of permutations")
    parser.add_argument("--num_quantiles",default=10,type=int,help="Number of quantile bins for confounding factors")
    parser.add_argument("--subset_tissues",action="store_true",help="Indicate if user wishes to subset tissues from qtl_directory.")
    parser.add_argument("--tissue_names",type=str,help="File containing tissues to keep if --subset_tissues is present.")
    parser.add_argument("--keep_null_variants",action="store_true",help="print out matched null variant IDs.")
    parser.add_argument("--keep_pvalue_matrix",action="store_true",help="print out p-values from matched null variants.")
    parser.add_argument("-O","--output_directory",type=str,help="Directory to place output files")
    parser.add_argument("-e","--exp_label",type=str,help="Unique name or detail on the analysis that will be added to output directory and output file names.")

    #GeneEnrich Input options
    parser.add_argument("-E","--GeneEnrich_input",action="store_true",help="Generate input files for eGeneEnrich")
    parser.add_argument("--eGenes_directory",type=str,help="Secondary directory to generate eGeneEnrich Input files for independent eqtls")
    parser.add_argument("--eGenes_file_name",type=str,help="File name for additionally eGenes files when preparing GeneEnrich input files.")
    return parser.parse_args()

def retrieve_current_date():
    """
    Retrieves the date.
    Date is in the form of MonthDay_Year, where month is abbreviated
    and time is based on a 24-hour clock (e.g. Mar21_2019_13:45:20).

    directory date will not include time, only MonthDay_Year.
    """
    date = datetime.datetime.now().strftime("%b%d_%Y_%H:%M:%S")
    month_date = datetime.datetime.now().strftime("%b%d_%Y")
    return (date,month_date)

def retrieve_directory():
    """
    two directories up from main scripts
    """
    return(os.path.abspath(os.path.join(os.getcwd(),"..")))

def create_output_directory_name(directory,month_date,exp_label=""):
    """
    e.g. directory: Output_Dir_QTLEnrich_exp_label_Mar21_2019
    """
    if not exp_label:
        return(directory+"/Output_QTLEnrich_"+month_date)
    return(directory+"/Output_QTLEnrich_"+exp_label+"_"+month_date)

def make_directory(directory_name):
    """
    internal function to create directory

    equivalent to mkdir in linux
    makes directory if not already exists
    """
    try:
        os.mkdir(directory_name)
    except FileExistsError:
        pass
    return(directory_name)

def create_output_files_directory(output_directory):
    """
    creates output_files directory
    """
    output_files_directory_name = str(output_directory)+"output_files"
    output_files_directory = make_directory(output_files_directory_name)
    return(output_files_directory+"/")

def create_output_directory(args,month_date,exp_label=""):
    """
    This directory should always be one directory above the main src directory.

    e.g. directory: output_directory_eqtlEnrich_v2_exp_label_Mar21_2019
    """
    if args.output_directory and args.output_directory[-1] == "/":
        return(args.output_directory)
    elif args.output_directory and args.output_directory[-1] != "/":
        return(args.output_directory+"/")
    else:
        directory=retrieve_directory()
        output_directory_name = create_output_directory_name(directory,month_date,exp_label)
        output_directory = make_directory(output_directory_name)
        return(output_directory+"/")

def create_results_directory(output_directory):
    """
    creats results directory within qtlEnrich directory in which to place qtlEnrich outputs
    """
    results_directory_name = output_directory+"results"
    #create directory
    results_directory = make_directory(results_directory_name)
    return (results_directory+"/")

def create_log_filename(output_directory,date,gwas_trait,qtl_type,exp_label=""):
    """
    Creates log filename
    """
    if not exp_label:
        return(output_directory+"QTLEnrich_"+gwas_trait+"_"+qtl_type+"_"+date+".log")
    return(output_directory+"QTLEnrich_"+exp_label+"_"+gwas_trait+"_"+qtl_type+"_"+date+".log")

def parse_exp_label(exp_label):
    """
    return empty string if none
    """
    if exp_label is None:
        return("")
    return(exp_label)

def extract_qtlenrich_arguments(args):
    """
    parses inputted arguments
    """
    #parse argument name
    argument_name = ["--"+argument for argument in vars(args)]

    #parse argument parameters
    #if none, set to empty string
    argument_parameter = [getattr(args,argument) for argument in vars(args)]
    argument_parameter = ["" if parameter is None else str(parameter) for parameter in argument_parameter]

    return(argument_name,argument_parameter)

def write_initial_logfile(output_directory,date,gwas_trait,qtl_type,arguments,parameters,exp_label=""):
    """
    write arguments and their parameters to log file
    """
    log_filename = create_log_filename(output_directory,date,gwas_trait,qtl_type,exp_label)

    with open(log_filename,"w") as f:
        f.write("QTLEnrich\n")
        f.write((date)+"\n\n")
        f.write("Running QTLEnrich on the following arugments:\n")
        for item in zip(arguments,parameters):
            f.write("\t".join(item)+"\n")
        f.write("\nOutput and log files will be placed in "+output_directory+"\n")
    return(log_filename)

def create_output_filename(results_directory,date,gwas_trait,qtl_type,independent_ranking,exp_label):
    """
    creates filename for qtlEnrich output table
    """
    if exp_label:
        if (qtl_type == "dapg_independent") or (qtl_type == "conditional_independent"):
            return(results_directory+"QTLEnrich_"+qtl_type+"_"+str(independent_ranking)+"_"+gwas_trait+"_"+exp_label+"_"+date+".tsv")
        return(results_directory+"QTLEnrich_"+qtl_type+"_"+gwas_trait+"_"+exp_label+"_"+date+".tsv")
    return(results_directory+"QTLEnrich_"+qtl_type+"_"+gwas_trait+"_"+date+".tsv")

def raise_error_gencode_file(args):
    """
    --subset_genes and --compute_tss_distance require --gencode_file
    """
    if (args.subset_genes or args.compute_tss_distance) and not args.gencode_file:
        raise argparse.ArgumentError("--gencode_file for both --compute_tss_distance and --subset_genes")
        sys.exit(1)
    else:
        pass

def subset_qtl_files(qtl_files,file_name,tissue_name):
    """
    selects tissues from tissue_name file
    """
    with open(tissue_name,"r") as f:
        tissues_keep = [x.strip() for x in f]
    #remove potential empty strings
    tissues_keep = [x for x in tissues_keep if x]

    #subset qtl_files
    return([x for x in qtl_files if re.sub(file_name,"",x) in tissues_keep])

def check_independent_GeneEnrich_parameters(args):
    """
    for independent eQTLs, needs eGenes directory and eGenes filenames
    """
    if not args.eGene_directory or not args.eGenes_file_name:
        raise argparse.ArgumentError("GeneEnrich input file parsing for independent eQTLs requires both --eGenes_directory and --eGenes_file_name")
        sys.exit(1)
    else:
        pass


