#!/usr/bin/env python3

"""
Python script for generating qqplots from QTLEnrich output files

Author: Andrew Hamel
Date: 24 January 2020

Written in the Segre lab at Massachusetts Eye and Ear, 
Department of Ophthalmology, Harvard Medical School, Boston, MA

Required Python3, scipy, numpy, pandas, and matplotlib 

Inputs:

Required inputs:

1. significant qtl file (outputted from QTLEnrich)
  - tab-delimited file comprised of two columns: variant and gwas_p_value
2. either GWAS file or pvalue matrix (output from QTLEnrich)
s module requires at least one of the below inputs  - GWAS file is tab-delimited and must include a gwas_p_value column
3. gwas_background or matched_null_background
  - Option to plot gwas background or matched null background. 
  - include both if want both plots (two separate plots)

Optional inputs

1. num_permutations
  - indicates number of permutations from null pvalue matrix to plot (default is first 100)
2. confidence_interval
  - option to shade confidence interval of null pvalues
3. upper_bound_ci
  - Select upper bound confidence interval (default is 97.5 percentile)
4. lower_bound_ci
  - Select lower bound confidence interval (default is 2.5 percentile)
5. upper_bound_clip
  - set upper bound on qq-plot 
6. trait_label
  - label for gwas trait
7. qtl_type
  - type of QTLs in the plot
8. tissue_label
  - type of tissue in this plot
9. show_title
  - option to include title in plot
10. qqplot_title
  - option to input own title
11. output_directory
  - place plots in directory of choosing. Defaults to creating one automatically
12. exp_label
  - add label to filenames
13. xinches
  - number of horizontal inches on plot
14. y_inches
  - number of vertical inches
15. dpi
  - resolution

"""

import matplotlib
matplotlib.use('Agg')
import os
import os.path
import datetime
import re
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def create_gwas_background_filename(output_directory,date,exp_label,qtl_type,tissue_label,trait_label):
    """
    filename for qqplot with gwas background
    """
    if exp_label:
        return(output_directory+"QTLEnrich_qqplot_"+trait_label+"_"+tissue_label+"_"+qtl_type+"_"+exp_label+"_gwas_backgr_"+date)
    return(output_directory+"QTLEnrich_qqplot_"+trait_label+"_"+tissue_label+"_"+qtl_type+"_gwas_backgr_"+date+)

def create_null_background_filename(output_directory,date,exp_label,qtl_type,tissue_label,trait_label):
    """
    filename for qqplot with null background
    """
    if exp_label:
        return(output_directory+"QTLEnrich_qqplot_"+trait_label+"_"+tissue_label+"_"+qtl_type+"_"+exp_label+"_null_backgr_"+date)
    return(output_directory+"QTLEnrich_qqplot_"+trait_label+"_"+tissue_label+"_"+qtl_type+"_null_backgr_"+date)

def set_clip(df,upper_bound_clip):
    """
    sets upper bound clip
    """
    if upper_bound_clip == 0:
        return(df)
    else:
        df = df.clip(upper = upper_bound_clip)
        return(df)    

def add_expected_pvalue_headers(df):
    """
    Add expected pvalue column and multiply p-values by -log10
    """
    df = df.sort_values("gwas_p_value",ascending=True)
    df["exp"] = np.linspace(start = 1/len(df),stop = 1,num=len(df),endpoint=True)
    df["obs_log10"] = -1*np.log10(df["gwas_p_value"])
    df["exp_log10"] = -1*np.log10(df["exp"])
    df = df[["obs_log10","exp_log10"]]
    return(df)

def format_title(trait_label,title,show_title=False):
    """
    return title for qqplot
    """
    if show_title:
        if not title:
            return(" ".join(trait_label))
        return(" ".join(title))
    return("")

def format_signifcant_qtl_label(qtl_type,tissue_label):
    return(" ".join(tissue_label)+" "+" ".join(qtl_type))

def prepare_gwas(gwas_file):
    """
    Clean gwas Dataset to account for differences in
    notation of chromosomes and select relevant columns.
    Position will be treated as integers
    """
    gwas = pd.read_csv(gwas_file,sep='\t')
    gwas['gwas_p_value'] = gwas['gwas_p_value'].replace(0.0,1e-280)
    return(gwas)

def parse_qtl(qtl_file):
    """
    Filters on qvalue if in header.
    Merges with GWAS
    """
    qtl_set = pd.read_csv(qtl_file,sep='\t')
    return(qtl_set)

def parse_pvalue_null_matrix(pvalue_null_matrix_file,num_permutations=100):
    """
    selects number of permutations
    """
    pvalue_null_matrix = pd.read_csv(pvalue_null_matrix_file,sep=",")
    pvalue_null_matrix = pvalue_null_matrix.rename(columns={"Unnamed: 0":"variant_id"}).drop("variant_id",axis=1)
    pvalue_null_matrix = pvalue_null_matrix.iloc[:,0:num_permutations]
    return(pvalue_null_matrix)

def retrieve_current_date():
    """
    Retrieves the date.
    Date is in the form of MonthDay_Year, where month is abbreviated
    MonthDay_Year, e.g. Mar24_2020
    """
    date = datetime.datetime.now().strftime('%b%d_%Y')
    return(date)

def retrieve_directory():
    """
    two directories up from main scripts
    """
    return(os.path.abspath(os.path.join(os.getcwd(),'..')))

def create_output_directory_name(directory,date,exp_label=""):
    """
    e.g. directory: Output_QTLEnrich_exp_label_Mar21_2019
    """
    if not exp_label:
        return(directory+'/Output_QTLEnrich_'+date)
    return(directory+'/Output_QTLEnrich_'+exp_label+'_'+date)

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

def make_plotting_directory(output_directory_name):
    """
    within output directory, sub directory entitled "plots"
    """
    if output_directory_name[-1] != "/":
        plot_directory_name = output_directory_name+"/"+"plots"
    else:
        plot_directory_name = output_directory_name+"plots"
    plot_directory = make_directory(plot_directory_name)
    return(plot_directory+"/")

def create_output_directory(args,date,exp_label=""):
    """
    This directory should always be one directory above the main src directory.

    e.g. directory: output_directory_exp_label_Mar21_2019
    """
    if args.output_directory and args.output_directory[-1] == '/':
        return(args.output_directory)
    elif args.output_directory and args.output_directory[-1] != '/':
        return(args.output_directory+'/')
    else:
        directory=retrieve_directory()
        output_directory_name = create_output_directory_name(directory,date,exp_label)
        output_directory = make_directory(output_directory_name)
        plot_directory = make_plotting_directory(output_directory_name)
        return(plot_directory+'/')

def parse_exp_label(exp_label):
    """
    return empty string if none
    """
    if exp_label is None:
        return('')
    return(exp_label)

def raise_error_qq_plot_options(args):
    """
    requires either --gwas_background or --matched_null_background
    """
    if (not args.gwas_background) and (not args.matched_null_background):
        raise argparse.ArgumentError("Did not include which type of qqplot to generate. Requires either --gwas_background or --matched_null_background.")
        sys.exit(1)
    else:
        pass

def abline(m, b):
    """
    Plot a line from slope and intercept
    m: slope,b: intercept
    """
    axes = plt.gca()
    x = np.array(axes.get_xlim())
    x[0] = 0
    y = b + m * x
    plt.plot(x, y,linewidth=3,color="darkgray",zorder=0)

def generate_qqplot_gwas_background(gwas_background_filename,significant_qtl,gwas,qqplot_title,signifcant_qtl_label,upper_bound_clip = 0, x_inches=12, y_inches=10, dpi=100):
    """
    Produces qqplot consisting of significant QTLs and full GWAS
    """

    #add expected pvalues
    significant_qtl_expected = add_expected_pvalue_headers(significant_qtl)
    significant_qtl_expected = set_clip(significant_qtl_expected,upper_bound_clip)

    gwas_expected = add_expected_pvalue_headers(gwas,upper_bound_clip)
    gwas_expected = set_clip(gwas_expected,upper_bound_clip)

    print("plotting")

    #for dynamic resolution scaling
    x_pixels = int(x_inches*dpi)
    y_pixels = int(y_inches*dpi)

    fig, ax = plt.subplots(figsize=(x_inches, y_inches),dpi=dpi)

    x_gwas = list(gwas_expected["exp_log10"])
    y_gwas = list(gwas_expected["obs_log10"])

    x_qtl = list(significant_qtl_expected["exp_log10"])
    y_qtl = list(significant_qtl_expected["obs_log10"])

    #For redundant point removal
    gwas_expected['x_gwas_scaled'], gwas_expected['y_gwas_scaled'] = pixel_rounder(x_pixels, y_pixels, x_gwas, y_gwas)
    gwas_expected = gwas_expected[['x_gwas_scaled','y_gwas_scaled']]
    gwas_expected = gwas_expected.drop_duplicates()
    
    significant_qtl_expected["x_qtl_scaled"], significant_qtl_expected["y_qtl_scaled"] = pixel_rounder(x_pixels, y_pixels, x_qtl, y_qtl)
    significant_qtl_expected = significant_qtl_expected[["x_qtl_scaled","y_qtl_scaled"]]
    significant_qtl_expected = significant_qtl_expected.drop_duplicates()

    ax.scatter(gwas_expected['x_gwas_scaled'], gwas_expected['y_gwas_scaled'],c="k",label="Full GWAS",zorder=5)
    abline(m=1,b=0)
    ax.scatter(significant_qtl_expected["x_qtl_scaled"], significant_qtl_expected["y_qtl_scaled"],c="red",label=signifcant_qtl_label,zorder=10)
    ax.set_xlabel(r'EXP $\mathregular{log_{10}}$(GWAS P-value)',fontsize=30)
    ax.set_ylabel(r'OBS $\mathregular{log_{10}}$(GWAS P-value)',fontsize=30)
    ax.tick_params(labelsize=20)

    #title
    ax.axes.set_title(qqplot_title,fontsize=30)
    #legend
    ax.legend(handletextpad=0.1,prop={"size":20},loc=2,markerscale=2)

    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    print("saving")
    fig.savefig(gwas_background_filename+".png")
    fig.savefig(gwas_background_filename+".pdf")
    plt.close(fig)    

def sort_pvalue_null_matrix(pvalue_null_matrix):
    """
    sorts each column of pvalue matrix independently
    """
    pvalues = pvalue_null_matrix.values
    pvalues.sort(axis=0)
    pvalues = pvalues[::-1]
    new_pvalue_null_matrix_dataframe = pd.DataFrame(pvalues,index=pvalue_null_matrix.index,columns=pvalue_null_matrix.columns).T
    return(new_pvalue_null_matrix_dataframe)

def extract_percentiles(pvalue_null_matrix_sorted,upper_percentile=97.5,lower_percentile=2.5):
    """
    computes upper and lower bound percentiles for each rank
    """
    dict_lower = {}
    dict_upper = {}
    for i in pvalue_null_matrix_sorted:
        rank = pvalue_null_matrix_sorted[[i]]
        dict_lower[i] = np.percentile(rank[i],q=lower_percentile)
        dict_upper[i] = np.percentile(rank[i],q=upper_percentile)

    pvalue_confidence_interval_lower = pd.DataFrame(dict_lower.items(),columns = ["permut","gwas_p_value"]).drop("permut",axis=1)
    pvalue_confidence_interval_upper = pd.DataFrame(dict_upper.items(),columns = ["permut","gwas_p_value"]).drop("permut",axis=1)

    return(pvalue_confidence_interval_lower,pvalue_confidence_interval_upper)

def extract_confidence_interval(upper,lower,pvalue_null_matrix):
    """
    returns dataframe of upper and lower confidence intervals for each rank
    """
    #sort pvalue matrix columns separately
    pvalue_null_matrix_sorted = sort_pvalue_null_matrix(pvalue_null_matrix)

    #obtain lower and upper percentiles
    pvalue_confidence_interval_lower,pvalue_confidence_interval_upper = extract_percentiles(pvalue_null_matrix_sorted,upper,lower)

    return(pvalue_confidence_interval_lower,pvalue_confidence_interval_upper)

def interpolate_shading(df):
    """
    Creates linear function based on confidence intervals
    """
    x_axis_values = np.array(df["exp_log10"])
    y_axis_values = np.array(df["obs_log10"])

    x_new = np.linspace(x_axis_values.min(),x_axis_values.max(),500) 
    f = interp1d(x_axis_values, y_axis_values, kind='linear')
    y_smooth=f(x_new)
    df_new = pd.DataFrame({"x":x_new,"y":y_smooth})
    return(df_new)

def round_to(pixels, my_range):
    if abs(my_range) > 1:
        round_to = len(str(round(abs(pixels)))) - len(str(round(abs(my_range))))+1
    else:
        leading_zeroes = abs(int(np.log10(abs(my_range))))
        round_to = leading_zeroes +len(str(round(abs(pixels))))+1
    return(round_to)

def pixel_rounder(x_pixels, y_pixels, x_values, y_values):
    """
    x_pixels: number of pixels on x-axis
    y_pixels: number of pixels on y-axis
    x_values: values on x-axis
    y_values: values on y-axis
    """
    x_values = list(x_values)
    y_values = list(y_values)

    #x and y ranges
    x_range = max(x_values) - min(x_values)
    y_range = max(y_values) - min(y_values)

    #figure out how many to round to for x & y separately
    x_to = round_to(x_pixels, x_range)
    y_to = round_to(y_pixels, y_range)
    #do the rounding
    new_x_values = np.round(x_values, x_to)
    new_y_values = np.round(y_values, y_to)
    return(new_x_values, new_y_values)

def plot_sampled_null(ax,pvalue_null_matrix,upper_bound_clip,x_pixels,y_pixels):
    """
    plots sampled null pvalues to qqplot
    """
    #sampled null
    pvalue_null_matrix_separated = pvalue_null_matrix.iloc[:,0]
    pvalue_null_matrix_separated = pd.DataFrame(pvalue_null_matrix_separated)
    pvalue_null_matrix_separated = pvalue_null_matrix_separated.rename(columns={"P0":"gwas_p_value"})
    pvalue_null_matrix_separated_expected = add_expected_pvalue_headers(pvalue_null_matrix_separated)
    ax.scatter(pvalue_null_matrix_separated_expected["exp_log10"],pvalue_null_matrix_separated_expected["obs_log10"],label="Matched null variants",c="lightgray")

    full_df = pd.DataFrame()

    for i in pvalue_null_matrix.iloc[:,1:]:
        pvalue_null_matrix_separated = pvalue_null_matrix[[i]]
        pvalue_null_matrix_separated = pvalue_null_matrix_separated.rename(columns={i:"gwas_p_value"})
        pvalue_null_matrix_separated_expected = add_expected_pvalue_headers(pvalue_null_matrix_separated)
        pvalue_null_matrix_separated_expected = set_clip(pvalue_null_matrix_separated_expected,upper_bound_clip)

        full_df = full_df.append(pvalue_null_matrix_separated_expected)

    #For redundant point removal
    x = full_df["exp_log10"]
    y = full_df["obs_log10"]
    full_df['x_scaled'],full_df['y_scaled'] = pixel_rounder(x_pixels, y_pixels, x, y)
    full_df = full_df[["x_scaled","y_scaled"]]
    full_df = full_df.drop_duplicates()

    ax.scatter(full_df['x_scaled'],full_df['y_scaled'],c="lightgray")

def plot_confidence_intervals_with_shading(ax,pvalue_null_matrix,upper,lower):
    """
    adds confidence intervals to qqplot and shades in regions
    """
    pvalue_confidence_interval_lower,pvalue_confidence_interval_upper = extract_confidence_interval(upper,lower,pvalue_null_matrix)
    pvalue_confidence_interval_lower_expected = add_expected_pvalue_headers(pvalue_confidence_interval_lower)
    pvalue_confidence_interval_upper_expected = add_expected_pvalue_headers(pvalue_confidence_interval_upper)

    pvalue_shade_lower = interpolate_shading(pvalue_confidence_interval_lower_expected)
    pvalue_shade_upper = interpolate_shading(pvalue_confidence_interval_upper_expected)

    plt.fill_between(pvalue_shade_lower["x"],pvalue_shade_lower["y"],pvalue_shade_upper["y"],color="lightgray",label="95% CI",alpha=0.5,zorder=10)

def generate_qqplot_null_background(null_background_filename,significant_qtl,pvalue_null_matrix,qqplot_title,significant_qtl_label,upper,lower,confidence_interval=False,upper_bound_clip=0,x_inches=12, y_inches=10, dpi=100):

    #for dynamic resolution scaling
    x_pixels = int(x_inches*dpi)
    y_pixels = int(y_inches*dpi)

    #add expected pvalues
    significant_qtl_expected = add_expected_pvalue_headers(significant_qtl)
    significant_qtl_expected = set_clip(significant_qtl_expected,upper_bound_clip)

    x_qtl = list(significant_qtl_expected["exp_log10"])
    y_qtl = list(significant_qtl_expected["obs_log10"])

    #For redundant point removal
    significant_qtl_expected["x_qtl_scaled"], significant_qtl_expected["y_qtl_scaled"] = pixel_rounder(x_pixels, y_pixels, x_qtl, y_qtl)
    significant_qtl_expected = significant_qtl_expected[["x_qtl_scaled","y_qtl_scaled"]]
    significant_qtl_expected = significant_qtl_expected.drop_duplicates()

    #plotting
    fig, ax = plt.subplots(figsize=(x_inches, y_inches),dpi=dpi)

    #confidence interval if needed
    if confidence_interval:
        plot_confidence_intervals_with_shading(ax,pvalue_null_matrix,upper,lower)

    else:
    #if no confidence interval, plot matched null points
        plot_sampled_null(ax,pvalue_null_matrix,x_pixels,y_pixels)

    #diagonal
    abline(m=1,b=0)

    #title
    ax.axes.set_title(qqplot_title,fontsize=30)

    #significan qtl
    ax.scatter(significant_qtl_expected["x_qtl_scaled"],significant_qtl_expected["y_qtl_scaled"],c="red",label=significant_qtl_label,zorder=20)
    ax.set_xlabel(r'EXP $\mathregular{log_{10}}$(GWAS P-value)',fontsize=30)
    ax.set_ylabel(r'OBS $\mathregular{log_{10}}$(GWAS P-value)',fontsize=30)
    ax.tick_params(labelsize=20)

    #legend
    handles,labels = ax.get_legend_handles_labels()

    if confidence_interval:
        handles = [handles[-1],handles[0]]
        labels = [labels[-1],labels[0]]

    if not confidence_interval:
        handles = [handles[-1],handles[0]]
        labels = [labels[-1],labels[0]]

    ax.legend(handles,labels,handletextpad=0.1,prop={"size":20},loc=2,markerscale=2)

    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    print("saving")
    fig.savefig(null_background_filename+".png")
    fig.savefig(null_background_filename+".pdf")
    plt.close(fig)
 
def parse_label(label):
    """
    Joins spaces
    """
    label = " ".join(label)
    return(re.sub(" ","",label))

def parse_args():
    """
    Arguments to pass to generate qq-plot.
    """
    parser = argparse.ArgumentParser(description='Generate qq-plot from QTLEnrich output.')

    parser.add_argument('-Q','--qtl_file',type=str,help="qtl file containing significant variants",required=True)
    parser.add_argument('-G','--gwas_file', type=str, help='gwas trait file of interest')
    parser.add_argument('-M','--pvalue_null_matrix',type=str,help="pvalue matrix output from qtlenrich")
    parser.add_argument('--gwas_background',action="store_true",help='Plot qq-plot of significantl QTLs against gwas background.')
    parser.add_argument('--matched_null_background',action="store_true",help="Plot qq-plot of significan QTLs against matched null variants.")
    parser.add_argument('--num_permutations',type=int,default=100,help="Number of permutations user wishes to plot. QTLEnrich outputs first 1000. Default is 100.")
    parser.add_argument('--confidence_interval',action="store_true",help="Add confidence intervals to --matched_null_background.")
    parser.add_argument('--upper_bound_ci',type=float,default=97.5,help="Upper bound confidence interval. Default: 97.5 percentile")
    parser.add_argument('--lower_bound_ci',type=float,default=2.5,help="Lower bound confidence interval. Default: 2.5 percentile")
    parser.add_argument('--upper_bound_clip',type=int,default=0,help='set an upper bound on qq-plot. Default is 0, indicating no upper bound.')
    parser.add_argument('--trait_label',type=str,nargs="+",default="Trait",help="Denotes name of trait for GWAS study, e.g. 'Standing Height'. Default is 'Trait'.")
    parser.add_argument('--qtl_type',type=str,nargs="+",default="Significant QTLs",help="Denotes name of significant qtls, e.g. 'eQTLs'. Default is 'Significant QTLs'.")
    parser.add_argument('--tissue_label',type=str,nargs="+",default="Tissue",help="Denotes name of tissue, e.g. 'Whole Blood'. Default is 'Tissue'.")
    parser.add_argument('--show_title',action="store_true",help="User indicates the addition of title to qqplot.")
    parser.add_argument('--qqplot_title',type=str,nargs="+",help="If user chooses own title.")
    parser.add_argument('-O','--output_directory',type=str,help='Directory to place output files')
    parser.add_argument('-e','--exp_label',type=str,help='Unique name or detail on the analysis that will be added to output directory and output file names.')

    #For dynamic resolution scaling
    parser.add_argument('--x_inches', type=float, default=(12.),help='Number of horizontal inches on the plot.')
    parser.add_argument('--y_inches', type=float, default=(10.),help='Number of vertical inches on the plot.')
    parser.add_argument('--dpi', type=float, default=100.,help='Display pixels per inch.')
    return parser.parse_args()

if __name__ == "__main__":

    args = parse_args()
    #check specific arguments
    raise_error_qq_plot_options(args)
    #retrieve dates
    date = retrieve_current_date()
    #parse exp label
    exp_label = parse_exp_label(args.exp_label)
    #parse other labels
    trait_label = parse_label(args.trait_label)
    tissue_label = parse_label(args.tissue_label)
    qtl_type = parse_label(args.qtl_type)
    #Create output directory
    output_directory = create_output_directory(args,date,exp_label)
    #format qqplot title
    qqplot_title = format_title(args.trait_label,args.qqplot_title,args.show_title)
    #format significant qtl label
    significant_qtl_label = format_signifcant_qtl_label(args.qtl_type,args.tissue_label)
    #parse significant qtl and gwas files
    print("parsing significant qtl")
    significant_qtl = parse_qtl(args.qtl_file)

    #qqplot with significant qtls and full gwas
    if args.gwas_background:
        print("Preparing to plot qq-plot with gwas background.")
        print("parsing gwas")
        gwas = prepare_gwas(args.gwas_file)

        #create filename
        gwas_background_filename = create_gwas_background_filename(output_directory,date,exp_label,qtl_type,tissue_label,trait_label)

        generate_qqplot_gwas_background(gwas_background_filename,significant_qtl,gwas,qqplot_title,significant_qtl_label,args.upper_bound_clip, args.x_inches, args.y_inches, args.dpi)

    #qqplot with significant qtls and matched null gwas p-values
    if args.matched_null_background:
        print("Preparing to plot qq-plot with matched null background")
        print("parsing p-value matrix")
        pvalue_null_matrix_parsed = parse_pvalue_null_matrix(args.pvalue_null_matrix,args.num_permutations)

        #create filename
        null_background_filename = create_null_background_filename(output_directory,date,exp_label,qtl_type,tissue_label,trait_label)

        generate_qqplot_null_background(null_background_filename,significant_qtl,pvalue_null_matrix_parsed,qqplot_title,significant_qtl_label,args.upper_bound_ci,args.lower_bound_ci,args.confidence_interval,args.upper_bound_clip, args.x_inches, args.y_inches, args.dpi)
