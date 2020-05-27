#!/usr/bin/env Rscript

# To generate plots to analyze QTLEnrich outputs

print("Loading libraries.  Requires ggplot2, dplyr, readr, and cowplot.")

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(cowplot))

print("Libraries were loaded effectively.")

source("plotting_forest_boxplot_functions.R")
source("parsing_forest_boxplot_functions.R")
source("os_sys_forest_boxplot_functions.R")

args <- commandArgs(trailingOnly = TRUE)

#qtlenrich input files
qtlenrich_file1 <- qtlenrich_file1_argument(args)
qtlenrich_file2 <- qtlenrich_file2_argument(args)

check_qtlenrich_plotting_options(args)
print("Required options were found. Proceeding through argument parsing.")

#label for qtlenrich output files
label_file1 <- label_file1_argument(args)
label_file2 <- label_file2_argument(args)

#significant option
significant_option_boolean <- significant_option_cutoff(args) 
significant_cutoff <- identify_significant_cutoff(args,significant_option_boolean)
number_significant_tissues <- identify_number_significant_tissues(args)

#bar,box,violin plot option
bar_box_violin_option <- retrieve_bar_violin_bar_option(args)
bar_box_violin_header <- retrieve_bar_violin_bar_header(bar_box_violin_option) 
bar_box_violin_axes_label <- retrieve_bar_box_violin_axes_label(bar_box_violin_header)

#remove tissues list
tissues_list_to_remove <- remove_tissues(args)

#alternate tissues file to remove certain tissues from plots
gtex_tissue_names <- parse_alternate_gtex_tissues(args)

#output directory
date <- retrieve_date()
exp_label <- parse_exp_label(args)
output_directory <- create_output_directory(args,exp_label,date)

#generating plots
print("Arguments parsed. Now generating plots.")
print("Parsing QTLEnrich files.")

#qtlenrich file1
qtlenrich_file1_parsed <- parse_qtlenrich_file(qtlenrich_file1,label_file1,gtex_tissue_names,significant_cutoff,number_significant_tissues,tissues_list_to_remove)

#qtlenrich file2
if(!is.na(match("--paired_plot",args))){
  qtlenrich_file2_parsed <- parse_qtlenrich_file2(qtlenrich_file1_parsed,qtlenrich_file2,label_file2,gtex_tissue_names)
  qtlenrich_files <- rbind(qtlenrich_file1_parsed,qtlenrich_file2_parsed)
}

#aes strings for special cases
qtlenrich_file1_parsed$adjusted_sorted_2 <- reorder(qtlenrich_file1_parsed$tissue_site_detail,qtlenrich_file1_parsed$`Adjusted_Fold_Enrichment`)
adjusted_sorted_2 <- "adjusted_sorted_2"
type <- "type"

#dynamic sizes for multiple tissues
if (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 5){
  plot_width <- 8
  plot_height <- 6
}

#forest plot

if(!is.na(match("--forest_plot",args)) & is.na(match("--paired_plot",args))){
  print("Generating forest plot for single trait.")

if (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 10){
  plot_width <- 6
  plot_height <- 4
  } else {
    plot_width <- 14
    plot_height <- 10
  }

  forest_plot <- generate_single_forest_plot(qtlenrich_file1_parsed)
  forest_plot_filename <- create_single_forest_plot_filename(output_directory,label_file1,exp_label,date)
  ggsave(forest_plot_filename,plot=forest_plot,device="pdf",width=plot_width,height=plot_height,units="in") #14.10

}

if( !is.na(match("--forest_plot",args)) & !is.na(match("--paired_plot",args))){
  print("Generating forest plot for two qtlenrich output files..")


  if (length(unique(qtlenrich_files$`Tissue_(QTL)`)) <= 5){
      plot_width <- 10
      plot_height <- 6
   }else if ( (length(unique(qtlenrich_files$`Tissue_(QTL)`)) > 5) & (length(unique(qtlenrich_files$`Tissue_(QTL)`)) <= 20)){
      plot_width <- 14
      plot_height <- 10
     } else if ((length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) > 20) & (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 35)) {
      plot_width <- 14
      plot_height <- 10
     } else {
      plot_width <- 20
      plot_height <- 45
    }


  forest_plot <- generate_paired_forest_plot(qtlenrich_files,label_file1,label_file2)
  paired_forest_plot_filename <- create_paired_forest_plot_filename(output_directory,label_file1,label_file2,exp_label,date)
  ggsave(paired_forest_plot_filename,plot=forest_plot,device="pdf",width=plot_width,height=plot_height,units="in")

}
#bar plot
if(!is.na(match("--bar_plot",args)) & is.na(match("--paired_plot",args))){
  print("Generating bar plot for single trait.")

if (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 10){
  plot_width <- 6
  plot_height <- 4
  } else {
    plot_width <- 14
    plot_height <- 10
  }

  bar_plot <- generate_single_bar_plot(qtlenrich_file1_parsed,adjusted_sorted_2,bar_box_violin_header,bar_box_violin_axes_label)
  bar_plot_filename <- create_single_bar_plot_filename(output_directory,label_file1,exp_label,date)
  ggsave(bar_plot_filename,plot=bar_plot,device="pdf",width=plot_width,height=plot_height,units="in")

}

if((!is.na(match("--bar_plot",args))) & !is.na(match("--paired_plot",args))){
  print("Generating bar plot for two qtlenrich output files..")

  if (length(unique(qtlenrich_files$`Tissue_(QTL)`)) <= 2){
    plot_width <- 14 #16,12
    plot_height <- 10
    } else {
      plot_width <- 20
      plot_height <- 15
    }

  bar_plot <- generate_paired_bar_plot(qtlenrich_files,label_file1,label_file2,type,bar_box_violin_header,bar_box_violin_axes_label)
  paired_bar_plot_filename <- create_paired_bar_plot_filename(output_directory,label_file1,label_file2,exp_label,date)
  ggsave(paired_bar_plot_filename,plot=bar_plot,device="pdf",width=plot_width,height=plot_height,units="in")

}

#side by side option not paired
if(!is.na(match("--bar_plot",args)) & !is.na(match("--forest_plot",args)) & is.na(match("--paired_plot",args)) & !is.na(match("--side_by_side",args))){

  print("Generating forest and bar plots side by side.")

if (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 5){
  plot_width <- 10
  plot_height <- 8
  } else if ((length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) > 5) & (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 10)){
  plot_width <- 14
  plot_height <- 10
  }
    else {
    plot_width <- 30
    plot_height <- 20
  }

  side_by_side_plot <- generate_unpaired_plot_side_by_side(qtlenrich_file1_parsed,adjusted_sorted_2,bar_box_violin_header,bar_box_violin_axes_label)
  side_by_side_filename <- create_unpaired_side_by_side_plot_filename(output_directory,label_file1,exp_label,date)
  ggsave(side_by_side_filename,plot=side_by_side_plot,device="pdf",width=plot_width,height=plot_height,units="in") 

}

#side by side option paired
if(!is.na(match("--bar_plot",args)) & !is.na(match("--forest_plot",args)) & !is.na(match("--paired_plot",args)) & !is.na(match("--side_by_side",args))){

   print("Generating paired forest and bar plots side by side.")

  if (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 5){
    plot_width <- 14
    plot_height <- 6
    } else if ((length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) > 5) & (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 20)){
    plot_width <- 16
    plot_height <- 20
    } else {
    plot_width <- 30
    plot_height <- 40
    }

  side_by_side_plot <- generate_paired_plot_side_by_side(qtlenrich_files,label_file1,label_file2,type,bar_box_violin_header,bar_box_violin_axes_label)
  side_by_side_filename <- create_paired_side_by_side_plot_filename(output_directory,label_file1,label_file2,exp_label,date)
  ggsave(side_by_side_filename,plot=side_by_side_plot,device="pdf",width=plot_width,height=plot_height,units="in")

}

#box plot
if(!is.na(match("--box_plot",args)) & is.na(match("--paired_plot",args))){
  print("Generating box plot for single trait.")
  qtlenrich_file1_parsed$type <- label_file1

  box_plot <- generate_single_box_plot(qtlenrich_file1_parsed,label_file1,type,bar_box_violin_header,bar_box_violin_axes_label)
  box_plot_filename <- create_single_box_plot_filename(output_directory,label_file1,exp_label,date)

  ggsave(box_plot_filename,plot=box_plot,device="pdf",width=5,height=4,units="in")

}

if((!is.na(match("--box_plot",args))) & !is.na(match("--paired_plot",args))){
  print("Generating bar plot for two qtlenrich output files..")

  box_plot <- generate_paired_box_plot(qtlenrich_files,label_file1,label_file2,type,bar_box_violin_header,bar_box_violin_axes_label)
  paired_box_plot_filename <- create_paired_box_plot_filename(output_directory,label_file1,label_file2,exp_label,date)
  ggsave(paired_box_plot_filename,plot=box_plot,device="pdf",width=6,height=5,units="in")

}

#violin plot
if(!is.na(match("--violin_plot",args)) & is.na(match("--paired_plot",args))){
  print("Generating violin plot for single trait.")

  qtlenrich_file1_parsed$type <- label_file1
  violin_plot <- generate_single_violin_plot(qtlenrich_file1_parsed,label_file1,type,bar_box_violin_header,bar_box_violin_axes_label)
  violin_plot_filename <- create_single_violin_plot_filename(output_directory,label_file1,exp_label,date)

  ggsave(violin_plot_filename,plot=violin_plot,device="pdf",width=4,height=3,units="in")

}

if((!is.na(match("--violin_plot",args))) & !is.na(match("--paired_plot",args))){
  print("Generating violin plot for two qtlenrich output files..")

  violin_plot <- generate_paired_violin_plot(qtlenrich_files,label_file1,label_file2,type,bar_box_violin_header,bar_box_violin_axes_label)
  paired_violin_plot_filename <- create_paired_violin_plot_filename(output_directory,label_file1,label_file2,exp_label,date)
  ggsave(paired_violin_plot_filename,plot=violin_plot,device="pdf",width=6,height=5,units="in")

}

paste("Plots are finished. Placing plots in output directory:",output_directory,sep=" ")
