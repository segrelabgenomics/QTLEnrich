#plotting functions for qtlenrich_plotting_module

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(readr))

#### PARSING ARGUMENT FUNCTIONS ####
#qtlenrich input file
qtlenrich_file1_argument <- function(args){
if(is.na(match("--qtlenrich_file1",args))){
  stop("qtlenrich_file1 not found. Please include.")
}else {
  qtlenrich_file1_input <- match("--qtlenrich_file1",args)+1
  return(args[qtlenrich_file1_input])
  }
}

qtlenrich_file2_argument <- function(args){
if(is.na(match("--qtlenrich_file2",args))){
  print("qtlenrich_file2 not found. Assuming generating plot for one qtlenrich_file.")
  qtlenrich_file2 <- ""
  return(qtlenrich_file2)
}else {
  qtlenrich_file2_input <- match("--qtlenrich_file2",args)+1
  return(args[qtlenrich_file2_input])
  }
}

#check otions
check_qtlenrich_plotting_options <- function(args){
  if(is.na(match("--forest_plot",args)) & is.na(match("--bar_plot",args)) & is.na(match("--box_plot",args)) & is.na(match("--violin_plot",args))){
  stop("User did not indicate which plot to generate. Please add --forest_plot, --bar_plot, or --box_plot")
}
  if((!is.na(match("--paired_plot",args))) & is.na(match("--qtlenrich_file2",args))){
    stop("To generate paired_plot, two separate outputs from qtlenrich are required. It is also recommended, though not required, user adds --label_file2 for propering labeling of figures.")
  }  
}


#qtlenrich labels
label_file1_argument <- function(args){
  if(is.na(match("--label_file1",args))){
    print("label_file1 not found. Defaulting to qtlenrich_file1.")
    label_file1 <- "qtlenrich_file1"
    return(label_file1)
  }else{
    label_file1_input <- match("--label_file1",args)+1
    return(args[label_file1_input])
  }
}

label_file2_argument <- function(args){
  if(is.na(match("--label_file2",args))){
    print("label_file2 not found. Defaulting to qtlenrich_file2.")
    label_file2 <- "qtlenrich_file2"
    return(label_file2)
  }else{
  label_file2_input <- match("--label_file2",args)+1
  return(args[label_file2_input])
  }
}

#significant option
significant_option_cutoff <- function(args){
  if(is.na(match("--significant_option",args))){
    print("User did not indicate to plot only significant tissues.  If this is a mistake, please add --significant_option to command.")
    significant_option_boolean <- FALSE
    return(significant_option_boolean)
  }else{
     print("Significant option indicated. Will denote significant tissue cutoff during generation of plots.")
     significant_option_boolean <- TRUE
     return(significant_option_boolean)
  }
}

identify_significant_cutoff <- function(args,significant_option_boolean){
  if(isFALSE(significant_option_boolean) & is.na(match("--significant_tissues_cutoff",args))){
    significant_cutoff <- NA
    return(significant_cutoff)
  }else if(isTRUE(significant_option_boolean) & is.na(match("--significant_tissues_cutoff",args))){
    print("User wishes to plot only significant tissues. However, cutoff was not given. Defaulting to bonferroni cutoff.")
    significant_cutoff <- "bonferroni"
    return(significant_cutoff)
   }else{
      significant_cutoff_input <- match("--significant_tissues_cutoff",args)+1
      significant_cutoff <- args[significant_cutoff_input]
      significant_cutoff <- as.double(significant_cutoff)
      return(significant_cutoff)
    }
}

identify_number_significant_tissues <- function(args){
  if(!is.na(match("--significant_option",args)) &  is.na(match("--subset_significant_tissues",args))){
    print("User did not select the top number of significant tissues. Will plot all significant tissues.")
    number_significant_tissues <- NA
    return(number_significant_tissues)
  }else if(!is.na(match("--subset_significant_tissues",args))){
     number_significant_tissues_input <- match("--subset_significant_tissues",args)+1
     return(as.numeric(args[number_significant_tissues_input]))
   }
}

#remove tissues option
remove_tissues <- function(args){
  if(is.na(match("--remove_subset_tissues",args))){
    tissues_list_to_remove <- NA
  }else{
    tissues_list_file <- args[match("--remove_subset_tissues",args)+1]
    tissues_list_to_remove <- suppressMessages(read_csv(tissues_list_file,col_names=FALSE)) %>%
      as.data.frame()
    colnames(tissues_list_to_remove) <- "Tissue_(QTL)"
    return(tissues_list_to_remove) 
   }
}

#box, violin, bar plot option
retrieve_bar_violin_bar_option <- function(args){
  if(is.na(match("--bar_box_violin_option",args))){
    bar_box_violin_option <- "pi1_trait"
  }
  else{
    bar_box_violin_option <- args[match("--bar_box_violin_option",args)+1]
  }
}

retrieve_bar_violin_bar_header <- function(bar_box_violin_option){
  if (bar_box_violin_option == "pi1_trait"){
    return("Estimated_total_num_trait_associations")
  }else if (bar_box_violin_option == "EmpPi1"){
   return("Empirical_Pi1")
  }else if (bar_box_violin_option == "adjusted_fold"){
   return("Adjusted_Fold_Enrichment")
   }else {
     stop("incorrect option for bar, violin, or box plot was inputted. Options include true_trait, EmpPi1, or adjusted_fold.")
   }
}

retrieve_bar_box_violin_axes_label <- function(bar_box_violin_header){
  if (bar_box_violin_header == "Estimated_total_num_trait_associations"){
     return("Estimated # trait associations based on EmpPi1")
  } else if (bar_box_violin_header == "Empirical_Pi1"){
    return("Empirical Pi1")
  }else{
     return("Adjusted Fold-Enrichment")
   }
}
