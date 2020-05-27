
#plotting functions for qtlenrich_plotting_module

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(readr))

#### OS,SYS FUNCTIONS ####
retrieve_date <- function(){
  return(format(Sys.Date(),"%b%d_%Y"))
}

#output directory
parse_exp_label <- function(args){
  if(is.na(match("--exp_label",args))){
    exp_label <- ""
    return(exp_label)
  }else{
     exp_label_input <- match("--exp_label",args)+1
     return(args[exp_label_input])
   }
}

create_output_directory_filename <- function(exp_label,current_date){
  if (exp_label == ""){
  return(paste("../Output_QTLEnrich",current_date,sep="_"))}
  else {
  return(paste("../Output_QTLEnrich",exp_label,current_date,sep="_"))
}
}

format_output_directory <- function(output_directory){
if(substr(output_directory,nchar(output_directory),nchar(output_directory)) != "/"){
  return(paste(output_directory,"/",sep=""))
}else{
  return(output_directory)
  }
}

create_output_directory <- function(args,exp_label,current_date){
  if(is.na(match("--output_directory",args))){
    output_directory <- create_output_directory_filename(exp_label,current_date)
    if(!dir.exists(output_directory)){
    dir.create(output_directory)}
    plot_directory <- paste(output_directory,"/plots",sep="")
    if(!dir.exists(plot_directory)){ 
    dir.create(plot_directory)}
    plot_directory <- format_output_directory(plot_directory)
    return(plot_directory)
}else{
  output_directory_input <- match("--output_directory",args)+1
  output_directory <- args[output_directory_input]
  output_directory <- format_output_directory(output_directory)
  return(output_directory)
}
}





