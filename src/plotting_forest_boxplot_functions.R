#plotting functions for qtlenrich_plotting_module

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(cowplot))

#### PLOTTING FUNCTIONS ####
parse_alternate_gtex_tissues <- function(args){
  if(is.na(match("--gtex_tissue_names",args))){
    gtex_tissue_names <- NA
    return(gtex_tissue_names)
  }else{
    gtex_tissues_file_input <- match("--gtex_tissue_names",args)+1
    gtex_tissues_file <- args[gtex_tissues_file_input]
    gtex_tissue_names <- suppressMessages(read_tsv(gtex_tissues_file)) %>%
      select(tissue_id,tissue_site_detail) %>% 
      rename(`Tissue_(QTL)` = tissue_id)
    return(gtex_tissue_names)
   }
}

merge_with_alternate_gtex_tissues <- function(tissue_names,qtlenrich_result){
  qtlenrich_new_tissue_names <- inner_join(tissue_names,qtlenrich_result,by="Tissue_(QTL)")
  return(qtlenrich_new_tissue_names)
}

parse_qtlenrich_file <- function(qtlenrich_file1,label_file1,gtex_tissue_names,significant_cutoff,number_significant_tissues,tissues_list_to_remove){

  qtlenrich_file1_parsed <- suppressMessages(read_tsv(qtlenrich_file1)) %>%
    mutate(enrichment_p_value_log_10 = -1*log10(`Enrichment_P_value`),
           type = label_file1)

  qtlenrich_file1_parsed$enrichment_p_value_log_10 <- round(qtlenrich_file1_parsed$enrichment_p_value_log_10,digits=2)

  if(is.data.frame(gtex_tissue_names)){
    qtlenrich_file1_parsed <- merge_with_alternate_gtex_tissues(gtex_tissue_names,qtlenrich_file1_parsed)
  }

  if(significant_cutoff == "bonferroni"){
    significant_cutoff <- 0.05/(nrow(qtlenrich_file1_parsed))
  } 

  qtlenrich_file1_parsed <- qtlenrich_file1_parsed %>%
    filter(`Enrichment_P_value` <= significant_cutoff)

  #define negation
  `%!in%` = Negate(`%in%`)

  if(is.data.frame(tissues_list_to_remove)){
    qtlenrich_file1_parsed <- qtlenrich_file1_parsed %>%
      filter(`Tissue_(QTL)` %!in% tissues_list_to_remove$`Tissue_(QTL)`)
  }


  if(!is.na(number_significant_tissues)){
   qtlenrich_file1_parsed <- qtlenrich_file1_parsed %>%
     top_n(n=number_significant_tissues,wt = `Adjusted_Fold_Enrichment`)
  }

  #sort by adjusted fold-enrichment
  qtlenrich_file1_parsed <- qtlenrich_file1_parsed %>%
   mutate(adjusted_sorted = reorder(qtlenrich_file1_parsed$"tissue_site_detail",-qtlenrich_file1_parsed$"Adjusted_Fold_Enrichment"))
  return(qtlenrich_file1_parsed)
}

parse_qtlenrich_file2 <- function(qtlenrich_file1_parsed,qtlenrich_file2,label_file2,gtex_tissue_names){

  qtlenrich_file2_parsed <- suppressMessages(read_tsv(qtlenrich_file2)) %>%
    mutate(enrichment_p_value_log_10 = -1*log10(`Enrichment_P_value`),
           type = label_file2)

  qtlenrich_file2_parsed$enrichment_p_value_log_10 <- round(qtlenrich_file2_parsed$enrichment_p_value_log_10,digits=2)

  if(length(gtex_tissue_names) != 0){
    qtlenrich_file2_parsed <- merge_with_alternate_gtex_tissues(gtex_tissue_names,qtlenrich_file2_parsed)
  }

  qtlenrich_file2_parsed <- qtlenrich_file2_parsed %>%
      filter(`Tissue_(QTL)` %in% qtlenrich_file1_parsed$`Tissue_(QTL)`)

  #sort by adjusted fold-enrichment
  qtlenrich_file2_parsed <- qtlenrich_file2_parsed %>%
   mutate(adjusted_sorted = reorder(qtlenrich_file2_parsed$"tissue_site_detail",-qtlenrich_file2_parsed$"Adjusted_Fold_Enrichment"))
  return(qtlenrich_file2_parsed)
}

generate_single_forest_plot <- function(qtlenrich_file1_parsed){

  if (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 10){
    font_size <- 12
   } else if ((length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) > 10) & (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 35)) {
    font_size <- 24
   } else {
     font_size <- 16
   }

  max_limit <- max(qtlenrich_file1_parsed$`Upper_bound_95%_CI`)+0.1

  forest_plot <- ggplot(qtlenrich_file1_parsed,
       aes(x=reorder(tissue_site_detail,`Adjusted_Fold_Enrichment`),y=`Adjusted_Fold_Enrichment`))+
  geom_point(aes(size=enrichment_p_value_log_10))+
  geom_pointrange(aes(ymin = `Lower_bound_95%_CI`, ymax = `Upper_bound_95%_CI`))+
  geom_hline(yintercept = 1,linetype="f8")+
  ylab("Adjusted Fold-Enrichment")+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.line.y=element_blank(),
        axis.text.y=element_text(size=font_size), #24
        axis.text.x=element_text(size=font_size),
        axis.title.x=element_text(size=font_size),
        legend.title = element_text(size=font_size),
        legend.text=element_text(size=font_size),
        legend.spacing.x=unit(0.5, "cm"),
        legend.spacing.y=unit(0.5, "cm"))+
  guides(size = guide_legend(title="-log10(P-value)"))+
  ylim(0.8,max_limit)+
  coord_flip()

}

generate_paired_forest_plot <- function(qtlenrich_files,label_file1,label_file2){

  if (length(unique(qtlenrich_files$`Tissue_(QTL)`)) <= 5){
     font_size <- 18
   }else if ( (length(unique(qtlenrich_files$`Tissue_(QTL)`)) > 5) & (length(unique(qtlenrich_files$`Tissue_(QTL)`)) <= 20)){
    font_size <- 24
   } else if ((length(unique(qtlenrich_files$`Tissue_(QTL)`)) > 20) & (length(unique(qtlenrich_files$`Tissue_(QTL)`)) <= 35)) {
    font_size <- 24 
   } else {
     font_size <- 30
   }


  qtlenrich_files$type <- factor(qtlenrich_files$type,levels=c(label_file2,label_file1)) 
 
  num_rows <- length(unique(qtlenrich_files$`Tissue_(QTL)`))

  max_limit <- max(qtlenrich_files$`Upper_bound_95%_CI`)+0.1

  forest_plot <- ggplot(qtlenrich_files,
       aes(x=type,y=`Adjusted_Fold_Enrichment`,col=type))+ #size = enrichment_p_value_log_10
  geom_point(aes(size=enrichment_p_value_log_10))+
  geom_pointrange(aes(ymin = `Lower_bound_95%_CI`, ymax = `Upper_bound_95%_CI`))+
  ylab("Adjusted Fold-Enrichment")+
  facet_wrap(~adjusted_sorted,strip.position = "left",nrow=num_rows,scales="free_y")+
  theme_classic()+
  scale_color_manual(breaks= c(label_file1,label_file2),
                     labels=c(label_file1,label_file2),
                     values = c("gray","black"))+
  guides(size = guide_legend(title="-log10(P)"),
         col = guide_legend(title="QTL",override.aes=list(size=1.8)))+
  geom_hline(yintercept = 1,linetype="f8")+
  theme(axis.text.x=element_text(size=font_size), #40 to 32 (minus 8)
        axis.title.x=element_text(size=font_size),
        legend.title = element_text(size=font_size),
        legend.text=element_text(size=font_size,margin = margin(r = 2, unit = 'in')),
        legend.spacing.x=unit(0.5, "cm"),
        legend.spacing.y=unit(0.5, "cm"),
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(2, "cm"),
        axis.line.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title=element_text(size=font_size),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,size=font_size),
        panel.spacing = unit(0, "mm"),
        strip.background = element_blank())+
  ylim(0.8,max_limit)+
  coord_flip()
  return(forest_plot)
}

generate_single_bar_plot <- function(qtlenrich_file1_parsed,adjusted_sorted_2,bar_box_violin_header,bar_box_violin_axes_label){

  if (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 10){
    font_size <- 12
   } else if ((length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) > 10) & (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 35)) {
    font_size <- 24
   } else {
     font_size <- 16
   }

  bar_plot <- ggplot(qtlenrich_file1_parsed,
                     aes_string(adjusted_sorted_2,bar_box_violin_header))+
  geom_bar(stat="identity",fill="#56B4E9")+
  ylab(bar_box_violin_axes_label)+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size=font_size),
        axis.text.x = element_text(size=font_size-4),
        axis.title.x = element_text(size=font_size-2))+
  coord_flip()
}

generate_paired_bar_plot <- function(qtlenrich_files,label_file1,label_file2,type,bar_box_violin_header,bar_box_violin_axes_label){
  qtlenrich_files$type <- factor(qtlenrich_files$type,levels=c(label_file2,label_file1))

  if (length(unique(qtlenrich_files$`Tissue_(QTL)`)) <= 2){
    font_size <- 26
   } else if ((length(unique(qtlenrich_files$`Tissue_(QTL)`)) > 2) & (length(unique(qtlenrich_files$`Tissue_(QTL)`)) <= 35)) {
    font_size <- 32
   } else {
     font_size <- 16
   }

  num_rows <- length(unique(qtlenrich_files$`Tissue_(QTL)`))

  paired_bar_plot <- ggplot(qtlenrich_files,
       aes_string(x=type,y=bar_box_violin_header,fill=type))+
   geom_bar(stat = "identity",position=position_dodge(),width=0.5)+
   ylab(bar_box_violin_axes_label)+
   facet_wrap(~adjusted_sorted,strip.position = "left",nrow=num_rows,scales="free_y")+
   theme_classic()+
   theme(axis.text.x=element_text(size=font_size),
        axis.title.x=element_text(size=font_size),
        legend.title = element_text(size=font_size),
        legend.text=element_text(size=font_size),
        legend.spacing.x=unit(0.5, "cm"),
        legend.spacing.y=unit(0.5, "cm"),
        axis.line.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title=element_text(size=font_size),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,size=font_size),
        panel.spacing = unit(0, "mm"),
        strip.background = element_blank())+
   guides(fill = guide_legend(title="QTL"))+
   scale_fill_manual(breaks= c(label_file1,label_file2),
                     labels=c(label_file1,label_file2),
                     values = c("gray","black"))+
   expand_limits(y=1)+
   coord_flip()
  return(paired_bar_plot)
}

generate_single_forest_plot_side_by_side <- function(qtlenrich_file1_parsed){

  if (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 5){
    font_size <- 8
   }
   else if ((length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) > 5) & (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 10))  {
    font_size <- 12
   } else if ((length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) > 10) & (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 35)) {
    font_size <- 24
   } else {
     font_size <- 16
   }

  max_limit <- max(qtlenrich_file1_parsed$`Upper_bound_95%_CI`)+0.1

  forest_plot <- ggplot(qtlenrich_file1_parsed,
       aes(x=reorder(tissue_site_detail,`Adjusted_Fold_Enrichment`),y=`Adjusted_Fold_Enrichment`))+
  geom_point(aes(size=enrichment_p_value_log_10))+
  geom_pointrange(aes(ymin = `Lower_bound_95%_CI`, ymax = `Upper_bound_95%_CI`))+
  ylab("Adjusted Fold-Enrichment")+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.line.y = element_blank(),
        axis.text.y=element_text(size=font_size),
        axis.text.x=element_text(size=font_size),
        axis.title.x=element_text(size=font_size),
        legend.title = element_text(size=font_size),
        legend.text=element_text(size=font_size),
        legend.spacing.x=unit(0.5, "cm"),
        legend.spacing.y=unit(0.5, "cm"))+
  guides(size = guide_legend(title="-log10(P-value)"))+
  geom_hline(yintercept = 1,linetype=5)+
  ylim(0.8,max_limit)+
  coord_flip()
}

generate_single_bar_plot_side_by_side <- function(qtlenrich_file1_parsed,adjusted_sorted_2,bar_box_violin_header,bar_box_violin_axes_label){

  if (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 5){
    font_size <- 8
    width_size <- 0.5
   }
   else if ((length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) > 5) & (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 10))  {
    font_size <- 12
    width_size <- 0.8
   } else if ((length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) > 10) & (length(unique(qtlenrich_file1_parsed$`Tissue_(QTL)`)) <= 35)) {
    font_size <- 24
    width_size <- 0.8
   } else {
     font_size <- 16
     width_size <- 0.8
   }

  bar_plot <- ggplot(qtlenrich_file1_parsed,
                     aes_string(adjusted_sorted_2,bar_box_violin_header))+
  geom_bar(stat="identity",fill="#56B4E9",width=width_size)+
  ylab(bar_box_violin_axes_label)+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=font_size),
        axis.title.x = element_text(size=font_size))+
  coord_flip()
}

generate_unpaired_plot_side_by_side <- function(qtlenrich_file1_parsed,adjusted_sorted_2,bar_box_violin_header,bar_box_violin_axes_label){

  forest_plot <- generate_single_forest_plot_side_by_side(qtlenrich_file1_parsed)

  bar_plot <- generate_single_bar_plot_side_by_side(qtlenrich_file1_parsed,adjusted_sorted_2,bar_box_violin_header,bar_box_violin_axes_label)

side_by_side <- plot_grid(forest_plot,bar_plot)
return(side_by_side)
}

generate_paired_bar_plot_side_by_side <- function(qtlenrich_files,label_file1,label_file2,type,bar_box_violin_header,bar_box_violin_axes_label){

  if (length(unique(qtlenrich_files$`Tissue_(QTL)`)) <= 5){
    font_size <- 8
   } else if ( length(unique(qtlenrich_files$`Tissue_(QTL)`)) > 5 &  length(unique(qtlenrich_files$`Tissue_(QTL)`)) <= 20){
    font_size <- 12
   } else {
     font_size <- 24
   }

  qtlenrich_files$type <- factor(qtlenrich_files$type,levels=c(label_file2,label_file1))

  num_rows <- length(unique(qtlenrich_files$`Tissue_(QTL)`))

  paired_bar_plot <- ggplot(qtlenrich_files,
       aes_string(x=type,y=bar_box_violin_header,fill=type))+
   geom_bar(stat = "identity",position=position_dodge(),width=0.7)+
   ylab(bar_box_violin_axes_label)+
   facet_wrap(~adjusted_sorted,strip.position = "left",nrow=num_rows,scales="free_y")+
   theme_classic()+
   theme(axis.text.x=element_text(size=font_size),
        axis.title.x=element_text(size=font_size),
        legend.title = element_text(size=font_size),
        legend.text=element_text(size=font_size),
        legend.spacing.x=unit(0.5, "cm"),
        legend.spacing.y=unit(0.5, "cm"),
        axis.line.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title=element_text(size=font_size),
        strip.text.y = element_blank(),
        panel.spacing = unit(0, "mm"),
        strip.background = element_blank())+
   guides(fill = guide_legend(title="QTL"))+
   scale_fill_manual(breaks= c(label_file1,label_file2),
                     labels=c(label_file1,label_file2),
                     values = c("gray","black"))+
   coord_flip()
  return(paired_bar_plot)
}

generate_paired_forest_plot_side_by_side <- function(qtlenrich_files,label_file1,label_file2){

  if (length(unique(qtlenrich_files$`Tissue_(QTL)`)) <= 5){
    font_size <- 8
   } else if ( length(unique(qtlenrich_files$`Tissue_(QTL)`)) > 5 &  length(unique(qtlenrich_files$`Tissue_(QTL)`)) <= 20){
    font_size <- 12
   } else {
     font_size <- 24
   }

  qtlenrich_files$type <- factor(qtlenrich_files$type,levels=c(label_file2,label_file1))

  num_rows <- length(unique(qtlenrich_files$`Tissue_(QTL)`))

  max_limit <- max(qtlenrich_files$`Upper_bound_95%_CI`)+0.1

  forest_plot <- ggplot(qtlenrich_files,
       aes(x=type,y=`Adjusted_Fold_Enrichment`,col=type))+ #size = enrichment_p_value_log_10
  geom_point(aes(size=enrichment_p_value_log_10))+
  geom_pointrange(aes(ymin = `Lower_bound_95%_CI`, ymax = `Upper_bound_95%_CI`))+
  ylab("Adjusted Fold-Enrichment")+
  facet_wrap(~adjusted_sorted,strip.position = "left",nrow=num_rows,scales="free_y")+
  theme_classic()+
  scale_color_manual(breaks= c(label_file1,label_file2),
                     labels=c(label_file1,label_file2),
                     values = c("gray","black"))+
  guides(size = guide_legend(title="-log10(P)"),
         col = guide_legend(title="QTL",override.aes=list(size=1.8)))+
  geom_hline(yintercept = 1,linetype=5)+
  theme(axis.text.x=element_text(size=font_size),
        axis.title.x=element_text(size=font_size),
        legend.title = element_text(size=font_size),
        legend.text=element_text(size=font_size),
        legend.spacing.x=unit(1, "cm"),
        legend.spacing.y=unit(0.5, "cm"),
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(2, "cm"),
        axis.line.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title=element_text(size=font_size-8),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,size=font_size),
        panel.spacing = unit(0, "mm"),
        strip.background = element_blank())+
  ylim(0.8,max_limit)+
  coord_flip()

  return(forest_plot)
}


generate_paired_plot_side_by_side <- function(qtlenrich_files,label_file1,label_file2,type,bar_box_violin_header,bar_box_violin_axes_label){

forest_plot <- generate_paired_forest_plot_side_by_side(qtlenrich_files,label_file1,label_file2)

bar_plot <- generate_paired_bar_plot_side_by_side(qtlenrich_files,label_file1,label_file2,type,bar_box_violin_header,bar_box_violin_axes_label)

side_by_side <- plot_grid(forest_plot,bar_plot)
return(side_by_side)
}

#function to compute max y-axis limit for box and violin plots
compute_max_limit <- function(qtlenrich_data,bar_box_violin_header){
  if (bar_box_violin_header == "Estimated_total_num_trait_associations"){
    max_limit <- max(qtlenrich_data[,bar_box_violin_header])+20
    return(max_limit)
  }else{
    max_limit <- max(qtlenrich_data[,bar_box_violin_header])+0.1
    return(max_limit)
  }
}
#box plots
generate_single_box_plot <- function(qtlenrich_file1_parsed,label_file1,type,bar_box_violin_header,bar_box_violin_axes_label){

max_limit <- compute_max_limit(qtlenrich_file1_parsed,bar_box_violin_header)

boxplot <- ggplot(qtlenrich_file1_parsed,
       aes_string(x=type,y =bar_box_violin_header,fill=type))+
  geom_boxplot(width=0.4)+
  theme_classic()+
  xlab("QTL")+
  ylab(bar_box_violin_axes_label)+
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10))+
  scale_fill_manual(breaks = c(label_file1),
                     labels= c(label_file1),
                     values = c("#696969"),
                    guide=FALSE)+
  ylim(0,max_limit)
}

generate_paired_box_plot <- function(qtlenrich_files,label_file1,label_file2,type,bar_box_violin_header,bar_box_violin_axes_label){

max_limit <- compute_max_limit(qtlenrich_files,bar_box_violin_header)

qtlenrich_files$type <- factor(qtlenrich_files$type,levels=c(label_file1,label_file2))

boxplot <- ggplot(qtlenrich_files,
       aes_string(x=type,y =bar_box_violin_header,fill=type))+
  geom_boxplot()+
  theme_classic()+
  ylab(bar_box_violin_axes_label)+
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10))+
  scale_fill_manual(breaks = c(label_file1,label_file2),
                     labels= c(label_file1,label_file2),
                     values = c("#696969","#D3D3D3"),
                    guide=FALSE)+
  ylim(0,max_limit)

}

#violin plot
generate_single_violin_plot <- function(qtlenrich_file1_parsed,label_file1,type,bar_box_violin_header,bar_box_violin_axes_label){

max_limit <- compute_max_limit(qtlenrich_file1_parsed,bar_box_violin_header)

violin_plot <- ggplot(qtlenrich_file1_parsed,
       aes_string(x=type,y=bar_box_violin_header,fill=type))+
  geom_violin(width=0.4)+
  stat_summary(fun.y=median, geom="point", size=2, color="black")+
  theme_classic()+
  ylab(bar_box_violin_axes_label)+
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10))+
  scale_fill_manual(breaks = c(label_file1),
                     labels= c(label_file1),
                     values = c("#696969"),
                    guide=FALSE)+
  ylim(0,max_limit)
return(violin_plot)
}


generate_paired_violin_plot <- function(qtlenrich_files,label_file1,label_file2,type,bar_box_violin_header,bar_box_violin_axes_label){

qtlenrich_files$type <- factor(qtlenrich_files$type,levels=c(label_file1,label_file2))
max_limit <- compute_max_limit(qtlenrich_files,bar_box_violin_header)

violinplot <- ggplot(qtlenrich_files,
       aes_string(x=type,y=bar_box_violin_header,fill=type))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(width=.1)+
  #stat_summary(fun.y=median, geom="point", size=2, color="black")+
  theme_classic()+
  ylab(bar_box_violin_axes_label)+
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10))+
  scale_fill_manual(breaks = c(label_file1,label_file2),
                     labels= c(label_file1,label_file2),
                     values = c("#696969","#D3D3D3"),
                    guide=FALSE)+
  ylim(0,max_limit)

}

#output filenames

create_single_forest_plot_filename <- function(output_directory,label_file1,exp_label,date){
  filename <- paste(output_directory,"QTLEnrich",sep="")
  if (exp_label == ""){
  filename <- paste(filename,label_file1,"forest_plot",date,sep="_")}
  else{
  filename <- paste(filename,label_file1,exp_label,"forest_plot",date,sep="_")
}
  filename <- paste(filename,"pdf",sep=".")
  return(filename)
}

create_paired_forest_plot_filename <- function(output_directory,label_file1,label_file2,exp_label,date){
  filename <- paste(output_directory,"QTLEnrich",sep="")
  if (exp_label == ""){
  filename <- paste(filename,label_file1,label_file2,"paired_forest_plot",date,sep="_")}
  else{
  filename <- paste(filename,label_file1,label_file2,exp_label,"paired_forest_plot",date,sep="_")}
  filename <- paste(filename,"pdf",sep=".")
  return(filename)
}

create_single_bar_plot_filename <- function(output_directory,label_file1,exp_label,date){
  filename <- paste(output_directory,"QTLEnrich",sep="")
  if (exp_label == ""){
  filename <- paste(filename,label_file1,"bar_plot",date,sep="_")}
  else{
  filename <- paste(filename,label_file1,exp_label,"bar_plot",date,sep="_")
}
  filename <- paste(filename,"pdf",sep=".")
  return(filename)
}

create_paired_bar_plot_filename <- function(output_directory,label_file1,label_file2,exp_label,date){
  filename <- paste(output_directory,"QTLEnrich",sep="")
  if (exp_label == ""){
  filename <- paste(filename,label_file1,label_file2,"paired_bar_plot",date,sep="_")}
  else{
  filename <- paste(filename,label_file1,label_file2,exp_label,"paired_bar_plot",date,sep="_")
}
  filename <- paste(filename,"pdf",sep=".")
  return(filename)
}

create_unpaired_side_by_side_plot_filename <- function(output_directory,label_file1,exp_label,date){

  filename <- paste(output_directory,"QTLEnrich",sep="")
  if (exp_label == ""){
  filename <- paste(filename,label_file1,"unpaired_side_by_side_plot",date,sep="_")}
  else{
  filename <- paste(filename,label_file1,exp_label,"unpaired_side_by_side_plot",date,sep="_")
}
  filename <- paste(filename,"pdf",sep=".")
  return(filename)
}

create_paired_side_by_side_plot_filename <- function(output_directory,label_file1,label_file2,exp_label,date){
  filename <- paste(output_directory,"QTLEnrich",sep="")
   if (exp_label == ""){
  filename <- paste(filename,label_file1,label_file2,"paired_side_by_side_plot",date,sep="_")}
  else{
  filename <- paste(filename,label_file1,label_file2,exp_label,"paired_side_by_side_plot",date,sep="_")
}
  filename <- paste(filename,"pdf",sep=".")
  return(filename)
}

create_single_box_plot_filename <- function(output_directory,label_file1,exp_label,date){
  filename <- paste(output_directory,"QTLEnrich",sep="")
  if (exp_label == ""){
  filename <- paste(filename,label_file1,"box_plot",date,sep="_")}
  else{
  filename <- paste(filename,label_file1,exp_label,"box_plot",date,sep="_")
}
  filename <- paste(filename,"pdf",sep=".")
  return(filename)
}

create_paired_box_plot_filename <- function(output_directory,label_file1,label_file2,exp_label,date){
  filename <- paste(output_directory,"QTLEnrich",sep="")
  if (exp_label == ""){
  filename <- paste(filename,label_file1,label_file2,"paired_box_plot",date,sep="_")}
  else{
  filename <- paste(filename,label_file1,label_file2,exp_label,"paired_box_plot",date,sep="_")
}
  filename <- paste(filename,"pdf",sep=".")
  return(filename)
}

create_single_violin_plot_filename <- function(output_directory,label_file1,exp_label,date){
  filename <- paste(output_directory,"QTLEnrich",sep="")
  if (exp_label == ""){
  filename <- paste(filename,label_file1,"violin_plot",date,sep="_")}
  else{
  filename <- paste(filename,label_file1,exp_label,"violin_plot",date,sep="_")
}
  filename <- paste(filename,"pdf",sep=".")
  return(filename)
}

create_paired_violin_plot_filename <- function(output_directory,label_file1,label_file2,exp_label,date){
  filename <- paste(output_directory,"QTLEnrich",sep="")
  if (exp_label == ""){
  filename <- paste(filename,label_file1,label_file2,"paired_violin_plot",date,sep="_")}
  else{
  filename <- paste(filename,label_file1,label_file2,exp_label,"paired_violin_plot",date,sep="_")
}
  filename <- paste(filename,"pdf",sep=".")
  return(filename)
}

