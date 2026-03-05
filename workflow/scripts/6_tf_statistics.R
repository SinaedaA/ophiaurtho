#!/usr/bin/env Rscript
###################
#### LIBRARIES ####
###################
# Use pacman to load libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr,
  tidyr,
  ggplot2,
  viridis, 
  argparser,
  rio
)

#########################
#### PARSE ARGUMENTS ####
#########################
p <- arg_parser("Compute statistics for all TFs")
p <- add_argument(p, "--tf_classification", help = "Relative path to the table containing the classification of TFs.")
p <- add_argument(p, "--cog2seqid", help = "Relative path to the COG to SEQID conversion (a table containing all the seqids containined in each COG).")
p <- add_argument(p, "--outdir", help = "Relative path to the output directory.")

argv <- parse_args(p)

###################
#### LOAD DATA ####
###################
tf_class <- rio::import(argv$tf_classification, header = TRUE)
cog2seqid <- rio::import(argv$cog2seqid, header = TRUE)

##################################
#### STATS FOR EACH TF FAMILY ####
###################################
tf_fams <- unique(tf_class$Family_name)
full_stats <- plyr::ldply(tf_fams, function(x){
    seqids <- tf_class[tf_class$Family_name == x, "geneid"]
    cogs_x <- dplyr::filter(cog2seqid, id %in% seqids) %>% filter(!is.na(cog)) %>% pull(cog) %>% unique
    stats_cog_x <- dplyr::filter(cog2seqid, cog %in% cogs_x) %>% 
		dplyr::mutate(Part_of_family = ifelse(id %in% seqids, 1, 0)) %>% 
		dplyr::group_by(cog) %>% 
		dplyr::mutate(n_Identified = sum(Part_of_family), n_Genes_in_COG = n(), Percentage = (n_Identified/n_Genes_in_COG)*100, Family = x) %>% 
		dplyr::select(-id, -Part_of_family) %>% unique() 
    return(stats_cog_x)
})
write.table(full_stats, paste0(argv$outdir, "/tf_full_stats.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

########################
#### PLOT THIS DATA ####
########################
mean_perc_identified <- full_stats %>% dplyr::group_by(Family) %>% dplyr::summarise(mean_perc = mean(Percentage))
mean_n_genes <- full_stats %>% dplyr::group_by(Family) %>% dplyr::summarise(mean_n_genes = mean(n_Genes_in_COG))

mean_perc_id_plot <- ggplot(mean_perc_identified, aes(x = reorder(Family, mean_perc), y = mean_perc, fill = Family)) +
    geom_bar(stat = 'identity', alpha = .5) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(
        x = "TF Family",
        y = "Mean percentage of genes identified as TFs in COGs",
        title = "Mean percentage of genes identified in COGs for each TF family"
    ) 
ggsave(paste0(argv$outdir, "/Mean_percentage_COG_genes_identified_TFs.pdf"), mean_perc_id_plot, height = 8, width = 18)
