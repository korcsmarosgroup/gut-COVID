# R script to create heatplots of receptors vs immune cells, showing the total number of interactions
#
# Author: Agatha Treveil, January 2021
#
# Input: ligand-receptor stats 2 files output from ligand_receptor_stats.R. tab delimited.
#        Columns: ligand	receptor	ligand_cell_source	ligand_dir	target_cell_counts	receptor_cell_target
#
# Output: Heatplots per epithelial cell type (split up and downregulated). The same but only receptors on >5 immune cells.

# Capture  messages and errors to a file.
zz <- file("gutcovid.out", open = "a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: receptor_imm_cell_heatplots.R\n")

##### Set up #####

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 6){
  stop("Wrong number of command line input parameters. Please check.")
}

library(tidyverse)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(viridis)
#setwd('/Volumes/group-gc/Projects/CovidOmix/stream4_gutcovid/intercellular_networks')

# Input files - infected and bystander cells
lr_sep_inf <- read.csv(args[1], sep = "\t")
lr_sep_by <- read.csv(args[2], sep = "\t")

# Outfile for the normal heatplot
out_f_inf <- args[3]
out_f_by <- args[4]

# Outfile for the smaller heatplot
out_f_inf2 <- args[5]
out_f_by2 <- args[6]

##### Preprocess #####

preproc_data <- function(indata1) {
  
  # Expand receptor/immune cell types
  expand_data <- indata1 %>% separate_rows(receptor_cell_target, sep = ",")
  
  # Count number of interactions by collapsing ligands
  expand_data2 <- expand_data %>% group_by(receptor, receptor_cell_target, target_cell_counts, ligand_dir, ligand_cell_source) %>% 
    summarise(n_ints = n())
  
  return(expand_data2)
}

hplot <- function(indata2, out_f){
  
  # Split up and downregulated interactions
  indata_spl <- split(indata2, indata2$ligand_dir)
  
  # Iterate split tables
  for (name1 in names(indata_spl)){
    
    # get dataframe
    table1 <- indata_spl[[name1]]
    
    # # Get number of facets
    # num_f <- length(unique(table1$ligand_cell_source))
    # 
    # # heatmap without clustering with facet over epithelial cells
    # gg <- ggplot(table1, aes(x=receptor, y=receptor_cell_target, fill=n_ints)) +
    #   geom_tile()+
    #   #scale_fill_brewer(palette = "Purples") +
    #   scale_fill_gradientn(name="# interactions",colours=c("#DADAEB", "#4A1486"), values=scales::rescale(c(1,7))) +
    #   coord_fixed(ratio=0.7) +
    #   facet_wrap(~ligand_cell_source, nrow=num_f) +
    #   labs(title="Number of interactions\n") +
    #   theme(panel.border=element_blank()) +
    #   theme(legend.title.align=1) +
    #   theme(legend.position="bottom") +
    #   theme(text = element_text(size=5), axis.text.x = element_text(size = 4, angle=45, hjust=1), axis.ticks = element_blank(),
    #        axis.text.y = element_text(size=4)) +
    #   theme(strip.text.x = element_text(size = 5)) +
    #   theme(panel.background = element_rect(fill = "#ffffff"))+
    #  theme(panel.grid.major = element_blank(), panel.grid.minor = element_line())
    # 
    # # Create output filepath
    # out_p <- paste0(out_f, "facet_ileum_heatplot_receptor_vs_immcell_", name1, ".png")
    # 
    # # Save facet plot
    # ggsave(gg, filename = out_p, width=10, height =6)
    # 
    # Heatmap with clustering and without facet 
    
    # Iterate epithelial cell types
    table_spl <- split(table1, table1$ligand_cell_source)
    for (table2 in table_spl){
      
      # Reshape dataframe
      dat <- table2 %>% ungroup() %>% select(c(receptor, receptor_cell_target, n_ints)) %>% spread(receptor, n_ints, fill = 0) %>%
        column_to_rownames(var="receptor_cell_target")
     
      # Specify colours - purple gradient
      #bk1 <- c(0,0.4)
      #bk2 <- c(seq(0.5,7, by = 0.5))
      #bk <- c(bk1, bk2)
      #my_palette <- c(colorRampPalette(colors = c("white"))(n = length(bk1)),
      #                c(colorRampPalette(colors = c("#DADAEB", "#4A1486"))(n = length(bk2))))
      
      # Specify colours - intervals
      bk <- c(0,0.9,1.9,2.9,3.9,4.9,5.9,7)
      my_palette <- c('#4575b4','#91bfdb','#e0f3f8','#ffffbf','#fee090','#fc8d59','#d73027')
      
      # Create output filepath
      out_p <- paste0(out_f, "heatplot_receptor_vs_immcell_", table2$ligand_cell_source[1], "_", name1, ".png")
      
      # heatplot
      if (ncol(dat) <3) {
        clus_c = F
      } else {
        clus_c = T
      }
      if (nrow(dat) <3) {
        clus_r = F
      } else {
        clus_r = T
      }
      
      pheatmap(as.matrix(dat),
                      angle_col = 45,
                      color = my_palette,
                      breaks = bk,
                      cluster_rows = clus_r,
                      cellheight = 8,
                      cellwidth = 10,
                      border_color = "darkgrey",
                      fontsize = 8,
                      cluster_cols =clus_c,
                      treeheight_row = 0,
                      treeheight_col = 0,
                      main = table2$ligand_cell_source[1],
                      filename = out_p)
      
    }
  }
}

# Call function to preprocess data
lr_sep_inf2 <- preproc_data(lr_sep_inf)
lr_sep_by2 <- preproc_data(lr_sep_by)

# Call function to plot
hplot(lr_sep_inf2, out_f_inf)
hplot(lr_sep_by2, out_f_by)

# Filter for receptors targeting <5 cells
lr_sep_inf3 <- lr_sep_inf2 %>% filter(target_cell_counts <5)
lr_sep_by3 <- lr_sep_by2 %>% filter(target_cell_counts <5)

# Call function to plot <5
hplot(lr_sep_inf3, out_f_inf2)
hplot(lr_sep_by3, out_f_by2)
