# R script to create heatplots of ligand-receptor interactions - focusing on the number of cells each ligand-rec interaction are in
#
# Author: Agatha Treveil, November 2020
#
# Input: ligand-receptor stats 2 files output from ligand_receptor_stats.R. tab delimited.
#        Columns: ligand	receptor	ligand_cell_source	ligand_dir	target_cell_counts	receptor_cell_target
#
# Output: 

##### Set up #####

library(tidyverse)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(viridis)

#setwd("/Users/poletti/OneDrive\ -\ Norwich\ BioScience\ Institutes/Bioinformatics/Lung_run/ciliated_cells_moderate")

# Capture  messages and errors to a file.
zz <- file("gutcovid_lung.out", open="a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: 14_ligand_receptor_heatplots.R\n")

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 2){
  stop("Wrong number of command line input parameters. Please check.")
}

# Input files - infected and bystander cells
lr_sep_inf<- read.csv(args[1], sep = "\t")

# Outfile for the normal heatplot
out_f_inf <- args[2]

##### Preprocess #####

hplot <- function(indata, out_f){
  
  # Split up and downregulated interactions
  indata_spl <- split(indata, indata$ligand_dir)
  
  # Iterate split tables
  for (name1 in names(indata_spl)){
    
    # get dataframe
    table1 <- indata_spl[[name1]]
    
    # Get number of facets
    num_f <- length(unique(table1$ligand_cell_source))
  
    # heatmap without clustering with facet over epithelial cells
    gg <- ggplot(table1, aes(x=ligand, y=receptor, fill=target_cell_counts)) +
      geom_tile()+
      scale_fill_gradientn(name="# immune cells",colours=c("#EBEBEB", "#DADAEB", "#4A1486"), values=scales::rescale(c(0, 0.01, 1))) +
      coord_fixed(ratio=0.7) +
      facet_wrap(~ligand_cell_source, ncol=num_f) +
      labs(title="Number of immune cell targets for each epithelial cell ligand\n") +
      theme(panel.border=element_blank()) +
      theme(legend.title.align=1) +
      theme(legend.position="bottom") +
      theme(text = element_text(size=6), axis.text.x = element_text(size = 5, angle=45, hjust=1), axis.ticks = element_blank(),
           axis.text.y = element_text(size=5)) +
      theme(strip.text.x = element_text(size = 6)) +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_line())
    
    # Create output filepath
    out_p <- paste0(out_f, "facet_heatplot_num_cells_", name1, ".png")
    
    # Save facet plot
    ggsave(gg, filename = out_p)
    
    # Heatmap with clustering and without facet 
    # Empty list of plots
    #plot_list <- list()
    # Iterate epithelial cell types
    table_spl <- split(table1, table1$ligand_cell_source)
    for (table2 in table_spl){
      
      # Reshape dataframe
      dat <- table2 %>% select(c(ligand, receptor, target_cell_counts)) %>% spread(ligand, target_cell_counts, fill = 0) %>%
        column_to_rownames(var="receptor")
     
      # Specify colours
      bk1 <- c(0,0.1)
      bk2 <- c(seq(0.2,12, by = 0.1))
      bk <- c(bk1, bk2)
      my_palette <- c(colorRampPalette(colors = c("white"))(n = length(bk1)),
                      c(colorRampPalette(colors = c("#DADAEB", "#4A1486"))(n = length(bk2))))
      
      # Create output filepath
      out_p <- paste0(out_f, "ileum_heatplot_num_cells_", table2$ligand_cell_source[1], "_", name1, ".png")
      
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
                      cellheight = 6,
                      cellwidth = 8,
                      border_color = "grey",
                      fontsize = 6,
                      cluster_cols =clus_c,
                      main = table2$ligand_cell_source[1],
                      filename = out_p)
      
    }
  }
}

# Call function to plot
hplot(lr_sep_inf, out_f_inf)
#hplot(lr_sep_by, out_f_by)
