# R script to create heatplots of ligands vs immune cells, showing the sum of receptor expression
#
# Author: Agatha Treveil, January 2021
#
# Input: ligand-imm cell interaction files output from ligand_imm_cell_network_Creation.R. tab delimited. Already filtered for ligand cell type.

# Output: Heatplots per epithelial cell type (split up and downregulated).only for one ligand cell type - no facet.

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
message("\nStarting differential expression analysis script: 12_ligand_imm_cell_heatplots_receptor_expression.R\n")

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 2){
  stop("Wrong number of command line input parameters. Please check.")
}

# Input files - infected cells (specific cell type already filtered)
lr_sep<- read.csv(args[1], sep = "\t")

# Name of cell type (and infected or bystander) for output filename
ligand_celltype <- "ciliated_cells"

# Outfile for the normal heatplot
out_f <- args[2]

##### Preprocess #####

# Add ligand dir column
lr_sep <- lr_sep %>% mutate(ligand_dir = if_else(ligand_lfc >0, "upreg", "downreg"))

# PLotting function
hplot <- function(indata2, out_f){
  
  # Split up and downregulated interactions
  indata_spl <- split(indata2, indata2$ligand_dir)
  
  # Iterate split tables
  for (name1 in names(indata_spl)){
    
    # get dataframe
    table1 <- indata_spl[[name1]]
    
    # Heatmap with clustering and without facet 
      
    # Reshape dataframe
    dat <- table1 %>% ungroup() %>% select(c(ligand, receptor_cell_target, sum_receptor_exp)) %>% spread(ligand, sum_receptor_exp, fill = 0) %>%
        column_to_rownames(var="receptor_cell_target")
     
    ## Specify colours
    #bk1 <- c(0,0.4)
    #bk2 <- c(seq(0.5,85, by = 0.5))
    #bk <- c(bk1, bk2)
    #my_palette <- c(colorRampPalette(colors = c("white"))(n = length(bk1)),
                    #c(colorRampPalette(colors = c("#DADAEB", "#4A1486"))(n = length(bk2))))
      
    # Specify colours - intervals
    bk <- c(0,0.9,1.9,2.9,3.9,4.9,5.9,7)
    my_palette <- c('#4575b4','#91bfdb','#e0f3f8','#ffffbf','#fee090','#fc8d59','#d73027')
    
    # Create output filepath
    out_p <- paste0(out_f, "1_heatplot_ligand_vs_immcell_", ligand_celltype, "_", name1, ".png")
    
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
                    border_color = "darkgrey",
                    fontsize = 6,
                    cluster_cols =clus_c,
                    treeheight_row = 0,
                    treeheight_col = 0,
                    main = paste0(ligand_celltype, " sum of receptor expression "),
                    filename = out_p)

  }
}


# Call function to plot
hplot(lr_sep, out_f)
