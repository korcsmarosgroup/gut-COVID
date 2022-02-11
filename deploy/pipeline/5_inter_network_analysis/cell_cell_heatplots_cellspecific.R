# R script to create heatplots of ligand-receptor interactions - focusing on the number of lig-rec interaction per cell-cell pair
# Plotting only the immature enterocytes 2 data on two seperate plots (one for ileum and one for colon). Squared heatplot only (no moon plot)
#
# Author: Agatha Treveil, November 2020
#
# Input: ligand-receptor stats files output from ligand_receptor_stats.R. tab delimited. Seperate up and downregulation of epithelial DEGs.: 
#           one file for infected epithelial cells and one for bystander cells
#        Columns: ligand_cell_source	receptor_cell_target	ligand_dir	uniq_ints	uniq_receptors	uniq_ligands
#
# Output: Heatplot

# Capture  messages and errors to a file.
zz <- file("gutcovid.out", open = "a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: cell_cell_heatplots_cellspecific.R\n")

##### Set up #####

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 4){
  stop("Wrong number of command line input parameters. Please check.")
}

library(tidyverse)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(gggibbous)
library(ggplot2)

#setwd("/Users/treveila/Documents/covid/temp_intercellular")
#setwd("/Volumes/group-gc/Projects/CovidOmix/stream4_gutcovid/intercellular_networks")

# Input files - infected and bystander cells
lr_sep_inf<- read.csv(args[1], sep = "\t")
lr_sep_by <- read.csv(args[2], sep = "\t")

# Outfile for the normal heatplot
out1 <- args[3]

# Column order for output plot
col_order <- c("imm_entero2_inf_u", "imm_entero2_inf_d", "imm_enterocyte_2_byst_u") # ileum

if (args[4] == "colon"){
  col_order <- c("imm_enterocyte_2_inf_u", "imm_enterocyte_2_inf_d", "imm_enterocyte_2_byst_u","imm_enterocyte_2_byst_d") # colon
}

# List of names of cells to plot
cells <- c("imm_entero2_inf_u","imm_enterocyte_2_inf_u","imm_entero2_inf_d","imm_enterocyte_2_inf_d","imm_enterocyte_2_byst_u","imm_enterocyte_2_byst_d")

##### Preprocess #####

# Join the infected and uninfected cell types together
lr_sep_by <- lr_sep_by %>% mutate(ligand_cell_source = paste0(ligand_cell_source, "_byst")) %>% mutate(status = "bystander")
lr_sep_inf <- lr_sep_inf %>% mutate(ligand_cell_source = paste0(ligand_cell_source,"_inf")) %>% mutate(status = "infected")
lr_sep <- rbind(lr_sep_inf, lr_sep_by) %>%
  mutate(ligand_cell_source = str_replace(ligand_cell_source, "enterocytes_","entero")) %>%
  mutate(receptor_cell_target = str_replace(receptor_cell_target, "_cell","")) %>%
  mutate(ligand_cell_source = ifelse(ligand_dir == "upreg", paste0(ligand_cell_source, "_u"), paste0(ligand_cell_source, "_d"))) %>%
  arrange(desc(ligand_dir), status, ligand_cell_source)

# Filter for cells of interest
lr_sep <- lr_sep %>% filter(ligand_cell_source %in% cells)

# Order epithelial cells
lr_sep <- lr_sep %>% arrange(match(ligand_cell_source, col_order))

# Spread data
data <- lr_sep %>% select(ligand_cell_source, receptor_cell_target, uniq_ints) %>%
  acast(receptor_cell_target~ligand_cell_source, value.var="uniq_ints")

# Reorder columns
data <- data[,unique(lr_sep$ligand_cell_source)]

##### Squared heatplot #####

# Annotation ligand dir
ann1 <- lr_sep %>% select(c(ligand_cell_source, ligand_dir,status )) %>% unique()
rownames(ann1) <- ann1$ligand_cell_source 
ann1 <- ann1 %>% select(-c(ligand_cell_source))

# Remove "_inf" and "byst" from cell names
#new_cols <- str_replace(colnames(data),c("_inf_u|_byst_u|_inf_d|_byst_d"), "")

# Remove imm entoercyte label from cell names
new_cols <- str_replace(colnames(data),c("imm_entero2_|imm_enterocyte_2_"), "")
new_cols <- str_replace_all(new_cols,c("inf"= "infected", "byst"="bystander","_u"="_upreg","_d"="_downreg"))

# Plot breaks so colour same between different plots
breaksList = seq(0, 275, by = 1)

# Plot
x=pheatmap(as.matrix(data), 
          color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(length(breaksList)),
          cluster_rows = T,
          cluster_cols =F,
          fontsize = 9,
          angle_col = 45,
          #gaps_col = 6,#7
          annotation_col = ann1,
          breaks = breaksList,
          main = "Imm enterocyte 2 intercellular interactions", #ileal, colon
          labels_col= c(new_cols))

ggsave(x, filename = out1, width = 4.4, height = 5)

