# R script to create statistics summaries of ligand-receptor interactions
#
# Author: Agatha Treveil, November 2020
#
# Input: ligand-receptor file output from epithelial_organoid_immune_cells_interaction2.py. tab delimited. 
#        Columns: ligand, ligand_lfc, receptor, receptor_Expression (log2(tpm+1)), cell type ligand, cell type receptor
#
# Output: 2 Tab delimited text files listing cell-cell pairs- one considering up and down regulated and one not - include total interaction, total receptors and total ligands, mean ligand lfc (for all ints and for uniq ligands).
#         (considering lignad-receptor pairs)

##### Set up #####

# Capture  messages and errors to a file.
zz <- file("gutcovid_lung.out", open="a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: 6_ligand_receptor_stats.R\n")

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 4){
  stop("Wrong number of command line input parameters. Please check.")
}

library(tidyverse)
#setwd("/Users/poletti/OneDrive\ -\ Norwich\ BioScience\ Institutes/Bioinformatics/Lung_run/ciliated_cells_moderate")

# Input file
lr_ints <- read.csv(args[1], sep = "\t")

# Outfiles
out_tog <- args[2]
out_sep <- args[3]
out_lr <-  args[4]

##### Process - per cell-cell pair #####

# Add column to say if up or downregulated
lr_ints1 <- lr_ints %>% mutate(ligand_dir = ifelse(ligand_lfc >0, "upreg","downreg")) %>% unique()

## Not accounting for up and downreg sererately

# Collapse to count all ints per cell-cell #NB. the mean lfc here is the mean of the ligand lfc across all interactions - so some ligands will be counted multiple times
lr_ints_count1 <- lr_ints1 %>% group_by(ligand_cell_source, receptor_cell_target) %>% summarise(lfc_mean_all= mean(ligand_lfc), uniq_ints=n()) %>% mutate(ligand_dir = "both")

# Collapse to count receptors 
lr_ints_count2 <- lr_ints1 %>% select(-c(ligand, ligand_lfc)) %>% unique() %>% count(ligand_cell_source, receptor_cell_target) %>%
  rename("uniq_receptors"=n)

# Collapse to count ligands
lr_ints_count3 <- lr_ints1 %>% select(-c(receptor, receptor_expression)) %>% unique() %>% group_by(ligand_cell_source, receptor_cell_target) %>%
  summarise(lfc_mean= mean(ligand_lfc), uniq_ligands=n())

# Join last 3
lr_ints_tog <- list(lr_ints_count1, lr_ints_count2, lr_ints_count3) %>% reduce(full_join, by = c("ligand_cell_source", "receptor_cell_target"))

## Accounting for up and downreg sererately

# Collapse to count number of up and downreg interactions per cell-cell type
lr_ints_count4 <- lr_ints1 %>% group_by(ligand_cell_source, receptor_cell_target,ligand_dir) %>% summarise(lfc_mean_all= mean(ligand_lfc), uniq_ints=n())
                               
# Collapse to count receptors per up/down per cell-cell
lr_ints_count5 <- lr_ints1 %>% select(-c(ligand, ligand_lfc)) %>% unique() %>% count(ligand_cell_source, receptor_cell_target, ligand_dir) %>%
  rename("uniq_receptors"=n)

# Collapse to count ligand per up/down per cell-cell
lr_ints_count6 <- lr_ints1 %>% select(-c(receptor, receptor_expression)) %>% unique() %>% group_by(ligand_cell_source, receptor_cell_target, ligand_dir) %>%
  summarise(lfc_mean= mean(ligand_lfc), uniq_ligands=n())
  
# Join last 3
lr_ints_sep <- list(lr_ints_count4, lr_ints_count5, lr_ints_count6) %>% reduce(full_join, by = c("ligand_cell_source", "receptor_cell_target", "ligand_dir"))

# Save
write.table(lr_ints_tog, file= out_tog, sep = "\t", quote = F, row.names=F)
write.table(lr_ints_sep, file= out_sep, sep = "\t", quote = F, row.names=F)

# Clean
rm(lr_ints1, lr_ints_count1, lr_ints_count2, lr_ints_count3, lr_ints_count4, lr_ints_count5, lr_ints_count6)

##### Process - per ligand-receptor pair #####

# Add column to say if up or downregulated
per_lr <- lr_ints %>% mutate(ligand_dir = ifelse(ligand_lfc >0, "upreg","downreg")) %>% unique()

# Collapse by ligand-receptor pair - list the receptor cells
per_lr <- per_lr %>% select(-c(receptor_expression)) %>% 
  group_by(ligand, receptor, ligand_cell_source, ligand_dir) %>%
  summarise(target_cell_counts = n(), receptor_cell_target = paste0(receptor_cell_target, collapse=",")) 

# Save
write.table(per_lr, file= out_lr, sep = "\t", quote = F, row.names=F)

