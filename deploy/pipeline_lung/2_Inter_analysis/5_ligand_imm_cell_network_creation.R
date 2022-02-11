# Script to combine the stat table to create ligand -> imm cell network table
#
# Input: lignad_rec table output from epithelial_organoid_immune_cells_interaction_{colon/ileum}.py
#        cell type of interest for ligand source
#
# Output: Tab delimited table of ligand-> immcell interactions

##### Set up #####

library(dplyr)

# Capture  messages and errors to a file.
zz <- file("gutcovid_lung.out", open="a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: 5_ligand_imm_cell_network_creation.R\n")

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 2){
  stop("Wrong number of command line input parameters. Please check.")
}

#setwd("/Users/poletti/OneDrive\ -\ Norwich\ BioScience\ Institutes/Bioinformatics/Lung_run/ciliated_cells_moderate")

stats <- read.csv(args[1], sep = "\t")
ligand_cell <- "ciliated_cells"

##### Process #####

# filter for ligand cell source
stats <- stats %>% filter(ligand_cell_source == ligand_cell)

stats2 <- stats %>% select(c(ligand, ligand_lfc, receptor, receptor_expression, receptor_cell_target)) %>%
  group_by(ligand, ligand_lfc, receptor_cell_target) %>% 
  summarise(receptors = paste0(receptor, collapse = ","), n_receptors = n(), mean_receptor_exp = mean(receptor_expression), sum_receptor_exp = sum(receptor_expression))

# Save
write.table(stats2, file = args[2], sep = "\t", quote = F, row.names = F)
