# Script to get the number of ligands per intercellular network
#
# Author : Agatha Treveil
#
# Input: 'ligand_receptor_ints' networks created with scripts 
#         epithelial_organoid_immune_cells_interaction_colon{ileum}.py. One file for each
#         colon infected, coln bystander, ileum infected, ileum bystander
# Ouput: Table of # ligands for each intercellular network

# Capture  messages and errors to a file.
zz <- file("gutcovid.out", open = "a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: get_num_ligands_per_network.R\n")

#### Set up ####

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 5){
  stop("Wrong number of command line input parameters. Please check.")
}

library(dplyr)

#setwd("/Volumes/group-gc/Projects/CovidOmix/stream4_gutcovid/intercellular_networks/results/triana_2020/ligand_receptor_ints/3rd_set_imm_cells")

# List of input files
infiles <- c(args[1],
             args[2],
             args[3],
             args[4])

# Input names - same order as list of files
innames <- c("infected_colon","bystander_colon","infected_ileum","bystander_ileum")

#### Process ####

# Empty df
all_res <- data.frame()

# Iterate datasets
for (i in 1:length(infiles)){
  
  # Load the file
  whole_net <- read.csv(infiles[i], sep = "\t")
  
  # Remove unnecessary cols and group by cell types, counting number of unique ligands
  net <- whole_net %>% select(c(ligand, ligand_cell_source, receptor_cell_target)) %>%
    group_by(ligand_cell_source, receptor_cell_target) %>% 
    summarise(num_distict_ligands =n_distinct(ligand)) %>%
    mutate(condition = innames[i])
  
  # Merge to df
  all_res <- rbind(all_res, net)
  
}

# Save data
write.table(all_res, args[5], sep = "\t", quote=F, row.names = F)
