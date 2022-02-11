# Script to get the number of ligands per intercellular network
#
# Author : Agatha Treveil
#
# Input: 'ligand_receptor_ints' networks created with scripts 
#         epithelial_organoid_immune_cells_interaction_colon{ileum}.py. One file for each
#         colon infected, coln bystander, ileum infected, ileum bystander
# Ouput: Table of # ligands for each intercellular network

#### Set up ####

library(dplyr)

setwd("/Users/poletti/OneDrive\ -\ Norwich\ BioScience\ Institutes/Bioinformatics/Lung_run")

# List of input files
infiles <- c("3_inter_network_reconstruction/ligand_receptor_ints_lung_epithelial_immune_cells.txt")

# Input names - same order as list of files
innames <- c("squamousous_cells")

#### Process ####

# Empty df
all_res <- data.frame()

# Iterate datasets
for (i in 1:length(infiles)){
  
  # Load the file
  whole_net <- read.csv(infiles[i], sep = "\t")
  
  # Get ligand dir
  whole_net <- whole_net %>% mutate(ligand_dir = if_else(ligand_lfc > 0, "upreg", "downreg"))
  
  # Remove unnecessary cols and group by cell types, counting number of unique ligands
  net <- whole_net %>% select(c(ligand, ligand_dir, ligand_cell_source, receptor_cell_target)) %>%
    group_by(ligand_cell_source, receptor_cell_target, ligand_dir) %>% 
    summarise(num_distict_ligands =n_distinct(ligand)) %>%
    mutate(condition = innames[i])
  
  # Merge to df
  all_res <- rbind(all_res, net)
  
}

# Save data
write.table(all_res, "number_ligands_per_network_up_down.txt", sep = "\t", quote=F, row.names = F)
