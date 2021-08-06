# Script to combine the stat table to create ligand -> imm cell network table
#
# Input: lignad_rec table output from epithelial_organoid_immune_cells_interaction_{colon/ileum}.py
#        cell type of interest for ligand source
#
# Output: Tab delimited table of ligand-> immcell interactions

##### Set up #####

library(dplyr)

setwd("/Volumes/group-gc/Projects/CovidOmix/stream4_gutcovid/intercellular_networks/results/triana_2020/ligand_receptor_ints/3rd_set_imm_cells/colon/infected/network_figs/ligand_receptor_immcell/imm_enterocytes_2/")

stats <- read.csv("../../../ligand_receptor_ints_triana_colon_infec_smillie_colon_health.txt", sep = "\t")
ligand_cell <- "imm_enterocyte_2"

##### Process #####

# filter for ligand cell source
stats <- stats %>% filter(ligand_cell_source == ligand_cell)

stats2 <- stats %>% select(c(ligand, ligand_lfc, receptor, receptor_expression, receptor_cell_target)) %>%
  group_by(ligand, ligand_lfc, receptor_cell_target) %>% 
  summarise(receptors = paste0(receptor, collapse = ","), n_receptors = n(), mean_receptor_exp = mean(receptor_expression), sum_receptor_exp = sum(receptor_expression))

# Save
write.table(stats2, file = "colon_infec_imm_entero2_ligand_imm_cell_network.txt", sep = "\t", quote = F, row.names = F)
