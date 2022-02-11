# R script to create bar plot of number of DEGs for epithelial cells
#
# Author: Agatha Treveil, December 2020
#
# Input: Degs that are ligands in the epithelial-immune cell intercellular networks (generated using get_num_ligands_per_network.R script)
#
# Output: png bar plot numbe rof degs across all 4 condition (colon, ileum, bystander, infected)

##### Set up #####

library(tidyverse)
library(ggplot2)
library(tidytext)
setwd("/Users/poletti/OneDrive\ -\ Norwich\ BioScience\ Institutes/Bioinformatics/Lung_run/ciliated_cells")

# Input files
ligands <- read.csv("number_ligands_per_network_up_down.txt", sep ="\t")

ligands_by_group <- ligands %>%                                       
  group_by(ligand_cell_source, receptor_cell_target, ligand_dir, condition) %>%    
  summarise_at(vars(num_distict_ligands),              
               list(num_ligands_condition = sum))               

##### Plot #####

#fills <- rev(brewer.pal(2, "Set2"))

# Facet bar plot  
out_p <- ggplot(data = ligands_by_group, aes(x=receptor_cell_target, y=num_ligands_condition, fill=ligand_dir, label = num_ligands_condition)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~condition, scales = "free_y") +
  scale_x_reordered() + 
  geom_col(show.legend = FALSE) + 
  theme_bw() +
  scale_fill_manual(values=c("#04bcc4", "#fc746c")) +
  coord_flip() +
  ylab("Number of ligands in intercellular network") +
  xlab("Epithelial cell type") 

out_p

# Save
ggsave(file = "number_ligands_per_network_up_down_sep.png", width = 10, height=5)
