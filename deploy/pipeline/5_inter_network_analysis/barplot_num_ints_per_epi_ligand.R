# R script to create bar plot of number of interactions across all immune cells for each epithelial deg/ligand
#
# Author: Agatha Treveil, December 2020
#
# Input: All interactions output from the python scripts epithelial_organoid_immune_cells_interaction_ileum{colon}.py
#
# Output: Png plot

# Capture  messages and errors to a file.
zz <- file("gutcovid.out", open = "a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: barplot_num_ints_per_epi_ligand.R\n")

##### Set up #####

if (!requireNamespace("hrbrthemes", quietly = TRUE)) 
  install.packages("hrbrthemes")

library(tidyverse)
library(ggplot2)
library(hrbrthemes)
library(tidytext)
library(grid)
#setwd("/Volumes/group-gc/Projects/CovidOmix/stream4_gutcovid/intercellular_networks/results")

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 2){
  stop("Wrong number of command line input parameters. Please check.")
}

# Input files
in_file <- read.csv(args[1], sep = "\t")

# Outfile
plot_out <- args[2]

##### Process degs #####

# Up or downregulated ligand
ints <- in_file %>% mutate(ligand_dir = ifelse(ligand_lfc > 0, "upreg", "downreg"))

# Collapse per ligand-celltype, count num interactions, count num of imm cell types, keep ligand lfc
ints2 <- ints %>% group_by(ligand_cell_source, ligand, ligand_dir) %>%
  summarise(num_imm_cells = n_distinct(receptor_cell_target), num_ints = n()) %>%
  #arrange(ligand_cell_source, num_ints)
  ungroup() %>%
  mutate(ligand_cell_source = as.factor(ligand_cell_source),
         ligand = reorder_within(ligand, num_ints, ligand_cell_source))

# For cell type specific plots
ints3 <- ints %>% group_by(ligand_cell_source, ligand, ligand_dir) %>%
  summarise(num_imm_cells = n_distinct(receptor_cell_target), num_ints = n()) 

##### Plot facet #####

# Values used to transform data
coeff <- 10

# Plot data
plot_ligs <- ggplot(ints2, aes(x=ligand, fill = ligand_dir)) +
  geom_bar(aes(y=num_ints), stat="identity") +
  scale_fill_manual("legend", values = c("upreg" = "#F8766D", "downreg" = "#00BFC4"))+
  geom_point(aes(y=num_imm_cells*coeff), size=2, colour ="black") +
  #facet_grid(~ligand_cell_source, scales = "free_x", space='free') + # Good for the infected cells where some have way more ligands than others
  facet_wrap(~ligand_cell_source, scales = "free_x") + # Good for the bystanders where plots similar size between cell types
  scale_x_reordered() +
  scale_y_continuous(
    # Features of the first axis
    name = "Number of intercellular interactions",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./coeff, name="Number of immune cells targeted")) +
  theme_light()+
  theme(
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black")) +
  theme(legend.position = "bottom") +
  theme(
    strip.text.x = element_text(
      size = 5, color = "black", face = "bold"),
    strip.text.y = element_text(
      size = 5, color = "black", face = "bold")) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(size = 6))    

plot_ligs

# Save
ggsave(plot_ligs, file=file.path(paste0(plot_out,"num_ints_per_ligand.png")),width = 9, height=5)


##### Plot per epi cell type #####

# Get unique epithelial cell types
cell_types <- unique(ints3$ligand_cell_source)

for (cell in cell_types) {
  # Filter data for that epi cell type
  ints3_cell <- ints2 %>% filter(ligand_cell_source == cell)
  
  # Values used to transform data
  coeff <- 10
  
  # Plot data
  plot_ligs2 <- ggplot(ints3_cell, aes(x=ligand, fill = ligand_dir)) +
    scale_fill_manual("legend", values = c("upreg" = "#F8766D", "downreg" = "#00BFC4"))+
    geom_bar(aes(y=num_ints), stat="identity") +
    geom_point(aes(y=num_imm_cells*coeff), size=2, colour ="black") +
    scale_x_reordered() +
    scale_y_continuous(
      # Features of the first axis
      name = "Number of intercellular interactions",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~./coeff, name="Number of immune cells targeted")) +
    theme_light()+
    theme(
      axis.title.y = element_text(color = "black"),
      axis.title.y.right = element_text(color = "black")) +
    theme(legend.position = "bottom") +
    theme(
      strip.text.x = element_text(
        size = 5, color = "black", face = "bold"),
      strip.text.y = element_text(
        size = 5, color = "black", face = "bold")) +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          text = element_text(size = 6))    
  
  # Save
  ggsave(plot_ligs2, file=file.path(paste0(plot_out, "num_ints_per_ligand_",cell, ".png")),width = 5, height=2.5)
  
}
