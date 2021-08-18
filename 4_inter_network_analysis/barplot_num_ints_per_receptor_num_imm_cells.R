# R script to create bar plot of number of interactions across all epithelial cells for each immune cell receptor
# Not split up by immune cell
#
# Author: Agatha Treveil, December 2020
#
# Input: All interactions output from the python scripts epithelial_organoid_immune_cells_interaction_ileum{colon}.py
#
# Output: Png plots - per epithelial cell type. facter plot across epi cell types andall epi cell types together

##### Set up #####

library(tidyverse)
library(ggplot2)
library(hrbrthemes)
library(tidytext)
setwd("/Volumes/group-gc/Projects/CovidOmix/stream4_gutcovid/intercellular_networks/results/triana_2020/ligand_receptor_ints/3rd_set_imm_cells")

# Input files
#in_file <- read.csv("ileum/infected/ligand_receptor_ints_triana_ileum_infec_martin_si_uninfl.txt", sep = "\t")
#in_file <- read.csv("ileum/bystander/ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt", sep = "\t")
#in_file <- read.csv("colon/bystander/ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt", sep = "\t")
in_file <- read.csv("colon/infected/ligand_receptor_ints_triana_colon_infec_smillie_colon_health.txt", sep = "\t")

# Outfile dir
#plot_out <- "ileum/infected/plots/receptors/num_imm_cells/top25_receptors"
#plot_out <- "ileum/bystander/plots/receptors/num_imm_cells"
#plot_out <- "colon/bystander/plots/receptors/num_imm_cells"
plot_out <- "colon/infected/plots/receptors/num_imm_cells/top25_receptors"

##### Plot non facet function #####

plot_non_facet <- function(in_data, outfile, wid, hei, coeff, axis2){
  # Values used to transform data
  
  if (axis2 == TRUE){
    plot_ligs_tog <- ggplot(in_data, aes(x=reorder(receptor, num_ints), num_ints, fill = n_ligands)) +
      geom_bar(aes(y=num_ints), stat="identity") +
      scale_fill_continuous(low="lightblue", high="blue2",limits=c(0,7)) +
      geom_point(aes(y=num_imm_cells*coeff), size=2, colour ="black") +
      scale_x_reordered() +
      scale_y_continuous(
        # Features of the first axis
        name = "Number of intercellular interactions",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~./coeff, name="Number of immune cells targeted")) +
      theme_light()+
      theme(axis.title.y = element_text(color = "black"),
            axis.title.y.right = element_text(color = "black")) +
      theme(legend.position = "bottom") +
      theme(strip.text.x = element_text(size = 8, color = "black", face = "bold"),
            strip.text.y = element_text(size = 8, color = "black", face = "bold")) +
      theme(axis.text.x = element_text(angle=45, hjust=1),text = element_text(size = 8)) +
      xlab("Receptor")
  } else {
    plot_ligs_tog <- ggplot(in_data, aes(x=reorder(receptor, num_ints), num_ints, fill = n_ligands)) +
      geom_bar(aes(y=num_ints), stat="identity") +
      scale_fill_continuous(low="lightblue", high="blue2", limits=c(0,7)) +
      scale_x_reordered() +
      scale_y_continuous(name = "Number of intercellular interactions")+
      theme_light()+
      theme(axis.title.y = element_text(color = "black")) +
      theme(legend.position = "bottom") +
      theme(strip.text.x = element_text(size = 8, color = "black", face = "bold"),
            strip.text.y = element_text(size = 8, color = "black", face = "bold")) +
      theme(axis.text.x = element_text(angle=45, hjust=1),text = element_text(size = 8)) +
      xlab("Receptor")
    
  }
  
  # Save
  ggsave(plot_ligs_tog, file=outfile, width=wid, height=hei)
}

##### Process degs #####

# Split up and down regulated interactions
in_file <- in_file %>% mutate(ligand_dir = ifelse(ligand_lfc >0, "upreg","downreg"))
ints_spl <- split(in_file, in_file$ligand_dir)

for (name in names(ints_spl)){
  
  # Get interactions with that direction
  ints <- ints_spl[[name]]
  
  # Collapse per receptor, count num interactions, count num of imm cell types
  ints_sep <- ints %>% select(ligand_cell_source, receptor, receptor_cell_target) %>% group_by(ligand_cell_source, receptor) %>%
    summarise(num_imm_cells = n_distinct(receptor_cell_target), num_ints = n()) %>%
    #arrange(receptor_cell_target, num_ints)
    ungroup() %>%
    mutate(ligand_cell_source = as.factor(ligand_cell_source),
           receptor = reorder_within(receptor, num_ints, ligand_cell_source))
  
  ints_tog <- ints %>% group_by(receptor) %>%
    summarise(num_imm_cells = n_distinct(receptor_cell_target), num_ints = n(), n_ligands = n_distinct(ligand))
  
  #ints_ie2 <- ints %>% filter((ligand_cell_source == "imm_enterocyte_2")|(ligand_cell_source == "imm_enterocytes_2")) %>% group_by(receptor) %>%
  #  summarise(num_imm_cells = n_distinct(receptor_cell_target), num_ints = n(), n_ligands = n_distinct(ligand))
    
  ##### Plot facet #####
  
  
  # Plot data
  plot_ligs <- ggplot(ints_sep, aes(x=receptor, fill = num_imm_cells)) +
    geom_bar(aes(y=num_ints), stat="identity") +
    scale_fill_continuous(low="lightblue", high="darkblue") +
    #facet_grid(~ligand_cell_source, scales = "free_x", space='free') + # Good for the infected cells where some have way more ligands than others
    facet_wrap(~ligand_cell_source, scales = "free_x") + # Good for the bystanders where plots similar size between cell types
    scale_x_reordered() +
    scale_y_continuous(
      # Features of the first axis
      name = "Number of intercellular interactions")+
    theme_light()+
    theme(
      axis.title.y = element_text(color = "black"))+
    theme(legend.position = "bottom") +
    theme(
      strip.text.x = element_text(
        size = 5, color = "black", face = "bold"))+
    theme(axis.text.x = element_text(angle=45, hjust=1),
          text = element_text(size = 6))    
  
  #plot_ligs
  
  # Filename
  file_n <- paste0("facet_num_ints_per_receptor_num_imm_cells_", name, ".png")
  
  # Save
  ggsave(plot_ligs, file=file.path(plot_out,file_n),width = 30, height=8) 
  
  ##### Plot together #####
  
  file_n_sep <- paste0("together_num_ints_per_receptor_num_imm_cells_", name, ".png")
  
  plot_non_facet(ints_tog, file.path(plot_out, file_n_sep),10, 3, 5, TRUE)
  
  
  ##### Plot facet per epi cell type #####   #this is the one I need
  
  # Get unique epithelial cell types
  cell_types <- unique(ints_sep$ligand_cell_source)

  for (cell in cell_types) {
    # Filter data for that epi cell type
    ints_sep_cell <- ints %>% filter(ligand_cell_source == cell) %>% group_by(receptor) %>%
      summarise(num_imm_cells = n_distinct(receptor_cell_target), num_ints = n(), n_ligands = n_distinct(ligand)) #%>% top_n(25, num_ints) 
    # Get file path
    file_p <- paste0(cell, "_num_ints_per_receptor_num_imm_cells_", name, ".png")
    # plot
    plot_non_facet(ints_sep_cell, file.path(plot_out, file_p),8, 3, 4, TRUE)
  }
  
  
}
  
