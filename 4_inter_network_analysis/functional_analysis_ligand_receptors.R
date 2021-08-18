# R script to carry out functional analysis on lignad-receptor interactions between epitehlial cells and immune cells
# Author: Agatha Treveil, December 2020
#
# Input: ligand-receptor file output from epithelial_organoid_immune_cells_interaction2.py. tab delimited. 
#        Columns: ligand, ligand_lfc, receptor, receptor_Expression (log2(tpm+1)), cell type ligand, cell type receptor
#
# Output: 

##### Set up #####

library(tidyverse)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
setwd("/Volumes/group-gc/Projects/CovidOmix/stream4_gutcovid/intercellular_networks/results/triana_2020/ligand_receptor_ints/3rd_set_imm_cells")

# Input files
lr_ints <- read.csv("ileum/infected/ligand_receptor_ints_triana_ileum_infec_martin_si_uninfl.txt", sep = "\t")
#lr_ints <- read.csv("colon/infected/ligand_receptor_ints_triana_colon_infec_smillie_colon_health.txt", sep = "\t")
#lr_ints <- read.csv("colon/bystander/ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt", sep = "\t")
#lr_ints <- read.csv("ileum/bystander/ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt", sep = "\t")

backgr_lr_ints <- read.csv("../../../../data/ligand_receptor_interaction/lr_Network_Omnipath_23_09.tsv", sep = "\t")
  
# Out folder
#out_fold <- "ileum/bystander/functional_analysis"
out_fold <- "ileum/infected/functional_analysis"
#out_fold <- "colon/infected/functional_analysis"
#out_fold <- "colon/bystander/functional_analysis"

##### Functional analysis function #####

func_analysis <- function(in_list, backgr, label, out_dir){
  # Function to run overrepresentation analysis and save results
  # Runs GO and Reactome overrepresentation analysis ont he in_list usign backgr as the background.
  # Saves data table and dot plot to out_dir using label in filename
  
  # Check in list long enough
  if (length(in_list) >= 5){
    
    ##### Carry out normal GO gene overrepresentation analysis
    go1 <- enrichGO(gene = in_list, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", universe =backgr, ont = 'BP', qvalueCutoff = 0.05)
    
    # Plot and save GO output
    if (nrow(go1) > 0) {
     
       # Get dot plot
      dot_plot <- dotplot(go1, showCategory=10, orderBy="qvalue", font.size = 10)
      
      # Save dot plot
      filep <- file.path(out_dir,paste0("dot_go_overrep_", label, ".pdf"))
      pdf(filep, width=9)
      print(dot_plot)
      dev.off()
      
      # Save datatable
      filed <- file.path(out_dir,paste0("table_go_overrep_", label, ".txt"))
      write.table(as.data.frame(go1), file = filed, sep = "\t", quote = F, row.names = F)
    }
    
    ##### Carry out Reactome analysis
    
    # Convert to entrez genes
    in_list_e <- bitr(in_list, fromType='SYMBOL', toType='ENTREZID', OrgDb="org.Hs.eg.db")
    backgr_e <- bitr(backgr, fromType='SYMBOL', toType='ENTREZID', OrgDb="org.Hs.eg.db")
    
    # Carry out normal reactome gene overrepresentation analysis (can only use entrez :( )
    re1 <- enrichPathway(gene = in_list_e$ENTREZ, organism = "human", universe =as.character(backgr_e$ENTREZ), qvalueCutoff = 0.05)
    
    # Plot and save Reactome output
    if (nrow(re1) > 0) {
      
      # Get dot plot
      dot_plotr <- dotplot(re1, showCategory=10, orderBy="qvalue", font.size = 10)
      
      # Save dot plot
      filerp <- file.path(out_dir,paste0("dot_reactome_overrep_", label, ".pdf"))
      pdf(filerp)
      print(dot_plotr)
      dev.off()
      
      # Save datatable
      filerd <- file.path(out_dir,paste0("table_reactome_overrep_", label, ".txt"))
      write.table(as.data.frame(re1), file = filerd, sep = "\t", quote = F, row.names = F)
    }
  
  }
  
  
}
  
##### Process background lr interactions #####

backgr_lr_ints2 <- backgr_lr_ints %>% dplyr::select(to, from) %>% unique()

# Get lists of ligands, receptors and both
backgr_ligs <- unique(unlist(backgr_lr_ints2$from))
backgr_recs <- unique(unlist(backgr_lr_ints2$to))
backgr_both <- unique(unlist(backgr_lr_ints2))
  
##### Process + run lr interactions #####

# Split up and downregulated
lr_ints <- lr_ints %>% mutate(ligand_dir = ifelse(ligand_lfc >0, "upreg","downreg")) %>% unique()
lr_ints_spl <- split(lr_ints, f = lr_ints$ligand_dir)

# Iterate up and down regulated data
for (name in names(lr_ints_spl)){
  ints <- lr_ints_spl[[name]]
  
  # if no interactions then skip
  if (nrow(ints) == 0){
    next
  }
  
  # Empty vectors for ligand and receptors
  all_ligs <- c()
  all_recs <- c()
  
  # Split by epithelial cell type
  ints_cell <- split(ints, f=ints$ligand_cell_source)
  
  # Iterate cell types
  for (cell_name in names(ints_cell)){
    
    # Get ints for that cell
    ints_cellty <- ints_cell[[cell_name]]
    
    # if no interactions then skip
    if (nrow(ints_cellty) == 0){
      next
    }
    
    # Receptors and ligands together
    ints_tog <- ints_cellty %>% dplyr::select(ligand, receptor)
    rec_ligs <- unique(unlist(ints_tog))
    # Call functional analysis
    func_analysis(rec_ligs, backgr_both, paste0("recs+ligs_", cell_name, "_",name), file.path(out_fold,"ligands_and_receptors"))
    
    # Ligands only
    ligs <- unique(unlist(ints_tog$ligand))
    all_ligs <- c(all_ligs, ligs)
    # Call functional analysis
    func_analysis(ligs, backgr_ligs, paste0("ligs_", cell_name, "_",name), file.path(out_fold,"ligands"))
    
    # Receptors only
    recs <- unique(unlist(ints_tog$receptor))
    all_recs <- c(all_recs, recs)
    # Call functional analysis
    func_analysis(recs, backgr_recs, paste0("recs_", cell_name, "_",name), file.path(out_fold,"receptors"))
  }
  
  # Get all ligs and recs across cell types
  all_ligs <- unique(all_ligs)
  all_recs <- unique(all_recs)
  all_recs_ligs <- unique(c(all_ligs, all_recs))
  
  # Cell functional analysis
  func_analysis(all_recs_ligs, backgr_both, paste0("recs+ligs_all_",name), file.path(out_fold,"ligands_and_receptors"))
  func_analysis(all_ligs, backgr_ligs, paste0("ligs_all_",name), file.path(out_fold,"ligands"))
  func_analysis(all_recs, backgr_recs, paste0("recs_all_",name), file.path(out_fold,"receptors"))
}
  
