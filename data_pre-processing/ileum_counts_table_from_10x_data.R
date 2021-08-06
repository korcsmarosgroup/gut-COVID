# Get counts table from processed 10x data. Get DEGs for each cell type, sars infected vs mock.

# Input: Preprocessed files from GEO (GSE156760) - Rdata file
# Output: 2 tab delimited tables of average gene expression across each cell type/condition (using RNA assay): all cells and only infected cells.
#         DEG tables, filtered (p adj <= 0.05 and lfc 0.5). all cells vs mock, infected vs mock and bystander vs mock.

#### Set up ####

# Packages
library(Seurat)
library(MAST)
library(dplyr)

setwd("/Users/treveila/Documents/covid/updated_on_groupgc/stream4/")

load("data/triana_et_al_2020/COVID19_Oct.rda")

outdir <- "results/triana_et_al_2020/test/"

#### Get counts table all cell s#####

# Create copy
ileum <- Illeum_H_T

# Set identities as condition_celltype
Idents(object=ileum) <- paste0(Idents(object=ileum), "_", ileum$CellTypes)

# Get average expression per cluster for each gene (in the RNA assay)
ave_exp <- AverageExpression(ileum, assays="RNA")

# genes as column not rownames
ave_exp$RNA$gene <- row.names(ave_exp$RNA)
ave_exp$RNA <- ave_exp$RNA %>% select(gene, everything())
rownames(ave_exp$RNA) <- NULL

# Save
write.table(ave_exp$RNA, file = file.path(outdir,"expression","ileum_average_expression_per_celltype.txt"), sep = "\t", quote = F, row.names = F)


##### Get counts table infected cells vs uninfected cells #####

# Create copy
ileum_inf <- Illeum_H_T

# Set identities as condition_celltype
Idents(object=ileum_inf) <- paste0(Idents(object=ileum_inf), "_", ileum$CellTypes, "_", ileum$State)

# Get average expression per cluster for each gene (in the RNA assay)
ave_exp_inf <- AverageExpression(ileum_inf, assays="RNA")

# genes as column not rownames
ave_exp_inf$RNA$gene <- row.names(ave_exp_inf$RNA)
ave_exp_inf$RNA <- ave_exp_inf$RNA %>% select(gene, everything())
rownames(ave_exp_inf$RNA) <- NULL

# Save
write.table(ave_exp_inf$RNA, file = file.path(outdir, "expression","ileum_infected_bystander_cells_average_expression_per_celltype.txt"), sep = "\t", quote = F, row.names = F)


##### Get DEGs all cells #####

get_degs <- function(data, ident1a, ident1b, ident2a, ident2b, out1, out2){
  
  # Empty table fo DEG results
  all_degs <- data.frame()
  
  # Calculate DEGs between ileum 24 and mock for each cell type
  for (i in unique(ileum$CellTypes)){
    #print(i)
    
    identity1 <- paste0(ident1a, i, ident1b)
    identity2 <- paste0(ident2a, i, ident2b)
    
    tryCatch({
      # Get DEGs for all cell 
      degs <- FindMarkers(data, ident.1 = identity1, ident.2 = identity2, assay='RNA',test.use = "MAST")
      
      # Put gene names as column not row names
      degs$gene <- rownames(degs)
      
      # Add comparison as column
      degs <- degs %>% mutate(comparison = paste0(paste0(identity1, '_v_',identity2)))
      
      # Append to previous results
      all_degs <- rbind(all_degs, degs)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  # Filter for p adj <= 0.05 and lfc >= 1
  all_degs_f <- all_degs %>% dplyr::filter((p_val_adj <= 0.05)& ((avg_logFC >= 0.5)|(avg_logFC <= -0.5)))
  
  # Save
  write.table(all_degs, file = file.path(outdir, "degs", out1), sep = "\t", quote=F, row.names = F)
  write.table(all_degs_f, file = file.path(outdir, "degs", out2), sep = "\t", quote=F,row.names = F)
  
}

# DEGs between 24h and mock
get_degs(ileum, 'Illeum_24h_', "", 'Illeum_Mock_', "", "unfilt_degs_ileum_24h_v_mock_all_cells.txt", "filt_degs_ileum_24h_v_mock_all_cells_padj0.05_lfc0.5.txt")

# DEGS between 24h infected and mock
get_degs(ileum_inf, 'Illeum_24h_', "_Infected", 'Illeum_Mock_', "_Non-Infected", "unfilt_degs_ileum_24h_infected_v_mock_all_cells.txt", "filt_degs_ileum_24h_infected_v_mock_all_cells_padj0.05_lfc0.5.txt")

# DEGS between 24h bystander and mock
get_degs(ileum_inf, 'Illeum_24h_', "_Bystander", 'Illeum_Mock_', "_Non-Infected", "unfilt_degs_ileum_24h_bystander_v_mock_all_cells.txt", "filt_degs_ileum_24h_bystander_v_mock_all_cells_padj0.05_lfc0.5.txt")


##### Other useful #####

# Cell types annotation data
#ileum_H_T$CellTypes

# Timepoint of samples
#ileum_H_T$Cond

# Timepoint of samples
#ileum_H_T$State

# to get counts data use @assays$RNA@counts, But for the normalised data you would use @assays$RNA@data.

# Counts data - use RNA not SCT https://www.biostars.org/p/395951/
#ileum_H_T@assays[["RNA"]]@counts

## Get only infected 24h cells
#ileum_inf <- subset(ileum, subset = State == "Infected")
#ave_exp_inf <- AverageExpression(ileum_inf, assays="RNA")

