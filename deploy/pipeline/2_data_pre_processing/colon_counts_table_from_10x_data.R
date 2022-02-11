# Get counts table from processed 10x data. Get DEGs for each cell type, sars infected vs mock.

# Input: Preprocessed files from GEO (GSE156760) - Rdata file
# Output: 2 tab delimited tables of average gene expression across each cell type/condition (using RNA assay): all cells and only infected cells.
#         DEG tables, filtered (p adj <= 0.05 and lfc 0.5). all cells vs mock, infected vs mock and bystander vs mock.

#### Set up ####

# Capture  messages and errors to a file.
zz <- file("gutcovid.out", open="a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: colon_counts_table_from_10x_data.R\n")

# Packages
library(Seurat)
library(MAST)
library(dplyr)

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 1){
  stop("Wrong number of command line input parameters. Please check.")
}

input_file <- args[1]

load(input_file)

outdir <- "/home/gutcovid/pipeline/6_inter_output_data"

generate_average_expression <- function(seu_object, 
                                        output_file_path, 
                                        ident="Cluster",
                                        chunk_size=1000) {
  # Get the average expression from Seurat
  Idents(seu_object) <- ident
  avg_exp <- AverageExpression(seu_object, verbose = FALSE)$RNA
  avg_exp <- cbind(gene = rownames(avg_exp), avg_exp)
  l = split(1:nrow(avg_exp), ceiling(seq_along(1:nrow(avg_exp))/chunk_size))
  for(i in l){
    write.table(avg_exp[i,],
                file=output_file_path, 
                sep = "\t",
                row.names = FALSE,
                quote = FALSE,
                append = TRUE,
                col.names=!file.exists(output_file_path))
  }
  gc()
}

#### Get counts table all cell s#####

ident <- paste0(Idents(object=Colon_H_T), "_", Colon_H_T$CellTypes)
generate_average_expression(seu_object = Colon_H_T,
                            output_file_path = file.path(outdir, "expression", "colon_average_expression_per_celltype.txt"),
                            ident=ident)

##### Get counts table infected cells vs uninfected cells #####

ident <- paste0(Idents(object=Colon_H_T), "_", Colon_H_T$CellTypes, "_", Colon_H_T$State)
generate_average_expression(seu_object = Colon_H_T,
                            output_file_path = file.path(outdir, "expression", "colon_infected_bystander_cells_average_expression_per_celltype.txt"),
                            ident=ident)

##### Get DEGs all cells #####

get_degs <- function(data, ident1a, ident1b, ident2a, ident2b, out1, out2){
  
  # Empty table fo DEG results
  all_degs <- data.frame()
  
  # Calculate DEGs between colon 24 and mock for each cell type
  for (i in unique(Colon_H_T$CellTypes)){
    #print(i)
    
    identity1 <- paste0(ident1a, i, ident1b)
    identity2 <- paste0(ident2a, i, ident2b)
    
    # Get DEGs for all cell 
    degs <- FindMarkers(data, ident.1 = identity1, ident.2 = identity2, assay = 'RNA',test.use = "MAST")
    
    # Put gene names as column not row names
    degs$gene <- rownames(degs)
    
    # Add comparison as column
    degs <- degs %>% mutate(comparison = paste0(paste0(identity1, '_v_',identity2)))
    
    # Append to previous results
    all_degs <- rbind(all_degs, degs)
  }
  
  # Filter for p adj <= 0.05 and lfc >= 1
  all_degs_f <- all_degs %>% dplyr::filter((p_val_adj <= 0.05)& ((avg_log2FC >= 0.5)|(avg_log2FC <= -0.5)))
  
  # Save
  write.table(all_degs, file = file.path(outdir, "degs", out1), sep = "\t", quote = F, row.names = F)
  write.table(all_degs_f, file = file.path(outdir, "degs", out2), sep = "\t", quote = F,row.names = F)
  
}

# DEGs between 24h and mock
Idents(object=Colon_H_T) <- Colon_H_T$orig.ident
Idents(object=Colon_H_T) <- paste0(Idents(object=Colon_H_T), "_", Colon_H_T$CellTypes)
get_degs(Colon_H_T, 'Colon_24h_', "", 'Colon_Mock_', "", "unfilt_degs_colon_24h_v_mock_all_cells.txt", "filt_degs_colon_24h_v_mock_all_cells_padj0.05_lfc0.5.txt")

# DEGS between 24h infected and mock
Idents(object=Colon_H_T) <- Colon_H_T$orig.ident
Idents(object=Colon_H_T) <- paste0(Idents(object=Colon_H_T), "_", Colon_H_T$CellTypes, "_", Colon_H_T$State)
get_degs(Colon_H_T, 'Colon_24h_', "_Infected", 'Colon_Mock_', "_Non-Infected", "unfilt_degs_colon_24h_infected_v_mock_all_cells.txt", "filt_degs_colon_24h_infected_v_mock_all_cells_padj0.05_lfc0.5.txt")

# DEGS between 24h bystander and mock
get_degs(Colon_H_T, 'Colon_24h_', "_Bystander", 'Colon_Mock_', "_Non-Infected", "unfilt_degs_colon_24h_bystander_v_mock_all_cells.txt", "filt_degs_colon_24h_bystander_v_mock_all_cells_padj0.05_lfc0.5.txt")
