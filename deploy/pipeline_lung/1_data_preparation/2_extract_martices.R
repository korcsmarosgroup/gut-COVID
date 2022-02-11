library("optparse")
library(Matrix)
library(data.table)
library(tidyr)
library(Seurat)
library(tools)

# Capture  messages and errors to a file.
zz <- file("gutcovid_lung.out", open="a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: 2_extract_martices.R\n")

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 2){
  stop("Wrong number of command line input parameters. Please check.")
}


load_r_data_object <- function(file_name) {
	# Load RData Object to a new name
	if (file_ext(file_name) == "RData") {
		load(file_name)
		return(get(ls()[ls() != "file_name"]))
	} else if (file_ext(file_name) == "rds") {
		return(readRDS(file_name))
	}
}


generate_average_expression <- function(seu_object, output_file_path, ident="Cluster", add.ident=NA) {
	# Get the average expression from Seurat
	Idents(seu_object) <- ident
	avg_expression <- as.data.frame(AverageExpression(seu_object, group.by=c(ident, add.ident))$RNA)
	write.table(avg_expression, file = output_file_path, sep = '\t')
}


map_meta_data <- function(seu_object, meta_data_file_path) {
	#Map metadata to the seurat object
	meta_data <- read.delim(meta_data_file_path, sep="\t")
	rownames(meta_data) <- meta_data$NAME
	seu_object_rownames <- rownames(seu_object@meta.data)
	meta_data <- subset(meta_data, rownames(meta_data) %in% seu_object_rownames)
	meta_data_col_names <- colnames(meta_data)
	seu_object <- AddMetaData(object = seu_object, metadata = meta_data, col.name = meta_data_col_names)
	return(seu_object)
}


write_dgCMatrix_csv <- function(mat, filename, col1_name = "gene", chunk_size = 1000) {
  # Create inital matrix
  print("Making matrix...")
  mat <- Matrix::t(mat)
  row_names <- colnames(mat)
  col_names <- rownames(mat)
  n_row <- length(row_names)
  n_col <- length(col_names)
  n_chunks <- floor(n_row/chunk_size)
  print(paste("Total number of rows:", n_row))
  print(paste("Total number of columns:", n_col))
  print(paste("Number of chunks:", n_chunks))
  
  # Initial chunk
  chunk <- 1
  chunk_start <- 1 + chunk_size * (chunk - 1)
  chunk_end <- chunk_size * chunk
  print(paste0("Writing rows ",chunk_start," to ", chunk_end))
  chunk_mat <- t(as.matrix(mat[,chunk_start:chunk_end]))
  chunk_df <- cbind(data.frame(col1 = row_names[chunk_start:chunk_end]),as.data.frame(chunk_mat))
  names(chunk_df)[1] <- col1_name
  data.table::fwrite(chunk_df, file = filename, append = F)
  
  # chunkation over chunks
  for(chunk in 2:n_chunks) {
    chunk_start <- 1 + chunk_size * (chunk - 1)
    chunk_end <- chunk_size * chunk
    print(paste0("Writing rows ",chunk_start," to ", chunk_end))
    chunk_mat <- t(as.matrix(mat[,chunk_start:chunk_end]))
    chunk_df <- cbind(data.frame(col1 = row_names[chunk_start:chunk_end]),as.data.frame(chunk_mat))
    data.table::fwrite(chunk_df, file = filename, append = T)
  }
  
  # Remaining samples
  chunk_start <- (n_chunks*chunk_size + 1)
  chunk_end <- n_row
  print(paste0("Writing rows ",chunk_start," to ", chunk_end))
  chunk_mat <- t(as.matrix(mat[,chunk_start:chunk_end]))
  chunk_df <- cbind(data.frame(col1 = row_names[chunk_start:chunk_end]),as.data.frame(chunk_mat))
  data.table::fwrite(chunk_df, file = filename, append = T) 
}


extract_genes <- function(mat, gene_names, selected_gene_names_output_file_path, chunk_size=100) {
	# Extract a list of specific genes and transverse the matrix
	if (!(genes_names %in% rownames(mat))) {
		print("One or more gene names not found.")
	} else {
		print("All gene names found.")
		selected_gene_mat <- selected_gene_counts(t(mat[gene_names,]))
		write_dgCMatrix_csv(mat=selected_gene_mat, 
		                    output_file_path=selected_gene_names_output_file_path, 
		                    chunk_size=chunk_size)	
	}
	
}
 

run_smilie_et_al_extraction <- function(seu_object_path, meta_data_file_path, data_label, chunk_size=100, ident="Cluster", add_ident="Health") {
	# Main logic for the Smilie et al final data processing

	# Load Seurat Object
	print("Loading data...")
	seur <- load_r_data_object(file_name = seu_object_path)

	# Add meta data
	print("Adding metadata...")
	if(!is.na(meta_data_file_path)) {
		seur <- map_meta_data(seu_object = seur, meta_data_file_path = meta_data_file_path)
	}

	# Write out normalised counts
	print("Extracting and writing normalised counts...")
	mat <- GetAssayData(object = seur, slot = "data")
	normalised_counts_filename <- paste0(data_label, "_normalised_counts.csv")
	write_dgCMatrix_csv(mat=mat, filename=normalised_counts_filename, chunk_size = chunk_size)
	rm(mat)
	gc()

	print("Extracting and writing counts...")
	mat <- GetAssayData(object = seur, slot = "counts")
	counts_filename <- paste0(data_label, "_counts.csv")
	write_dgCMatrix_csv(mat=mat, filename=counts_filename, chunk_size = chunk_size)
	rm(mat)
	gc()

	# Generate average expression
	print("Generating average expression...")
	avg_expression_file_path <- paste0(data_label, "_average_epxression.tsv")
	generate_average_expression(seu_object=seur, output_file_path=avg_expression_file_path, ident=ident, add.ident=add_ident)

	print("Done.")
}


if (!interactive()) {
	#Get commandline inputs
	option_list <- list(make_option(c("-s", "--seurat"), action="store", default=NA, type='character', help="File path to seurat object."),
					    make_option(c("-m", "--metadata"), action="store", default=NA, type='character', help="File path to meta data."),
		                make_option(c("-d", "--datalabel"), action="store", default=NA, type='character', help="Name of the dataset."),
		                make_option(c("-c", "--chunksize"), action="store", default=100, type='integer', help="Chunk size based of available RAM."),
		                make_option(c("-i", "--ident"), action="store", default=NA, type='character', help="Main ident (typically the cluster)."),
		                make_option(c("-a", "--addident"), action="store", default=NA, type='character', help="Additional ident."))
	argv <- parse_args(OptionParser(option_list=option_list))
	run_smilie_et_al_extraction(seu_object_path=argv$seurat, 
								meta_data_file_path=argv$metadata, 
								data_label=argv$datalabel, 
								chunk_size=argv$chunksize,
								ident=argv$ident,
								add_ident=argv$addident)
}
run_smilie_et_al_extraction(
seu_object_path=args[1], 
meta_data_file_path=args[2], 
data_label="COVID_lung_dataset_nb_mainv2", 
chunk_size=1000,
ident="celltype",
add_ident="severity")
