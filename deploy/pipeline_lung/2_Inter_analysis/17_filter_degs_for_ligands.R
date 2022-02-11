library(dplyr)

#setwd("~/OneDrive - Norwich BioScience Institutes/Bioinformatics/Lung_run/ciliated_cells_moderate")

# Capture  messages and errors to a file.
zz <- file("gutcovid_lung.out", open="a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: 17_filter_degs_for_ligands.R\n")

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 4){
  stop("Wrong number of command line input parameters. Please check.")
}

ligands <- read.table(args[1], sep =",", head=T)
ligands <- ligands %>% rename(ligand="gene")

degs <- read.table(args[2], sep = "\t", head=T)

#filter degs for ligands

degs_filtered <- degs %>% inner_join(ligands, by="gene")
degs_filtered <- degs_filtered %>% rename(gene="Gene")

#save outputs

write.table(degs_filtered, args[3], sep="\t", quote=F, row.names=F)
write.table(degs_filtered, args[4], sep=",", quote=F, row.names=F)
