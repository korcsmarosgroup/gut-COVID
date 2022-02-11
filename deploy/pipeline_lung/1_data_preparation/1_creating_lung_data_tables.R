#Loading requeired packages
library(dplyr)
library(Seurat)
library(patchwork)
require(data.table)
library(ggplot2)
#Loading the SEURAT obejct 
#setwd("/Users/poletti/OneDrive\ -\ Norwich\ BioScience\ Institutes/Bioinformatics/Lung_run")

# Capture  messages and errors to a file.
zz <- file("gutcovid_lung.out", open="a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: 1_creating_lung_data_tables.R\n")

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 5){
  stop("Wrong number of command line input parameters. Please check.")
}

#Quick checks and log normalisation
loc <- readRDS(args[1])
VlnPlot(loc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
loc <- NormalizeData(loc, normalization.method = "LogNormalize", scale.factor = 10000)
loc <- FindVariableFeatures(loc, selection.method = "vst", nfeatures = 2000)


## Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(loc), 10)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

## plot variable features with and without labels
#plot1 <- VariableFeaturePlot(loc)
#plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

#A quick umap plot creation, because it is fun to do and beutiful
#loc <- ScaleData(loc)
#loc <- RunPCA(loc, npcs =50)
#loc <- RunUMAP(loc, dims = 1:50)
#DimPlot(object = loc, reduction = 'umap')

#Cheking the data in the amtrix, hether MAtthew's script will be able toa cces it.
mat <- GetAssayData(object = loc, slot = "data")

#We will need DEGS regarding the COVID vs healthy in the lung laveolar epithelial cells.
# In the lung therea re two type of aleveolar cells Type I - aka squamousus aleveolar cells
# and Type II alveolar cells aka secretory cells. THe cilliated cells are in the upper 
# airway all the rest are various immune cells if I udnerstand it correctly. 
#Iwill use the squamous and secreatory alveiolar cells. 
meta_data <- loc@meta.data
write.table(meta_data, args[2], quote = FALSE)

meta_data$covid_cell_type <- paste(meta_data$infection,meta_data$celltype)
meta_data$severity_cell_type <- paste(meta_data$severity,meta_data$celltype)
loc@meta.data <- meta_data

unique(meta_data$covid_cell_type)
unique(meta_data$severity_cell_type)
Idents(object = loc) <- loc@meta.data$'covid_cell_type'
Idents(object = loc) <- loc@meta.data$'severity_cell_type'

levels(loc)


#Finding differetially expressed gens
deg_lung_squamousous_alveolar_cells <- FindMarkers(loc, ident.1 = "critical Squamous",
                                                   ident.2 = "control Squamous")
deg_lung_secreatory_alveolar_cells <- FindMarkers(loc, ident.1 = "critical Secretory", 
                                                  ident.2 = "control Secretory")
deg_lung_cilliary_alveolar_cells <- FindMarkers(loc, ident.1 = "critical Ciliated", 
                                                  ident.2 = "control Ciliated")

write.table(deg_lung_squamousous_alveolar_cells, args[3], 
            sep="\t", quote = FALSE)
write.table(deg_lung_secreatory_alveolar_cells, args[4], 
            sep="\t", quote = FALSE)
write.table(deg_lung_secreatory_alveolar_cells, args[5], 
            sep="\t", quote = FALSE)
saveRDS(loc, file = args[6])
