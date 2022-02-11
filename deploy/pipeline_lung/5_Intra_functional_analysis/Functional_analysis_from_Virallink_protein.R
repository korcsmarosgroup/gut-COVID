#packages
library(tidyverse)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)

#set working directory
#setwd("/Users/poletti/OneDrive\ -\ Norwich\ BioScience\ Institutes/Bioinformatics/Lung_run")

# Capture  messages and errors to a file.
zz <- file("gutcovid_lung.out", open="a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: Functional_analysis_from_Virallink.R\n")

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 5){
  stop("Wrong number of command line input parameters. Please check.")
}

# Read background files
backgr_ppi <- read.csv(args[1], sep = "")  #here I have used the regular one for CARNIVAL and contextualiaed one for Virallink

# Read nodes file
nodes <- read.csv(args[2], sep = ",") #only import the one needed
nodes <- nodes %>% dplyr::rename(node = "name") 

#mirna only
nodes_mirna <- nodes %>% filter(source_node == "viralmirna") %>% dplyr::select(node)

#shared only
nodes_shared <- nodes %>% filter(source_node == "viralmirna;viralprotein") %>% dplyr::select(node)

#protein only
nodes_protein <- nodes %>% filter(source_node == "viralprotein") %>% dplyr::select(node)

# Preprocess background nodes

preproc_backgr_ppi <- function(backgr){
  # Get all nodes in ppi background network
  
  source <- backgr %>% dplyr::select(node = from) %>% unique()
  target <- backgr %>% dplyr::select(node = to) %>% unique()
  nodes <- rbind(source, target) %>% unique()
  
  return(nodes)
}

back_nodes_ppi <- preproc_backgr_ppi(backgr_ppi)

# Get all ppi nodes 
nodes_ppi <- nodes_protein %>% dplyr::select(node)  #repeat for mirna, shared or protein only
#nodes_ppi <- nodes_mirna %>% dplyr::select(node)
#nodes_ppi <- nodes_shared %>% dplyr::select(node)

###

# Run GO overenrichment analysis

# Carries out GO overrepresentation analysis.

# Convert to UNIPROT IDs
#nodes_ppi <- bitr(nodes_ppi$node, fromType='SYMBOL', toType='UNIPROT', OrgDb="org.Hs.eg.db") #USE FOR CARNIVAL

# Carry out normal GO gene overrepresentation analysis
#go1 <- enrichGO(gene = as.character(nodes_ppi$UNIPROT), OrgDb = "org.Hs.eg.db", keyType = "UNIPROT", universe =as.character(back_nodes_ppi$node), ont = 'BP', qvalueCutoff = 0.05) # USE FOR CARNIVAL 
go1 <- enrichGO(gene = as.character(nodes_ppi$node), OrgDb = "org.Hs.eg.db", keyType = "UNIPROT", universe =as.character(back_nodes_ppi$node), ont = 'BP', qvalueCutoff = 0.05) #use for ViralLink

# Remove redundancy of GO terms. Cutoff refers to similarity
go2 <- clusterProfiler::simplify(go1, by = "qvalue", select_fun=min)

# Get as dataframe
go2_df <- as.data.frame(go2)

# Get dot plot
dot_plot <- dotplot(go2, showCategory=10, orderBy="qvalue", font.size = 10)
# Save dot plot
filep <- file.path("GO_overrep_protein.pdf")
pdf(filep)
print(dot_plot)
dev.off()

write.table(go2_df, "GO_functional_merged_network_protein.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# Run reactome overenrichment analysis

# Carries out reactome overrepresentation analysis.

# Convert to entrez genes
net_e <- bitr(nodes_ppi$node, fromType='UNIPROT', toType='ENTREZID', OrgDb="org.Hs.eg.db") #use for ViralLink
#net_e <- bitr(nodes_ppi$SYMBOL, fromType='SYMBOL', toType='ENTREZID', OrgDb="org.Hs.eg.db") #use for CARNIVAL

back_nodes_e <- bitr(back_nodes_ppi$node, fromType='UNIPROT', toType='ENTREZID', OrgDb="org.Hs.eg.db")

# Carry out normal reactome gene overrepresentation analysis (can only use entrez :( )
re1 <- enrichPathway(gene = as.character(net_e$ENTREZ), organism = "human", universe =as.character(back_nodes_e$ENTREZ), qvalueCutoff = 0.05)

# Get as dataframe
re1_df <- as.data.frame(re1)
  
# Get dot plot
dot_plot2 <- dotplot(re1, showCategory=10, orderBy="qvalue", font.size = 10)
  
# Save dot plot
filep2 <- file.path("REACT_overrep_protein.pdf")
print(dot_plot2)
pdf(filep2, width=6.5)
dev.off()
  
write.table(re1_df, "REACT_functional_merged_network_protein.txt", quote = FALSE, row.names = FALSE, sep = "\t")
