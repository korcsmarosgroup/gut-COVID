## Libraries and functions

# Capture  messages and errors to a file.
zz <- file("gutcovid.out", open = "a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: colon_counts_table_from_10x_data.R\n")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("piano")

## Set working directory
setwd("/home/gutcovid/pipeline/8_intra_output_data/carnival/")

#We first load the required libraries and some important functions to analyse CARNiVAL output.

library(readr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(piano) #install
library(tidyr)
library(ggsci)
library(igraph)
library(gprofiler2) #install
library(plotly) #install

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 1){
  stop("Wrong number of command line input parameters. Please check.")
}


## Function to extract the nodes that appear in CARNIVAL network and the 
## background genes (all genes present in the prior knowledge network).
## It returns a list with two objects: the success and the background genes.

extractCARNIVALnodes <- function(CarnivalResults){
  
  CarnivalNetwork <- 
    as.data.frame(CarnivalResults$weightedSIF, stringsAsFactors = FALSE)
  
  colnames(CarnivalNetwork) <- c("source", "sign", "target", "Weight")
  
  ## We define the set of nodes interesting for our condition
  sucesses <- unique(c(gsub("_.*","",CarnivalNetwork$source), 
                       gsub("_.*","",CarnivalNetwork$target)))
  
  CarnivalAttributes <- as.data.frame(CarnivalResults$nodesAttributes, 
                                      stringsAsFactors = FALSE)
  
  ## We define the background as all the genes in our prior knowledge network.
  bg <- unique(gsub("_.*","",CarnivalAttributes$Node))
  
  return(list(sucesses = sucesses, bg = bg))
}

### Function to print a barplot with the enriched pathways.
BarplotEnrichment <- function(PathwaysSelect, Interesting_pathways){ 
  
  p <- ggplot(PathwaysSelect, aes(x = reorder(pathway, pvalue), 
                                  y = -log10(pvalue))) + 
    geom_bar(aes(fill = mean_stat), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
                         mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.text.x = element_text(size = 4, angle = 90, hjust = 1,
                                     colour = ifelse(levels(reorder(PathwaysSelect$pathway, 
                                                                    PathwaysSelect$pvalue)) %in% Interesting_pathways, 
                                                     "red", "grey40"),
                                     face = ifelse(levels(reorder(PathwaysSelect$pathway, 
                                                                  PathwaysSelect$pvalue)) %in% Interesting_pathways, 
                                                   "bold", "plain")),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    xlab("")
  return(p)
}    

##Reading CARNIVAL results
#We then read the CARNIVAL results generated in the previous script. We define two different 
#gene sets in order tor conduct the enrichment. The first set contains the nodes that appear in the 
#CARNIVAL output and are therefore relevant in the context of our input transcriptomic data. 
#The second set contains all the genes in our prior knowledge network which are used as the backgroud.

carnival_results_inhibition <- 
  readRDS(args[1])
nodes_carnival_inhibition <- extractCARNIVALnodes(carnival_results_inhibition)

# carnival_results_activation <- 
#   readRDS("../../8_intra_output_data/carnival/colon_carnival_results_top50tf_pleitropic_minsize15_activation.rds")
# nodes_carnival_activation <- extractCARNIVALnodes(carnival_results_activation)

# carnival_results_undefined <- 
#   readRDS("../../8_intra_output_data/carnival/colon_carnival_results_top50tf_pleitropic_minsize15_undefinedEffect.rds")
# nodes_carnival_undefined <- extractCARNIVALnodes(carnival_results_undefined)
