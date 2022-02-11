## Introduction

#The present script takes the results of the differential expression analysis on
#transcriptomics data from intestinal organoids treated with Sars-CoV-2 for 60 
#hours and the following tools: Progeny, Dorothea and CARNIVAL (see References). 

#The original transcriptomic dataset is coming from the following publication:
  #[Lamers et al. 2020](https://science.sciencemag.org/content/369/6499/50).

#The differential expression analysis was conducted by Martina Poletti (<Martina.Poletti@earlham.ac.uk>) 

### Getting Ready

# Capture  messages and errors to a file.
zz <- file("gutcovid.out", open = "a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: Classical_pipeline_viral_proteins.R\n")

#install requires packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("progeny")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dorothea")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("CARNIVAL")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 6){
  stop("Wrong number of command line input parameters. Please check.")
}

#We first load the required libraries and we define a function to export CARNIVAL
#results to Cytoscpae. 

library(tibble)
library(dplyr)
library(readr)
library(progeny) #to install
library(dorothea) #to install
library(CARNIVAL) #to install
library(ggplot2)
library(OmnipathR)
library(biomaRt) #to install
library(grid)
library(gridExtra)
library(magrittr)

## We also define a function to format the CARNIVAL output to cytoscape
OutputCyto <- function(CarnivalResults, outputFile) {
  CarnivalNetwork <- 
    as.data.frame(CarnivalResults$weightedSIF, stringsAsFactors = FALSE) %>%
    dplyr::mutate(Sign = as.numeric(Sign), Weight = as.numeric(Weight)) %>% 
    dplyr::mutate(Weight = Sign * Weight) %>%
    dplyr::select(Node1, Weight, Node2)
  
  CarnivalNetworkNodes <- 
    unique(c(CarnivalNetwork$Node1,CarnivalNetwork$Node2))
  
  CarnivalAttributes <- CarnivalResults$nodesAttributes %>% 
    as.data.frame() %>%
    dplyr::filter(Node %in% CarnivalNetworkNodes) %>%
    dplyr::mutate(NodeType = as.character(NodeType)) %>%
    dplyr::mutate(NodeType=if_else(NodeType =="", "I", NodeType))
  
  nameOutputNetwork <- paste0(outputFile, "Network.sif")
  nameOutputAttributes <-  paste0(outputFile, "Attributes.txt")    
  
  write.table(CarnivalNetwork, file = nameOutputNetwork,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
  
  write.table(CarnivalAttributes, file = nameOutputAttributes,
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}

#We also read the results from the differential expression analysis: 
  
## There is a duplicate RPL21. I am going to remove it in a first approach

all_genes_dea <- read.csv(args[1], sep = ",")
all_genes_dea <- all_genes_dea[,1:6] %>% distinct(Gene, .keep_all = TRUE)

## Results 

### Pathway activity estimation using Progeny

#We first estimate the pathway activity using the Progeny package. In particular,
#we compute the normalised enriched score (NES) of the different pathways by 
#running Progeny using the statistic from the differential express analysis.

all_genes_dea_stat <- all_genes_dea %>% 
  dplyr::select(Gene, log2FoldChange) %>% #log2FoldChange was "stat" before
  dplyr::filter(!is.na(log2FoldChange)) %>% #log2FoldChange was "stat" before
  column_to_rownames(var = "Gene") 
pathways_zscore <- t(progeny(as.matrix(all_genes_dea_stat), 
                             scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE))
colnames(pathways_zscore) <- "NES"

## I also need to run progeny in such a way to have values between 1 and -1 to
## use as CARNIVAL input

pathways_inputCarnival <- 
  t(progeny(as.matrix(all_genes_dea_stat), 
            scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE))
colnames(pathways_inputCarnival) <- "Activity"


#We now display the normalized enrichment scores (NES) in a bar plot.

pathways_zscore_df <- as.data.frame(pathways_zscore) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))
ggplot(pathways_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")


#It is quite shocking that the MAPK pathway in un-active. This pathway was 
#usually found to be more active after SARS-CoV-2 infection in other studies.

weight_matrix <- getModel("Human", top=100)
weight_matrix <- data.frame(names = row.names(weight_matrix), 
                            row.names = NULL, weight_matrix)
plots <- progenyScatter(all_genes_dea_stat %>% 
                          tibble::rownames_to_column(var = "GeneID"), weight_matrix)

#We therefore check which genes are contributing the most to the MAPK progeny 
#score. In these plots, we have the progeny scores for each pathway and stats
#of the different genes. In red, we will have the genes contributing to the 
#positive activation score and on blue the opposite. 

#plot.new() run this to create a white backggrund to save the plots

grid.draw(plots[[1]]$NFkB)

grid.draw(plots[[1]]$TNFa)

grid.draw(plots[[1]]$JAK.STAT)

### Transcription factor activity with Dorothea and Viper. 

#Now, we estimate the transcription factor (TF) activity using the dorothea 
#package. We select Dorothea interactions with confidence level A, B and C.


## We load Dorothea Regulons
#data(dorothea_hs, package = "dorothea")

#regulons <- dorothea_hs %>%
#dplyr::filter(confidence %in% c("A", "B","C"))

#we load our pre-processed dorothea regulons
regulons <- read.csv(args[2], sep="") %>% 
  rename(tf = from) %>% 
    rename(target = to)

#We run Viper using the statistic from the different expression analysis. First,
#we run it considering TF with at least 5 targets and with no correction for 
#pleiotropic regulation

all_genes_dea_stat <-  all_genes_dea %>% 
  dplyr::select(Gene, log2FoldChange) %>%  #log2FoldChange was "stat" before
  dplyr::filter(!is.na(log2FoldChange)) %>%  #log2FoldChange was "stat" before
  column_to_rownames(var = "Gene") %>%
  as.matrix()

tf_activities_stat <- 
  dorothea::run_viper(all_genes_dea_stat, regulons,
                      options =  list(minsize = 5, eset.filter = FALSE, 
                                      cores = 8, verbose = FALSE, nes = TRUE))

#We now display the top 25 normalized enrichment scores (NES) for the TF in a 
#bar plot.

tf_activities_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "Tf") %>%
  dplyr::rename(NES = "log2FoldChange") %>% #log2FoldChange was "stat" before
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(Tf = factor(Tf))

ggplot(tf_activities_top25,aes(x = reorder(Tf, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") + #here stat is the statistics, so it doesn't need to be changed
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")

#I run again Viper with a more conservative setup. The TFs need to regulate at 
#least 15 targets genes and I include the correction for pleiotropic regulation 

tf_activities_stat_pleitropic <- 
  dorothea::run_viper(all_genes_dea_stat, regulons,
                      options =  list(minsize = 15, eset.filter = FALSE, 
                                      cores = 8, verbose = FALSE, nes = TRUE, pleiotropy= TRUE))
tf_activities_top25_pleitropic <- tf_activities_stat_pleitropic %>%
  as.data.frame() %>% 
  rownames_to_column(var = "Tf") %>%
  dplyr::rename(NES = "log2FoldChange") %>% #log2FoldChange was "stat" before
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(Tf = factor(Tf))
ggplot(tf_activities_top25_pleitropic,aes(x = reorder(Tf, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")


### Running CARNIVAL

#We use the OmnipathR package to fetch the Omnipath database and generate the 
#prior knowledge network. We take the signed and directed protein-protein 
#interactions. 

#ia_omnipath <- import_omnipath_interactions() %>% as_tibble()
#ia_pwextra <- import_pathwayextra_interactions() %>% as_tibble()
#ia_kinaseextra <- import_kinaseextra_interactions() %>% as_tibble()

## We bind the datasets

#interactions2 <- as_tibble(
  #bind_rows(
    #ia_omnipath %>% mutate(type = 'ppi'),
    #ia_pwextra %>% mutate(type = 'ppi'),
    #ia_kinaseextra %>% mutate(type = 'ppi')))

#import our pre-processed OmniPath network
interactions <- read.csv(args[3], sep ="") 
  
#some of these steps may already been done on the imported omnipath network
#I wasn't sure so I repeated them

signed_directed_interactions <- 
  dplyr::filter(interactions, consensus_direction==1) %>%
  filter(consensus_stimulation == 1 | consensus_inhibition == 1) %>% 
  dplyr::mutate(sign = if_else(consensus_stimulation==1,1,-1))  %>%
  dplyr::select(source_genesymbol, sign,  target_genesymbol) %>%
  dplyr::rename(source ="source_genesymbol", target ="target_genesymbol")

carnival_pkn <- signed_directed_interactions %>%
  dplyr::distinct(source, target, .keep_all = TRUE)

all_source_nodes <- unique(carnival_pkn$source)
all_target_nodes <- unique(carnival_pkn$target)
all_nodes_network <- unique(c(all_source_nodes,all_target_nodes))

#We are going to explore the host-virus interactions from IntAct. 

host_viral_interactions <- 
  read.csv(args[4], sep = "\t")

#we convert from Uniprot to gene name (ADDITION - the script below did not work so I substituted it)

library("org.Hs.eg.db") 

host_viral_interactions$hgnc_symbol <- mapIds(org.Hs.eg.db, keys = host_viral_interactions$human_protein, keytype = "UNIPROT", column="SYMBOL")
host_viral_interactions_hgnc <- host_viral_interactions

## We translate the human uniprot symbols to hgnc. 
#ensembl <- useMart('ensembl', dataset="hsapiens_gene_ensembl")

# listAttributes(ensembl)
#uniprot_hgnc <- 
  #getBM(attributes=c("uniprotswissprot", "hgnc_symbol"),  
        #filters="uniprotswissprot",
        #values=unique(host_viral_interactions$human_protein), mart=ensembl)

#host_viral_interactions_hgnc <- 
  #dplyr::inner_join(host_viral_interactions, uniprot_hgnc, 
                    #by = c("human_protein" = "uniprotswissprot"))

## These are the number of interactions: 
nrow(host_viral_interactions_hgnc)

## Are these human proteins sources of any interaction in our prior knowledge
## network? Otherwise it does not make sense to include them as perturbartions 
## of the signaling. We assume that they are inhibitions hampering the regular
## signaling 

host_viral_interactions_hgnc_filter <- host_viral_interactions_hgnc %>% 
  dplyr::filter(hgnc_symbol %in% all_source_nodes) %>% 
  dplyr::select(viral_protein, hgnc_symbol) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(sign = -1) %>% 
  dplyr::rename(source = "viral_protein", target = "hgnc_symbol") %>%
  dplyr::select(source, sign, target)

#We finally merge the virus-host interactions with our prior knowledge network
#and we define the viral proteins as perturbations. We also prepare TF activity
#scores and progeny weights. In particular, we are going to select the top 50 TFs
#when running viper with the pleitropic correction and the TFs should have at 
#least 15 targets. 

## Prior Knowledge Network
carnival_pkn_hostvirus <- 
  bind_rows(carnival_pkn, host_viral_interactions_hgnc_filter)

## Perturbation (Viral proteins are active)
viral_proteins <- unique(host_viral_interactions_hgnc_filter$source)

viral_proteins_perturbation <- 
  data.frame(viral_proteins = viral_proteins, sign = 1) %>% 
  tibble::column_to_rownames(var = "viral_proteins") %>%
  t() %>% as.data.frame()

## Top TFs

tf_top50_pleitropic <- tf_activities_stat_pleitropic %>%
  as.data.frame() %>% 
  rownames_to_column(var = "TF") %>%
  dplyr::filter(TF %in% all_nodes_network) %>% 
  dplyr::rename(NES = "log2FoldChange") %>% # here log2FoldChange was "stat" before
  dplyr::top_n(50, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>%
  tibble::column_to_rownames(var = "TF")  %>%
  t() %>% as.data.frame()

## Progeny Weigths
progeny_weigths <- pathways_inputCarnival %>% t()

#And we run CARNIVAL: 

carnival_results_top50tf_pleitropic_minsize15 <-runCARNIVAL(
  solverPath="/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex",
  netObj=carnival_pkn_hostvirus,
  measObj=tf_top50_pleitropic,
  inputObj = viral_proteins_perturbation,
  #dir_name="output_directory/CARNIVAL",
  weightObj=progeny_weigths,
  #nodeID = 'gene',
  timelimit = 7200,
  solver = "cplex")

saveRDS(carnival_results_top50tf_pleitropic_minsize15, 
        file = args[5]

OutputCyto(carnival_results_top50tf_pleitropic_minsize15, 
           outputFile = args[6])


#Network with the results: Rectangles are the most active transcription factors 
#after infection and the inverse triangles are the perturbed nodes. Ellipses 
#are signaling intermediates proteins linking those perturbations and TFs. 
#Red means activation after infection and blue the opposite. 


#  I am also going to run CARNIVAL by assuming that the virus-host protein 
#interactions are activatory. 

# host_viral_interactions_hgnc_filter_viralActivation <- 
#   host_viral_interactions_hgnc_filter %>%
#   dplyr::mutate(sign = 1) 

# ## Prior Knowledge Network
# carnival_pkn_hostvirus <- 
#   bind_rows(carnival_pkn, host_viral_interactions_hgnc_filter_viralActivation)

# carnival_results_top50tf_pleitropic_minsize15_activation <-runCARNIVAL(
#   solverPath="/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex",
#   netObj=carnival_pkn_hostvirus,
#   measObj=tf_top50_pleitropic,
#   inputObj = viral_proteins_perturbation,
#   # dir_name="Carnival_Results",
#   weightObj=progeny_weigths,
#   # nodeID = 'gene',
#   timelimit = 7200,
#   solver = "cplex")

# saveRDS(carnival_results_top50tf_pleitropic_minsize15_activation, 
#         file = "/home/gutcovid/8_intra_output_data/carnival/carnival_proteins_results_top50tf_pleitropic_minsize15_activation.rds")
# OutputCyto(carnival_results_top50tf_pleitropic_minsize15_activation, 
#            outputFile="/home/gutcovid/8_intra_output_data/carnival/carnival_proteins_results_top50tf_pleitropic_minsize15_activation")
  

# #Finally, I am going to run CARNIVAL defining the perturbations (human proteins
# #that interact with the viral proteins) but without defining their effect 
# #(stimulation or inhibition). CARNIVAL will infer the effect. 

# human_proteins_undefined_perturbation <- 
#   host_viral_interactions_hgnc_filter %>% 
#   dplyr::mutate(sign = NaN) %>%
#   dplyr::select(target, sign) %>% 
#   dplyr::mutate(sign = as.numeric(sign)) %>% 
#   dplyr::distinct() %>% 
#   dplyr::filter(target %in% all_source_nodes) %>%
#   tibble::column_to_rownames(var = "target") %>% 
#   t() %>% as.data.frame()

# carnival_results_top50tf_pleitropic_minsize15_undefinedEffect <-runCARNIVAL(
#   solverPath="/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex",
#   netObj=carnival_pkn,
#   measObj=tf_top50_pleitropic,
#   inputObj = human_proteins_undefined_perturbation,
#   # dir_name="Carnival_Results",
#   weightObj=progeny_weigths,
#   # nodeID = 'gene',
#   timelimit = 14400,
#   solver = "cplex")

# saveRDS(carnival_results_top50tf_pleitropic_minsize15_undefinedEffect, 
#         file = "/home/gutcovid/8_intra_output_data/carnival/colon_carnival_proteins_results_top50tf_pleitropic_minsize15_undefinedEffect.rds")
# OutputCyto(carnival_results_top50tf_pleitropic_minsize15_undefinedEffect, 
#            outputFile="/home/gutcovid/8_intra_output_data/carnival/colon_carnival_proteins_results_top50tf_pleitropic_minsize15_undefinedEffect")
