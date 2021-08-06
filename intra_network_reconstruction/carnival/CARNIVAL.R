## Libraries and functions

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("piano")

## Set working directory
setwd("~/OneDrive - Norwich BioScience Institutes/PhD files/Bioinformatics/2020_Gut-Covid/results/Gut-Covid_run5/CARNIVAL")

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
  
  return(list(sucesses = sucesses, bg= bg))
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
  readRDS("output_directory/networks/colon_carnival_results_top50tf_pleitropic_minsize15.rds")
nodes_carnival_inhibition <- extractCARNIVALnodes(carnival_results_inhibition)

carnival_results_activation <- 
  readRDS("output_directory/networks/colon_carnival_results_top50tf_pleitropic_minsize15_activation.rds")
nodes_carnival_activation <- extractCARNIVALnodes(carnival_results_activation)

carnival_results_undefined <- 
  readRDS("output_directory/networks/colon_new_08012021/colon_carnival_results_top50tf_pleitropic_minsize15_undefinedEffect.rds")
nodes_carnival_undefined <- extractCARNIVALnodes(carnival_results_undefined)


##Reading Pathway data sets from MSigDB
#We downloaded from MSigDB https://www.gsea-msigdb.org/ the following dataset: c2.cp.v7.0.symbols.gmt. 
##C2 (curated gene sets) -> CP (canonical pathways) -> Gene symbols
#It contains several pathways from different resources and the genes that are known to be involved in 
#those pathways.

Pathway_signatures <- loadGSC("input_files/c2.cp.v7.2.symbols.gmt")

##Reading and formatting statistic from DEG
#We read the results from the differential expression analysis. The statistic of the genes will 
#be mapped later on in the different significant pathways. Maybe, this is not very informative in this context.

DEA_results <- read.csv("input_files/expression_data/unfilt_degs_colon_24h_infected_enterocytes_vs_mock.csv", sep = ",") %>%
  #dplyr::rename(Gene = X) %>% 
    dplyr::select(Gene, log2FoldChange) #before it was "stat" instead of log2FoldChange

## Parsed with column specification:
## cols(
##   Gene = col_character(),
##   baseMean = col_double(),
##   log2FoldChange = col_double(),
##   lfcSE = col_double(),
##   stat = col_double(),
##   pvalue = col_double(),
##   padj = col_double()
## )

mean_stat <- unlist(lapply(Pathway_signatures$gsc, function(x, DEA_results) {
  genes_matching <- x[which(x %in% DEA_results$Gene)]
  mean_genes <- dplyr::filter(DEA_results, Gene %in% genes_matching) %>%
    dplyr::pull(log2FoldChange) %>% mean(na.rm = TRUE) #before it was "stat" instead of log2FoldChange
  return(mean_genes)
}, DEA_results = DEA_results)) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "pathway") %>%
  dplyr::rename(mean_stat = ".")

### Trying with Carnival node activity
carnival_avg_node_act <- carnival_results_undefined$nodesAttributes %>%  #before there was carnival_results_inhibition... don't know if I could change it
  as.data.frame() %>%
  dplyr::mutate(AvgAct = as.numeric(AvgAct))

mean_stat_carni <- 
  unlist(lapply(Pathway_signatures$gsc, function(x, carnival_avg_node_act) {
    genes_matching <- x[which(x %in% carnival_avg_node_act$Node)]
    mean_genes <- dplyr::filter(carnival_avg_node_act, Node %in% genes_matching) %>%
      dplyr::pull(AvgAct) %>% mean(na.rm = TRUE)
    return(mean_genes)
  }, carnival_avg_node_act = carnival_avg_node_act)) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "pathway") %>%
  dplyr::rename(mean_stat = ".")

##Performing Enrichment Analysis and plotting the Results
#Using the Piano R package, we run a gene set analysis (GSA) based on a list of significant genes (CARNIVAL nodes) and a gene set collection (background). It uses Fisher’s exact test.

#First for CARNIVAL results when we assume that the viral-host proteins interactions are inhibitory:
  
  hyper_results_inhibition <- 
  runGSAhyper(genes = nodes_carnival_inhibition$sucesses, 
              # pvalues = rep(0,length(nodes_carnival_inhibition$sucesses)),
              pcutoff = 0.01, 
              universe = nodes_carnival_inhibition$bg, gsc = Pathway_signatures,
              gsSizeLim = c(4,Inf), adjMethod = "fdr")

  ## Warning in runGSAhyper(genes = nodes_carnival_inhibition$sucesses, pcutoff
## = 0.01, : there are genes in gsc that are not in the universe, these will be
## removed before analysis

## Analyzing the overrepresentation of 85 genes of interest in 2054 gene sets, using a background of
#6811 non-interesting genes.

  enriched_pathways_inhibition <- hyper_results_inhibition$resTab %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "pathway") %>% 
  as_tibble() %>% 
  dplyr::filter(`Adjusted p-value` <= 0.001) %>% 
  dplyr::inner_join(mean_stat) %>% 
  dplyr::select(pathway, `p-value`, `Adjusted p-value`, mean_stat) %>%
  dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
  dplyr::mutate(pathway = as.factor(pathway))

  ## Joining, by = "pathway"
interesting_pathways_inhibition <- c()

p_inhibition <- BarplotEnrichment(enriched_pathways_inhibition, 
                                  interesting_pathways_inhibition)

p_inhibition


#Then, for CARNIVAL when we assume that the viral-host proteins interactions are stimulations:
  
  hyper_results_activation <- 
  runGSAhyper(genes = nodes_carnival_activation$sucesses, 
              # pvalues = rep(0,length(nodes_carnival_activation$sucesses)),
              pcutoff = 0.01, 
              universe = nodes_carnival_activation$bg, gsc = Pathway_signatures,
              gsSizeLim = c(4,Inf), adjMethod = "fdr")

## Warning in runGSAhyper(genes = nodes_carnival_activation$sucesses, pcutoff
## = 0.01, : there are genes in gsc that are not in the universe, these will be
## removed before analysis

## Analyzing the overrepresentation of 86 genes of interest in 2054 gene sets, using a background of 
#6810 non-interesting genes.
  
enriched_pathways_activation <- hyper_results_activation$resTab %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "pathway") %>% 
  as_tibble() %>% 
  dplyr::filter(`Adjusted p-value` <= 0.001) %>% 
  dplyr::inner_join(mean_stat) %>% 
  dplyr::select(pathway, `p-value`, `Adjusted p-value`, mean_stat) %>%
  dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
  dplyr::mutate(pathway = as.factor(pathway))

## Joining, by = "pathway"
interesting_pathways_activation <- c()

p_activation <- BarplotEnrichment(enriched_pathways_activation, 
                                  interesting_pathways_activation)
p_activation

#Then, for CARNIVAL when we assume that the viral-host proteins interactions are undefined and we 
#let CARNIVAL infer the effect: (Note: For visualization purposes, we reduce here the significant 
#p-value to 0.0001)

hyper_results_undefined <- 
  runGSAhyper(genes = nodes_carnival_undefined$sucesses, 
              # pvalues = rep(0,length(nodes_carnival_activation$sucesses)),
              pcutoff = 0.01, 
              universe = nodes_carnival_undefined$bg, gsc = Pathway_signatures,
              gsSizeLim = c(4,Inf), adjMethod = "fdr")

## Warning in runGSAhyper(genes = nodes_carnival_undefined$sucesses, pcutoff =
## 0.01, : there are genes in gsc that are not in the universe, these will be
## removed before analysis

## Analyzing the overrepresentation of 141 genes of interest in 2054 gene sets, using a background 
#of 6741 non-interesting genes.

enriched_pathways_undefined <- hyper_results_undefined$resTab %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "pathway") %>% 
  as_tibble() %>% 
  dplyr::filter(`Adjusted p-value` <= 0.000001) %>% 
  dplyr::inner_join(mean_stat) %>% 
  dplyr::select(pathway, `p-value`, `Adjusted p-value`, mean_stat) %>%
  dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
  dplyr::mutate(pathway = as.factor(pathway))

## Joining, by = "pathway"
interesting_pathways_undefined <- c()

p_undefined <- BarplotEnrichment(enriched_pathways_undefined, 
                                 interesting_pathways_undefined)
p_undefined

##Comparing Enrichment between activation and inhibition
#Now I compare the enriched pathways that change the most between both assumptions: 
#host-viral interactions are inhibitory or stimulatory.

results_activation <-   hyper_results_activation$resTab %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "pathway") %>% 
  as_tibble() %>% 
  dplyr::mutate(LogpValue_activation = -log(`p-value`), 
                p_value_activation = `p-value`) %>% 
  dplyr::select(pathway, p_value_activation, LogpValue_activation) 

results_inhibition <-   hyper_results_inhibition$resTab %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "pathway") %>% 
  as_tibble() %>% 
  dplyr::mutate(LogpValue_inhibition = -log(`p-value`), 
                p_value_inhibition = `p-value`) %>% 
  dplyr::select(pathway, p_value_inhibition, LogpValue_inhibition) 

results_join <- 
  dplyr::inner_join(results_activation, results_inhibition) %>% 
  dplyr::filter(p_value_activation < 0.05 | p_value_inhibition < 0.05) %>% 
  dplyr::mutate(diff = LogpValue_activation - LogpValue_inhibition) %>% 
  dplyr::top_n(15, wt = abs(diff)) %>% 
  dplyr::arrange(desc(abs(diff)))

## Joining, by = "pathway"
results_join_longer <- results_join %>% 
  dplyr::select(pathway, diff, LogpValue_activation, LogpValue_inhibition) %>% 
  dplyr::rename(activation = "LogpValue_activation", 
                inhibition = "LogpValue_inhibition") %>% 
  pivot_longer(-c(pathway,diff), values_to = "LogPvalue", 
               names_to = "Interaction")  

point_plot_activation_vs_inhibtion <- 
  ggplot(results_join_longer, aes(reorder(pathway, abs(diff)), LogPvalue)) + 
  geom_point(aes(color = Interaction), size = 3) + 
  coord_flip() + 
  theme_minimal() + 
  theme(legend.position = "bottom",  legend.justification = "center") +
  scale_color_lancet() +
  theme(axis.text.x = element_text(hjust = 1, size =8.5),
        axis.text.y = element_text(size = 7),
        panel.grid.minor = element_blank()) + 
  xlab("Pathways") + ylab("-Log (p-value)") +
  geom_hline(yintercept = -log(0.05), linetype="dashed", 
             color = "#2F4F4F", size=0.5)

point_plot_activation_vs_inhibtion



##Clustering of CARNIVAL solution and enrichment.
#For the clustering part, we are going to take our CARNIVAL results when the effect of the 
#viral proteins on the host proteins is undefined. We first convert our CARNIVAL output network 
#into an igraph object. Then, we perform different types of clustering and select the one with the 
#largest modularity for further analysis.

carnival_results_igraph <- carnival_results_undefined$weightedSIF %>%
  as.data.frame() %>% dplyr::select(-Sign) %>% 
  dplyr::rename(weight = "Weight") %>% 
  igraph::graph_from_data_frame(directed = FALSE)

c_edgeBetweeness <- 
  cluster_edge_betweenness(carnival_results_igraph)

c_FastGreedy <- 
  cluster_fast_greedy(carnival_results_igraph, 
                      merges = TRUE, 
                      modularity = TRUE, 
                      membership = TRUE, 
                      weights = E(carnival_results_igraph)$weight)

c_Infomap <- 
  cluster_infomap(carnival_results_igraph, 
                  e.weights =  E(carnival_results_igraph)$weight, 
                  v.weights = NULL, nb.trials = 10, modularity = TRUE)

c_labelProp <- 
  cluster_label_prop(carnival_results_igraph, 
                     weights = E(carnival_results_igraph)$weight)

c_leadingEig <- 
  cluster_leading_eigen(carnival_results_igraph, 
                        weights = E(carnival_results_igraph)$weight)

c_Louvain<- 
  cluster_louvain(carnival_results_igraph, 
                  weights = E(carnival_results_igraph)$weight)

c_walktrap <- 
  cluster_walktrap(carnival_results_igraph, 
                   weights = E(carnival_results_igraph)$weight, 
                   steps = 4, merges = TRUE, modularity = TRUE, membership = TRUE)

modularity_df <- data.frame(method = c("EdgeBetweeness", "FastGreedy", 
                                       "Infomap","LabelPropagation", "Louvain","Walktrap"), 
                            modularity = c(max(c_edgeBetweeness$modularity), 
                                           max(c_FastGreedy$modularity), 
                                           max(c_Infomap$modularity), 
                                           max(c_labelProp$modularity),
                                           #max(c_leadingEig$modularity),
                                           max(c_Louvain$modularity),
                                           max(c_walktrap$modularity))) %>% 
  dplyr::arrange(desc(modularity))

modularity_df

##method modularity

#colon
##1         FastGreedy  0.4794261
##2            Louvain  0.4794261
##3            Infomap  0.4737659
##4 LeadingEigenvector  0.4725587
##5   LabelPropagation  0.4551593
##6           Walktrap  0.4531015
##7     EdgeBetweeness  0.3938351

#ileum -> I did not use the leadingEig method because it was giving me errors

##1       FastGreedy  0.5616948
##2          Louvain  0.5542273
##3          Infomap  0.5355343
##4   EdgeBetweeness  0.5270886
##5         Walktrap  0.5000915
##6 LabelPropagation  0.4446930

##We then selected the communities generated by the Fast Greedy method (highest modularity). We will perform enrichment 
##analysis on them.

## Number of genes in the different communities. 

table(c_FastGreedy$membership)

#colon
#1  2  3  4  5  6 
#15  8  6  8  8  2 

#ileum
#1  2  3  4  5  6  7  8 
#17  8  8 16 13  8  6  2 

## We perform an Enrichment Analysis for each Community: 

n <- length(unique(c_FastGreedy$membership))
gene_communities <- data.frame(Gene = c_FastGreedy$names, Community = c_FastGreedy$membership)

p <- list()

for (i in seq_len(n)){
  current_genes <- dplyr::filter(gene_communities, Community==i) %>%
    dplyr::pull(Gene)
  current_result <- gost(current_genes, user_threshold = 0.01, 
                         correction_method = "fdr", custom_bg = nodes_carnival_undefined$bg,
                         sources = c("GO","KEGG","REAC","WP")) 
  if (!is.null(current_result)){
    #currentfile <- 
    #    paste0("02_analysis_carnival_results_files/figure-gfm/enrichment_cluster_", i, ".png")
    p[[i]] <- gostplot(current_result, capped = FALSE, interactive = TRUE)
    
  }
}

##We finally show the enrichment results for the different clusters. 

#Enrichment Cluster 1:
  
dplyr::filter(gene_communities, Community == 1) %>% dplyr::pull(Gene)

#colon
#[1] "ATM"     "YAP1"    "CTNNA2"  "HIF1A"   "AKT1"    "ERBB2"   "MUL1"    "PTPRF"   "ARHGEF2" "IMPDH2"  "ITGA6"  
#[12] "TBK1"    "S1PR3"   "TEAD1"   "ETS1"   

#ileum
#[1] "NFKB1"   "FOS"     "RELA"    "EXTL3"   "MAPK1"   "CSNK2A1" "NRG1"    "BSG"     "TBK1"    "ARF6"   
#[11] "PNPLA2"  "MET"     "PLAT"    "S1PR3"   "PTPRJ"   "ETS1"    "CEBPA"  

p[[1]] 

##Enrichment Cluster 2:

dplyr::filter(gene_communities, Community == 2) %>% dplyr::pull(Gene)

#colon
#[1] "TP53BP1" "PRKCD"   "RHOA"    "GSK3B"   "CEBPB"   "PRKACA"  "E2F4"    "KLF5"   

#ileum
#[1] "NOTCH1" "SMAD4"  "SMAD2"  "TGFB2"  "SOX2"   "MEF2C"  "MEF2A"  "E2F4"  

p[[2]] 

##Enrichment Cluster 3:
dplyr::filter(gene_communities, Community == 3) %>% dplyr::pull(Gene)

#colon
#[1] "JAK2"  "STAT2" "STAT1" "GSK3A" "PUF60" "MYC"  

#ileum
#[1] "PPP2CA"  "TGFBR2"  "PPP2R1A" "MAPK12"  "JAK2"    "POFUT1"  "PPP2CB"  "ELF1"   

p[[3]] 

##Enrichment Cluster 4:

dplyr::filter(gene_communities, Community == 4) %>% dplyr::pull(Gene)

#colon
#[1] "CDK1"   "E2F1"   "TP53"   "SSBP1"  "PPP1CB" "PRKDC"  "PARP1"  "ATR"  

#ileum
#[1] "GSK3B"    "OS9"      "EGLN1"    "EPAS1"    "PRKCD"    "RHOA"     "ARNT"     "AKT1"     "IMPDH2"  
#[10] "EPHB3"    "TXN"      "TNFRSF1A" "HEY2"     "MARK2"    "KLF5"     "BHLHE40" 

p[[4]] 

##Enrichment Cluster 5:

dplyr::filter(gene_communities, Community == 5) %>% dplyr::pull(Gene)

#colon
#[1] "NFKB1" "RELA"  "NKRF"  "BLVRA" "IRAK1" "RIPK1" "F2RL1" "JUN" 

#ileum
# [1] "PRKACA" "ATF2"   "NR2C2"  "BMPR2"  "EP300"  "BMPR1A" "PRKAB1" "TGFB1"  "FASN"   "ATM"    "BLVRA"  "ATF3"  
#[13] "ESR1"  

p[[5]] 

##Enrichment Cluster 6:

dplyr::filter(gene_communities, Community == 6) %>% dplyr::pull(Gene)

#colon
#[1] "GDF15" "EGR1" 

#ileum
#[1] "STAT3"  "STAT2"  "STAT1"  "NOP2"   "PLAU"   "PTPN11" "SPI1"   "IRF2"  

p[[6]]  


dplyr::filter(gene_communities, Community == 7) %>% dplyr::pull(Gene)

#ileum
#[1] "ATR"   "E2F1"  "NF1"   "TP53"  "KMT2A" "TAF1" 

p[[7]]

dplyr::filter(gene_communities, Community == 8) %>% dplyr::pull(Gene)

#ileum
#[1] "GSK3A" "TCF4"  

p[[8]]

