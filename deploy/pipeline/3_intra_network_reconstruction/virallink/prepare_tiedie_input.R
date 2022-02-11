# Author: Agatha Treveil
# Date: June 2020
#
# Script to prepare the input files for TieDIE
#
# Input: Contextualised PPI network output from 'filter_network_expressed_genes.R'
#        TF-DEGs network file output from 'get_regulator_deg_network.R'
#        Differentially expressed gene file output from 'diff_expression_deseq2.R'
#        Viral protein- human binding protein file (provided - based on Gordon et al.)
#
# Output: pathway.sif - sif format contextualised PPI (Omnipath) network
#         upstream.input - text file containing the upstream genes with a weight and direction (human binding partners)
#         downstream.intput - text file containing the downstream genes with a weight and direction (TFs)

##### Set up #####

# Capture  messages and errors to a file.
zz <- file("virallink.out", open="a")
sink(zz, type="message", append = TRUE)
message("\nStarting tiedie input preparation: prepare_tiedie_input.R\n")

# Install required packages
if (!requireNamespace("tidyverse", quietly = TRUE)) 
install.packages("tidyverse")

# Load packages
library(tidyverse)
library(dplyr)

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 5){
  stop("Wrong number of command line input parameters. Please check.")
}

# outdir <- "/home/gutcovid/8_intra_output_data/virallink/ileum_protein"
# ppis <- read.csv("/home/gutcovid/8_intra_output_data/virallink/ileum_protein/2_process_a_priori_networks/dorothea_contextualised_network.txt", sep = "\t")
# tfs <- read.csv("/home/gutcovid/8_intra_output_data/virallink/ileum_protein/2_process_a_priori_networks/contextualised_regulator-deg_network.txt", sep = "\t")
# degs <- read.csv("/home/gutcovid/1_input_data/intracellular_networks/expression_data/filt_degs_ileum_24h_infected_enterocytes_vs_mock_padj0.05_lfc0.5.csv", sep = ",")
# hbps <- read.csv("/home/gutcovid/1_input_data/intracellular_networks/sarscov2-human_viral_ppis.txt", sep = "\t")


outdir <- args[5]

# Contextualised PPI network (omnipath) output from 'filter_network_expressed_genes.R'
ppis <- read.csv(args[1], sep = "\t")

# TF-DEG interactions output from 'get_regulator_deg_network.R'
tfs <- read.csv(args[2], sep = "\t")

# Filtered differential expression file output from 'diff_expression_deseq2.R'
degs <- read.csv(args[3], sep = ",")

# Human binding proteins of viral proteins
hbps <- read.csv(args[4], sep = "\t")

# Create output dir if required
path <- file.path(outdir, "3_network_diffusion", "input_files")
dir.create(path, showWarnings = FALSE, recursive=TRUE)


##### Process pathways #####

# Convert to sif format
ppis2 <- ppis %>% mutate(direction = ifelse(consensus_stimulation == "1", "stimulates>", "inhibits>")) %>%
  select(c(source_genesymbol,consensus_stimulation, target_genesymbol)) %>%
  unique()

# Save pathways
write.table(ppis2, file = file.path(path,"pathway.sif"), sep = "\t", col.names = F, row.names = F, quote = F)

##### Process upstream input #####

# Number of viral proteins bound to each human protein is the weight
# If sign of interaction given then use it, else assume all inhibitory ("-")
if ("sign" %in% colnames(hbps)){
  # Check values in sign column are "-" or "+"
  hbps_f <- hbps %>% dplyr::filter((sign == "-") | (sign == "+"))
  if (nrow(hbps_f) != nrow(hbps)){
    print("WARNING: Some of the viral-human binding protein interactions were disgarded as the values in the 'sign' column were not '+' or '-'.")
  }
  if (nrow(hbps_f) == 0){
    stop("ERROR: viral-human binding protein interactions do not have the correct values in 'sign' column. They should be '+' or '-'.")
  }
  hbps2 <- hbps %>% select(human_protein, sign) %>% dplyr::rename(direction=sign) %>% group_by(human_protein,direction) %>% summarise(n = n()) %>%
    select(human_protein, n, direction)
} else {
  hbps2 <- hbps %>% select(human_protein) %>% group_by(human_protein) %>% summarise(n = n()) %>% mutate(direction = "-")
}

# Save upstream data
write.table(hbps2, file = file.path(path,"upstream.input"), sep = "\t", col.names = F, row.names = F, quote = F)

##### Process downstream input #####

# (1/#targetgenes) * sum(lfc(targetgene)*signofint)

# Join the tf-deg network with the deg lfc values
tfs2 <- left_join(tfs, degs, by =c("target_genesymbol"="Gene")) %>% select(c(source_genesymbol, target_genesymbol, consensus_stimulation, log2FoldChange)) %>%
  mutate(lfc_sign = ifelse(consensus_stimulation == "1", log2FoldChange, -log2FoldChange))

# Get the sum of all lfc*sign values - and the number of target genes for each tf
tfs3 <- tfs2 %>% select(source_genesymbol, lfc_sign) %>% group_by(source_genesymbol) %>% summarise(sumof = sum(lfc_sign), n = n())

# Divide sumof by n and determine sign (based on sign of the value)
tfs4 <- tfs3 %>% mutate(final_val = sumof/n) %>% mutate(sign = ifelse((final_val >= 0), "+", "-")) %>%
  select(source_genesymbol, final_val, sign)

# Save downstream data
write.table(tfs4, file = file.path(path,"downstream.input"), sep = "\t", col.names = F, row.names = F, quote = F)

# reset message sink and close the file connection
sink(type="message")
close(zz)
