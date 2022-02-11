# Script to combine the stat2 table and the new categories of receptors to create a ligand-> receptor family network table
#
# Input: stat2 table output from ligand_receptor_stats.R and filtered for ligand cell type of interest
#        table of receptor - receptor family conversion created manually. excel file. columns: receptro, receptor_cat, receptor_cat_sh, count
# Output: Tab delimited table of ligand-> receptor category interactions

##### Set up #####

library(dplyr)
library(readxl)

#setwd("~/OneDrive - Norwich BioScience Institutes/Bioinformatics/Lung_run/ciliated_cells_moderate")

# Capture  messages and errors to a file.
zz <- file("gutcovid_lung.out", open="a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: 16_ligand_receptor_family_network_creation.R\n")

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 3){
  stop("Wrong number of command line input parameters. Please check.")
}

stats <- read.csv(args[1], sep = "\t")
receptor_cats <- read_excel(args[2])

##### Process #####

# Left join
stats_all <- left_join(stats, receptor_cats, by = c("receptor" = "receptor"))

# Collapse receptors
stats_all2 <- stats_all %>% select(-c(receptor_cell_target)) %>%
  group_by(ligand, ligand_cell_source, receptor_cat, receptor_cat_sh, receptor_cat_count = count) %>%
  summarize(mean_num_imm_cells = mean(target_cell_counts, na.rm = TRUE), num_receptors_targeted = n())

# Get sum of incoming interactions for each receptor category (sum num_receptors_targeted)
stats_all_2_col <- stats_all2 %>% group_by(receptor_cat_sh) %>% summarize(sum_incoming_receptors = sum(num_receptors_targeted)) %>%
  mutate(log_sum_incoming_receptors = log(sum_incoming_receptors))

# Get all receptors collapsed by category
receptor_cats2 <- receptor_cats %>% select(-c(receptor_cat_sh, count)) %>%
  group_by(receptor_cat) %>% summarise(receptors_in_group = paste0(receptor, collapse = ","))

# Join to stats
stats_all2 <- left_join(stats_all2, receptor_cats2)
stats_all2 <- left_join(stats_all2, stats_all_2_col)

# Save
write.table(stats_all2, file = args[3], sep = "\t", quote = F, row.names = F)
