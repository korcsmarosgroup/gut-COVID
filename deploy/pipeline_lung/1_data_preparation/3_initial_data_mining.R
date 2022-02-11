library(dplyr)

# Capture  messages and errors to a file.
zz <- file("gutcovid_lung.out", open="a")
sink(zz, type = "message", append = TRUE)
message("\nStarting differential expression analysis script: 3_initial_data_mining.R\n")

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 3){
  stop("Wrong number of command line input parameters. Please check.")
}

#setwd("~/OneDrive - Norwich BioScience Institutes/Bioinformatics/Lung_run/")

data <- read.delim(args[1], row.names=1)

data_sars <- data %>% select_if(grepl("_moderate", names(.)))

colnames(data_sars) <- gsub("_moderate", "", colnames(data_sars))

ciliated <- data_sars %>% select("Ciliated")

immune <- data %>% select("B.cell_control", "CTL_control", "MC_control", "MoD.Ma_control", "moDC_control",
                          "NKT_control", "NKT.p_control", "nrMa_control","pDC_control", "rMa_control", 
                          "Treg_control")

write.table(immune, file = args[2], sep = ",",quote=F,row.names=T)
write.table(ciliated, file = args[3], sep = ",",quote=F,row.names=T)
