This repository contains all the input files, pipelines and generated intracellular and extracellular network relative to the following publication

#### <i>Poletti et al. (2021) </i>, Reprogramming of the intestinal epithelial-immune cell interactome during SARS-CoV-2 infection

Specifically, the contained code is used to process the single cell RNAseq epithelial and immune data of epithelial and use it as input to build intracellular and intercellular epithelial-immune cell interaction networks upon SARS-CoV-2 infection. 

This repository is divided into different folders, whereby the suffix “intra” or “extra” indicates whether the pipeline refers to the intracellular or intercellular networks. For each set of folders, there is one folder for the input data, one for the network reconstruction part and one for the network analysis part. Additionally, one folder contains scrips used to process the epithelial and immune cell data prior to newtork reconstruction and analysis. Please note that some manual/non-scripted analysis was carried out alongside these scripts, for example the merging of the viral_protein and viral_miRNAs networks was carried out using the merge function in Cytoscape. Description of the different folders and scripts is provided below.

### Inter_input_data

This folder contains all the necessary epithelial and immune cell input data used for the construction of intercellular networks for bystaster/infected epithelial cells in both colon and ileum.

### data_pre-processing

This folder contains all the scripts used to process immune cell and epithelial single cell RNAseq data for consequent construction of intracellular and extracellular networks construction. 

* <i>{ileum|colon}_counts_table_from_10x_data.R</i>| Script used to extract R object from Triana et al., 2020, the average expression / counts for each sample as well as carry out differential expression analysis (using MAST).

* <i>intercellular_networks/scripts/from_leila/log2_zscore_filter.py</i> | This script is used to filter the counts data using gaussian filtering script to obtain just the ‘expressed’ genes.


### Inter_network_reconstruction

This folder contains all scripts used to build the intercellular networks connecting differentially expressed ligands on epithelial cells to receptors expressed on immune cells.

* <i>download_lr_extra_omnipathR.R</i> | Script to download omnipath ligand-receptor interactions table used to build the intercellular networks. 

* <i>epithelial_organoid_immune_cells_interaction_{ileum|colon}.py</i> | Script used to calculate all ligand-receptor interactions for each epithelial-immune cell pairs

* <i>ligand_receptor_stats.R</i> | Script used to create different tables based on epithelial cell-immune cell connections

* <i>ligand_imm_cell_network_creation.R</i> | Script that creates tables for further intercellular networks visualization. This script has been only applied for infected/bystander immature enterocytes 2 data. 

* <i>Ligand_receptor_family_network_creation.R</i> | This script is used to create the ligand-receptor family interaction table, which is used for the overview of the intercellular network of infected immature enterocytes to immune cells.


### Inter_network_analysis

This folder contains all scripts used to create different plots used to analyse the intercellular networks between bystander/infected immature enterocytes and immune cells in both colon and ileum.

* <i>barplot_num_ligands.R</i> | Script to create bar chart of number of epithelial ligands per cell type

* <i>Cell_cell_heatplots_cellspecific.R</i> | Script to create a heat plot showing the number of interactions between epithelial cells vs immune cell.

*	<i>barplot_num_ints_per_epi_ligand.R</i> | Script to create a bar plot showing the number of interactions across all immune cells for each epithelial ligand

* <i>barplot_num_ints_per_receptor_num_imm_cells.R</i> | Script to create a bar plot of number of interactions across all epithelial cells for each immune cell receptor.

* <i>ligand_receptor_heatplots.R</i> | Script to create a heatplot showing the number of interactions between ligands and receptors.

*	<i>functional_analysis_ligand_receptors.R</i> | Script to carry out the functional analysis of ligands and receptors involved in each intercellular interactions pairs.

* <i>ligand_imm_cell_heatplots_receptor_expression.R</i> | Script to plot ligands against immune cells where the colour of the box represents the sum of receptor expression values (accounting for both the number of receptors and the level of expression)

*	<i>receptor_imm_cell_heatplots.R</i> | Script to plot receptors against immune cells with the colour of the box representing the total number of interactions.

### Inter_output_data
This folder contains ligand-receptor interaction tables between infected/bystander epithelial and immune cell populations. 


### Intra_input_data

This folder contains all the necessary epithelial and immune cell input data used for the construction of intracellular networks (using Virallink and CARNIVAL) for infected immature enterocytes only and immune cells in both colon and ileum.

### Intra_network_reconstruction

## /virallink

Scripts in this folder use the ViralLink pipeline to connect SARS-CoV-2 proteins or miRNAs (run seprarately) to differentially expressed ligands of infected immature enterocytes in ileum and colon (run separately) and creates an intracellular signalling network.

The python wrapper can be used to automatically run all scripts in this folder: virallink.py

## /carnival

* <i>Classical_pipeline_viral_{proteins/miRNAs}.R</i>
* <i>CARNIVAL.R</i>

These two scripts use three tools (Progeny, Dorothea and CARNIVAL) to connect SARS-CoV-2 proteins or miRNAs (run separately) to differentially expressed ligands of infected immature enterocytes in the ileum and colon (run separately).

### Intra_network_analysis

* <i>Functional_analysis_merged_networks.R</i> | This script runs a GO/Reactome functional overrepresentation analysis of the protein-protein interaction (PPI) layer of the final intracellular networks.

### Intra_output_data

This folder contains the node table of the intracellular signalling network for infected epithelial cells upon SARS-CoV-2 infection.




