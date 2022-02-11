# Gut-Covid pipeline

## Overview

Here, we propose an integrated systems biology workflow that models how infection alters the intracellular signalling of epithelial cells and how this change impacts the systemic immune response through modified interactions between epithelial cells and local immune cell populations. 

As a proof-of-concept, we focused on the effect of SARS-CoV-2 infection on intestinal and upper-airway epithelial cells. To characterise the modified epithelial-immune interactome, we integrated intra- and intercellular networks with single-cell RNA-seq data from SARS-CoV-2 infected human ileal and colonic organoids as well as from infected airway ciliated epithelial cells from moderate COVID-19 patients. This integrated methodology has proven useful to point out specific epithelial-immune interactions driving inflammation during disease response, and propose relevant molecular targets to guide focused experimental analysis.

This workflow is not limited to SARS-CoV-2 infection and intestinal or upper airway epithelial cells, but can be applied to any epithelial-immune populations provided the right input data is available.

The pipeline is currently available as a series of R and Python scripts which can be run through two different methods:
With docker: the whole pipeline and/or the separate stages can be run from within a docker container which negates the need for local Python and R installations making the pipeline easily accessible (recommended).
Without docker: the whole pipeline (via a Python wrapper script) and/or the separate scripts can be run locally using local installations of Python and R with associated packages.

More detailed information about this project is available in the following paper:

>Poletti M., Treveil A. et al. [Reprogramming of the intestinal epithelial-immune cell interactome during SARS-CoV-2 infection] : https://doi.org/10.1101/2021.08.09.455656, BioRxiv, (2021)

NOTE: This github repository has 2 pipelines: 1) original Gut-COVID pipeline (in the *pipeline* folder) and 2) the pipeline with lung data (in the *pipeline_lung* folder).
NOTE: They will be merged together soon and we are planning to create a pipeline, where the user can define the tissue.

<img src="https://raw.githubusercontent.com/korcsmarosgroup/gut-COVID/main/Workflow_figure.png" width="600" height="477">


## Getting started - Original Gut-COVID pipeline

You can have any information about this pipeline inside the "pipeline" folder.

The input R object file for this pipeline is too big to store on github, so there is a link here, from where you can download it: https://nbicloud-my.sharepoint.com/:f:/g/personal/boharb_nbi_ac_uk/EpyKV4t0_jpGo1y1k76dcD0B7gfZf2Fe8ifUnx17l19d8w?e=ntroFH
Password is: gutcovid2022
For the original pipeline run, you should download the *COVID19_Oct.rds* file. After you downloaded it, you have to copy it to the *gutcovid/deploy/pipeline/input_file* folder, which is originally an empty folder (important, that you have to do it before starting the bash script, so not inside the docker container!).

To run the dockerised environment for the project, you have to run the following command:
```
bash gutcovid.sh
```
Then you will get, you should see something like this:
```
root@3c172830ba15:/home/gutcovid#
```
With the following command, you can enter to the original pipeline folder:
```
cd pipeline
```
If you want to use other files, which are not inside the repository, first of all, do not close the docker container! Then, open a new terminal tab and run the following command:
```
docker cp {your_file_path} gutcovid:/home/gutcovid/pipeline{destination_path_of_your_file}
```
Or if you want to save files from the container, first of all, do not close the docker container! Then, open a new terminal tab and run the following command:
```
docker cp gutcovid:/home/gutcovid/pipeline{destination_path_of_your_file} {your_file_path}
```
Before you start to execute this original Gut-COVID pipeline, you have to register to the following website: *https://www.ibm.com/academic/topic/data-science*, then download the CPlex tool with version 12.10 (filename: cplex_entserv1210_linux-x86-64.bin).
After you downloaded it, you have to copy it into the docker container, with the following command:
```
docker cp {your_file_path} gutcovid:/home/gutcovid/
```
Then, you should give permissions to this file, with the following command (don't forget that, you should be inside the docker container, in the */home/gutcovid/* folder):
```
chmod u+x cplex_entserv1210_linux-x86-64.bin
```
After that, you have to install this tool inside the docker container, with the following command:
```
./cplex_entserv1210_linux-x86-64.bin
```
During its installation, you have to type the following keys after each other:
```
2, ENTER, 1, ENTER, ENTER, ENTER
```
After these steps, when you get back the *root@3c172830ba15:/home/gutcovid#* prompt, then you can start the pipeline with the following command:
```
python3 gutcovid.py
```

This repository contains all the input files, pipelines and generated intracellular and extracellular networks related to the following publication:

*Poletti et al. (2021) , Reprogramming of the intestinal epithelial-immune cell interactome during SARS-CoV-2 infection*

Specifically, the contained code is used to process the single cell RNAseq epithelial and immune data and use it as input to build intracellular and intercellular epithelial-immune cell interaction networks affected upon SARS-CoV-2 infection.
This repository is divided into different folders, whereby the suffix “intra” or “extra” indicates whether the section relates to the intracellular or intercellular networks. For each set of folders, there is one folder for the input data, one for the network reconstruction part and one for the network analysis part. Additionally, one folder contains scripts used to process the epithelial and immune cell data prior to network reconstruction and analysis. Please note that some manual/non-scripted analysis was carried out alongside these scripts, for example the merging of the viral_protein and viral_miRNAs networks was carried out using the merge function in Cytoscape (Institute for Systems Biology, 2019). Description of the different folders and scripts is provided below.

## 1_Input_data

**/intercellular_networks**

This folder contains all the necessary epithelial and immune cell input data used for the construction of intercellular networks for bystander/infected epithelial cells in both colon and ileum (Triana et al., 2021a; Smillie et al., 2019; Martin et al., 2019)

•	*COVID19_Oct.rda* | R object contain gene expression data for colonic and ileal epithelial cell populations of organoids infected with SARS-CoV-2 from Triana et al., 2021a

•	*imm_average_expression_uc_smillie_et_al.tsv* 
•	*imm_average_expression_cd_martin_et_al.tsv*

Gene expression table for colonic and ileal immune cell populations from Smillie et al., 2019 and Martin et al., 2019, respectively.

•	*lr_Network_Omnipath_23_09.tsv* | A priori ligand-receptor interaction table downloaded from OmniPath (Turei et al., 2016; Turei et at., 2021)

•	*cells_immune_10x_data.xlsx* | Table summarizing all immune cell populations selected for the intercellular network creation.

•	*ileum_inf_immentero_immune_receptor_categories.xlsx*
•	*colon_inf_immentero_immune_receptor_categories.xlsx*

Tables categorizing all receptors on immune cells into bigger receptor categories.

**/intracellular_networks**

This folder contains all the necessary epithelial and immune cell input data used for the construction of intracellular networks (using Virallink and CARNIVAL (Treveil et al., 2021; Liu et al., 2019) for infected immature enterocytes only and immune cells in both colon and ileum.

•	*miRNA_to_human_19.txt* | This file contains the SARS-CoV-2 viral miRNA – human protein interaction table from Saçar Demirci and Adan, 2020

•	*sarscov2-human_viral_ppis.txt* | This file contains the SARS-CoV-2 viral protein – human protein interaction table downloaded from IntAct (Hermjakob et al., 2004; Orchard et al., 2014)

•	*sarscov2_protein_annotations.txt* | Annotation table for SARS-CoV-2 proteins 

*/intracellular_networks/a_priori_networks*

•	*dorothea_abc_signed_directed_dev.txt* | This file contains the a priori transcription factor (TF) – target gene (TG) interaction table downloaded from DoroThea (Garcia-Alonso et al., 2019).

•	*omnipath_signed_directed.txt* | This file contains the a priori protein-protein interaction (PPI) table downloaded from OmniPath (Türei et al., 2016).

**/intracellular_networks/expression_data**

•	*unfilt_degs_ileum_24h_infected_enterocytes_vs_mock.csv*
•	*unfilt_degs_colon_24h_infected_enterocytes_vs_mock.csv*

Unfiltered differentially expressed genes (DEGs) tables for SARS-CoV-2 infected immature enterocytes vs mock infected cells in ileal and colonic organoids upon infection from Triana et al., 2021.

•	*filt_degs_ileum_24h_infected_enterocytes_vs_mock_padj0.05_lfc0.5.csv*
•	*filt_degs_colon_24h_infected_enterocytes_vs_mock_padj0.05_lfc0.5.csv*

Filtered differentially expressed genes (DEGs) tables (p-value < 0.05) for SARS-CoV-2 infected immature enterocytes vs mock infected cells in ileal and colonic organoids upon infection from Triana et al., 2021.

•	*expressed_genes_ileum_infected.txt*
•	*expressed_genes_colon_infected.txt*

Tables of expressed genes (post-gaussian filtering) in ileal and colonic infected immature enterocyte populations from from Triana et al., 2021.

**/functional_analysis**

•	*reactome_annotations_uniprot.txt* | Reactome biological pathways and processes annotation table related to Uniprot proteins. This is an input file for the functional overrepresentation analysis of the intracellular and intercellular networks. 


## 2_data_pre-processing

This folder contains all the scripts used to process immune cell and epithelial single cell RNAseq data for consequent construction of intracellular and extracellular networks.

•	*{ileum|colon}_counts_table_from_10x_data.R* | Script used to extract R object from Triana et al., 2021, the average expression / counts for each sample as well as carry out differential expression analysis (using MAST).

•	log2_zscore_filter.py | This script is used to filter the counts data for epithelial and immune cells using gaussian filtering to obtain just the ‘expressed’ genes.

## 3_intra_network_reconstruction

**/virallink**

Scripts in this folder use the ViralLink pipeline to connect SARS-CoV-2 proteins or miRNAs (run separately) to differentially expressed ligands of infected immature enterocytes in ileum and colon (run separately) and creates an intracellular signalling network.
The python wrapper can be used to automatically run all scripts in this folder: virallink.py
/carnival

•	*Classical_pipeline_viral_{proteins/miRNAs}.R*
•	*CARNIVAL.R*

These two scripts use three tools (Progeny (Schubert et al., 2018), VIPER (Alvarez et al., 2016), and CARNIVAL (Liu et al., 2019)) to connect SARS-CoV-2 proteins or miRNAs (run separately) to differentially expressed ligands of infected immature enterocytes in the ileum and colon (run separately).

## 4_inter_network_reconstruction

This folder contains all scripts used to build the intercellular networks connecting differentially expressed ligands on epithelial cells to receptors expressed on immune cells.

•	*download_lr_extra_omnipathR.R* | Script to download omnipath ligand-receptor interactions table used to build the intercellular networks.

•	*epithelial_organoid_immune_cells_interaction_{ileum|colon}.py* | Script used to calculate all ligand-receptor interactions for each epithelial-immune cell pairs

•	*ligand_receptor_stats.R* | Script used to create different data tables based on epithelial cell-immune cell connections.

•	*ligand_imm_cell_network_creation.R* | Script that creates tables for further intercellular networks visualization. This script has been only applied for infected/bystander immature enterocytes 2 data.

•	*Ligand_receptor_family_network_creation.R* | This script is used to create the ligand-receptor family interaction table, which is used for the overview of the intercellular network of infected immature enterocytes to immune cells.

## 5_inter_network_analysis

This folder contains all scripts used to create different plots used to analyse the intercellular networks between bystander/infected immature enterocytes and immune cells in both colon and ileum.

•	*barplot_num_ligands.R* | Script to create bar chart of number of epithelial ligands per cell type.


•	*Cell_cell_heatplots_cellspecific.R* | Script to create a heat plot showing the number of interactions between epithelial cells and immune cell.

•	*barplot_num_ints_per_epi_ligand.R* | Script to create a bar plot showing the number of interactions across all immune cells for each epithelial ligand

•	*barplot_num_ints_per_receptor_num_imm_cells.R* | Script to create a bar plot of number of interactions across all epithelial cells for each immune cell receptor.

•	*ligand_receptor_heatplots.R* | Script to create a heat plot showing the number of interactions between ligands and receptors.

•	*functional_analysis_ligand_receptors.R* | Script to carry out functional analysis of ligands and receptors involved in each intercellular interaction pairs.

•	*ligand_imm_cell_heatplots_receptor_expression.R* | Script to plot ligands against immune cells where the colour of the box represents the sum of receptor expression values (accounting for both the number of receptors and the level of expression).

•	*receptor_imm_cell_heatplots.R* | Script to plot receptors against immune cells with the colour of the box representing the total number of interactions.

## 6_inter_output_data

This folder contains ligand-receptor interaction tables between infected/bystander epithelial and immune cell populations.

## 7_intra_network_analysis

•	*Functional_analysis_merged_networks.R* | This script runs a GO/Reactome functional overrepresentation analysis of the protein-protein interaction (PPI) layer of the final intracellular networks.

## 8_intra_output_data***

This folder contains the node table of the intracellular signalling network for infected epithelial cells upon SARS-CoV-2 infection.


## Getting started #2 - Pipeline with Lung data

You can have any information about this pipeline inside the "pipeline_lung" folder

The input R object file for this pipeline is too big to store on github, so there is a link here, from where you can download it: https://nbicloud-my.sharepoint.com/:f:/g/personal/boharb_nbi_ac_uk/EpyKV4t0_jpGo1y1k76dcD0B7gfZf2Fe8ifUnx17l19d8w?e=ntroFH
Password is: gutcovid2022
For the lung data related pipeline run, you should download the *covid_nbt_main.rds* file. After you downloaded it, you have to copy it to the *gutcovid/deploy/pipeline_lung/input_file* folder, which is originally an empty folder (important, that you have to do it before starting the bash script, so not inside the docker container!).

To run the dockerised environment for the project, you have to run the following command: *bash gutcovid.sh*
Then you will get, you should see something like this:
```
root@3c172830ba15:/home/gutcovid#
```
With the following command, you can enter to this version folder of the pipeline:
```
cd pipeline_lung
```
If you want to use other files, which are not inside the repository, first of all, do not close the docker container! Then, open a new terminal tab and run the following command:
```
docker cp {your_file_path} gutcovid:/home/gutcovid/pipeline_lung/{destination_path_of_your_file}
```
Or if you want to save files from the container, first of all, do not close the docker container! Then, open a new terminal tab and run the following command:
```
docker cp gutcovid:/home/gutcovid/pipeline{destination_path_of_your_file} {your_file_path}
```
NOTE: In this pipeline with the lung data, you do not have to install the CPlex tool, you can start the pipeline immediately with the following command:
```
python3 gutcovid_lung.py
```

This repository contains all the input files for the lung-related data, pipelines and generated intracellular and extracellular networks related to the following publication:

*Poletti et al. (2021) , Reprogramming of the intestinal epithelial-immune cell interactome during SARS-CoV-2 infection*

Specifically, the contained code is used to process the single cell RNAseq epithelial and immune data and use it as input to build intracellular and intercellular epithelial-immune cell interaction networks affected upon SARS-CoV-2 infection.
This repository is divided into different folders, whereby the suffix “intra” or “extra” indicates whether the section relates to the intracellular or intercellular networks. For each set of folders, there is one folder for the input data, one for the network reconstruction part and one for the network analysis part. Additionally, one folder contains scripts used to process the epithelial and immune cell data prior to network reconstruction and analysis. Please note that some manual/non-scripted analysis was carried out alongside these scripts, for example the merging of the viral_protein and viral_miRNAs networks was carried out using the merge function in Cytoscape (Institute for Systems Biology, 2019). Description of the different folders and scripts is provided below.

## Some further information about this lung data related pipeline

# Regarding the 4_Intra_create_merged_network folder

- Here, you have to take the two output networks from virallink (*3_Intra_analysis_ViralLink* folder) (viral RNA and proteins) and you have to import them into Cytoscape, following to merge them using the *Merge* function of Cytoscape. 
- NOTE: before merging, you have to add an additional column called 'source_node' (viral protein OR viralmirna) for the node table, which in the merged network will allow to remember if the node belonged to the virallink mirna or virallink protein network.
- After merging, you have to export the node table of the merged network, which will be used to run the functional analysis(5_functional_analysis) in the next step.

# Regarding the 6_Intra_Inter_summary_figure folder

This network is the summary, where you have to connect the intracellular merged network with the intercellular one.

## References list

Alvarez, M.J., Shen, Y., Giorgi, F.M., Lachmann, A., Ding, B.B., Ye, B.H., and Califano, A. (2016). Functional characterization of somatic mutations in cancer using network-based inference of protein activity. Nat. Genet. 48, 838–847.

Garcia-Alonso, L., Holland, C.H., Ibrahim, M.M., Turei, D., and Saez-Rodriguez, J. (2019). Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Res. 29, 1363–1375.

Hermjakob, H., Montecchi-Palazzi, L., Lewington, C., Mudali, S., Kerrien, S., Orchard, S., Vingron, M., Roechert, B., Roepstorff, P., Valencia, A., et al. (2004). IntAct: an open source molecular interaction database. Nucleic Acids Res. 32, D452-5.

Institute for Systems Biology, 2019. Cytoscape, Available at: https://www.cytoscape.org.

Liu, A., Trairatphisan, P., Gjerga, E. et al. From expression footprints to causal pathways: contextualizing large signaling networks with CARNIVAL. npj Syst Biol Appl 5, 40 (2019). https://doi.org/10.1038/s41540-019-0118-z

Orchard, S., Ammari, M., Aranda, B., Breuza, L., Briganti, L., Broackes-Carter, F., Campbell, N.H., Chavali, G., Chen, C., del-Toro, N., et al. (2014). The MIntAct project - IntAct as a common curation platform for 11 molecular interaction databases. Nucleic Acids Res. 42, D358-63.

Smillie, C.S., Biton, M., Ordovas-Montanes, J., Sullivan, K.M., Burgin, G., Graham, D.B., Herbst, R.H., Rogel, N., Slyper, M., Waldman, J., et al. (2019). Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis. Cell 178, 714-730.e22.

Saçar Demirci, M.D., and Adan, A. (2020). Computational analysis of microRNA-mediated interactions in SARS-CoV-2 infection. PeerJ 8, e9369.

Schubert, M., Klinger, B., Klünemann, M., Sieber, A., Uhlitz, F., Sauer, S., Garnett, M.J., Blüthgen, N., and Saez-Rodriguez, J. (2018). Perturbation-response genes reveal signaling footprints in cancer gene expression. Nat. Commun. 9, 20.

Treveil, A., Bohar, B., Sudhakar, P., Gul, L., Csabai, L., Olbei, M., Poletti, M., Madgwick, M., Andrighetti, T., Hautefort, I., et al. (2021). ViralLink: An integrated workflow to investigate the effect of SARS-CoV-2 on intracellular signalling and regulatory pathways. PLoS Comput. Biol. 17, e1008685.

Triana, S., Metz-Zumaran, C., Ramirez, C., Kee, C., Doldan, P., Shahraz, M., Schraivogel, D., 
Gschwind, A.R., Sharma, A.K., Steinmetz, L.M., et al. (2021a). Single-cell analyses reveal SARS-CoV-2 interference with intrinsic immune response in the human gut. Mol. Syst. Biol. 17, e10232.

Triana, S., Stanifer, M.L., Boulant, S., and Alexandrov, T. (2021b). COVID19_July.rda. Figshare.

Türei, D., Korcsmáros, T., and Saez-Rodriguez, J. (2016). OmniPath: guidelines and gateway for literature-curated signaling pathway resources. Nat. Methods 13, 966–967.

Türei, D., Valdeolivas, A., Gul, L., Palacio-Escat, N., Klein, M., Ivanova, O., Ölbei, M., Gábor, A., Theis, F., Módos, D., et al. (2021). Integrated intra- and intercellular signaling knowledge for multicellular omics analysis. Mol. Syst. Biol. 17, e9923.
