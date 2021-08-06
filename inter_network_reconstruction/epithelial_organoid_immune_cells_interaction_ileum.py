# Predicting the effect of SARS-CoV2 infected epithelial cells (DEGs) on immune cells in small intestine
# Input files (each have gene symbols):
    # 1. Epithelial DEGs (healthy vs. SARS-CoV2 infected cells)
    # 2. Expressed genes in different immune cell types
        # Log2 based z-score filtration (log2_zscore_filter.py) - genes with low expression values have been discarded
        # Must be log2(tpm) values - this script converts to log2(tpm+1)
    # 3, Ligand-receptor connections from OmniPath
        # Ligand - receptor extra interactions involved as well (download_lr_extra_omnipathR.R)
# Authors: Lejla Gul, Dezso Modos, Agatha Treveil

##### Import and load #####

# Import libraries
import pandas as pd
import numpy as np
from pprint import pprint

# Importing expression information from healthy and diseased conditions in different cells as pandas dataframes
df_imm = pd.read_csv('../data/sc_data/filtered_by_zscore/CD_small_intest/martin_2019/filtered_genes_log2_CD_martin_et_al.csv',
                     index_col = 0) # contains all of the cells and conditions

#df_epi_deg = pd.read_csv('../data/infected_epithelial_cells/triana_2020/degs/filt_degs_ileum_24h_infected_v_mock_all_cells_padj0.05_lfc0.5.txt',
#                         index_col=5, sep = "\t")
df_epi_deg = pd.read_csv('../data/infected_epithelial_cells/triana_2020/degs/filt_degs_ileum_24h_bystander_v_mock_all_cells_padj0.05_lfc0.5.txt',
    index_col=5, sep="\t")

# specify ligand receptor interactions table
lr_file_p = '../data/ligand_receptor_interaction/lr_Network_Omnipath_23_09.tsv'

# specify output filepaths
out_file1 =  "../results/triana_2020/ligand_receptor_ints/3rd_set_imm_cells/ileum/bystander/ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt"
#out_file1 =  "../results/triana_2020/ligand_receptor_ints/3rd_set_imm_cells/ileum/infected/ligand_receptor_ints_triana_ileum_infec_martin_si_uninfl.txt"

# Specify which comparisons we want from the epithelial degs table
comparisons_degs = {'imm_enterocyte_2': 'Illeum_24h_Inmature Enterocyte 2_Bystander_v_Illeum_Mock_Inmature Enterocyte 2_Non-Infected',
                 'imm_enterocyte_1': 'Illeum_24h_Inmature Enterocyte 1_Bystander_v_Illeum_Mock_Inmature Enterocyte 1_Non-Infected',
                 'cycling_ta': 'Illeum_24h_Cycling TA_Bystander_v_Illeum_Mock_Cycling TA_Non-Infected',
                 'ta': 'Illeum_24h_TA_Bystander_v_Illeum_Mock_TA_Non-Infected',
                 'secretory_ta': 'Illeum_24h_Secretory TA_Bystander_v_Illeum_Mock_Secretory TA_Non-Infected',
                 'enterocyte_1': 'Illeum_24h_Enterocyte 1_Bystander_v_Illeum_Mock_Enterocyte 1_Non-Infected',
                 'stem': 'Illeum_24h_Stem Cells_Bystander_v_Illeum_Mock_Stem Cells_Non-Infected',
                 'goblet': 'Illeum_24h_Goblet Cells_Bystander_v_Illeum_Mock_Goblet Cells_Non-Infected',
                 'enteroendocrine': 'Illeum_24h_Enteroendocrine cells_Bystander_v_Illeum_Mock_Enteroendocrine cells_Non-Infected'}

# comparisons_degs = {'imm_enterocytes_2': 'Illeum_24h_Inmature Enterocyte 2_Infected_v_Illeum_Mock_Inmature Enterocyte 2_Non-Infected',
#                 'imm_enterocytes_1': 'Illeum_24h_Inmature Enterocyte 1_Infected_v_Illeum_Mock_Inmature Enterocyte 1_Non-Infected',
#                 'cycling_ta': 'Illeum_24h_Cycling TA_Infected_v_Illeum_Mock_Cycling TA_Non-Infected',
#                 'ta': 'Illeum_24h_TA_Infected_v_Illeum_Mock_TA_Non-Infected',
#                 'enterocytes_1': 'Illeum_24h_Enterocyte 1_Infected_v_Illeum_Mock_Enterocyte 1_Non-Infected'}

# Specify which cells we want from the immune expression table
imm_cells = {'DC1': 'RNA.DC1_uninflamed', 'DC2': 'RNA.DC2_uninflamed',
             'Treg': 'RNA.Tregs_uninflamed', 'mast_cell': 'RNA.Mast.cells_uninflamed',
             'ILC1': 'RNA.ILC1_uninflamed','ILC3': 'RNA.ILC3_uninflamed',
             'resident_macrophage': 'RNA.Resident.macrophages_uninflamed','naive_Tcell': 'RNA.Naive.T.cells_uninflamed',
             'inf_macrophage': 'RNA.Inf..Macrophages_uninflamed',
             'activated_Tcell': 'RNA.Highly.activated.T.cells_uninflamed',
             'Tcyto': 'RNA.Cytotoxic.T.cells_uninflamed','Tmem': 'RNA.Central.Memory.T.cells_uninflamed',
             'Trm':'RNA.Trm_uninflamed', 'Trm_Th17':'RNA.Th17.Trm_uninflamed',
             'CD8_Trm':'RNA.CD8.Trm_uninflamed','CD8_Trm_cyto':'RNA.Cytotoxic.CD8.Trm_uninflamed',
             'gd_Tcell':'RNA.T..gd._uninflamed',
             'mem_Bcell':'RNA.Memory.B.cells_uninflamed', 'naive_Bcell':'RNA.Naive.B.cells_inflamed',
             'TFH_like':'RNA.TFH.like_uninflamed','plasmablast':'RNA.Plasmablasts_uninflamed',
             'IgA_plasma':'RNA.IgA.plasma.cells_uninflamed','IgG_plasma':'RNA.IgG.plasma.cells_uninflamed',
             'IgM_plasma':'RNA.IgM.plasma.cells_uninflamed'}

##### Set up ligand receptor dictionaries #####

# create dictionaries for source-target interactions and their annotations from OmniPath
interactions = {}
interaction_annotation = {}

# OmniPathR version: Bioconductor version 3.12 - developer version (https://bioconductor.org/packages/devel/bioc/html/OmnipathR.html)
# Ligand - receptor extras
with open (lr_file_p) as interaction:
    interaction.readline()
    for line in interaction:
        line = line.strip().split('\t')
        # For each ligand, add as key to dictionary if not already there
        if line[0] not in interactions:
            interactions[line[0]] = set()
        # Add receptor for this ligand as one value in a set (whole set is the dictionary value). Sets can't have duplicates
        interactions[line[0]].add(line[1])
        # create a similar dictionary (single ligand -> multiple receptors), where ligand and receptor written (key and value of dict are tupules)
        if (line[0],'Ligand') not in interaction_annotation:
            interaction_annotation[(line[0], 'Ligand')] = set()
        interaction_annotation[(line[0],'Ligand')].add((line[1], 'Receptor'))

#print(interactions)

##### Set up epithelial cell DEGs #####

# dict comprehension to get take the comparisons_deg dict, keep the key static but use the value to filter the comparison col of the degs table and extract just the lfc column
epi_degs = {tup[0]: df_epi_deg[df_epi_deg.comparison.eq(tup[1])]['avg_logFC'] for tup in comparisons_degs.items()}

##### Set up immune cell expression #####

# Convert expression to log2(tpm+1)
df_imm = df_imm.replace('NaN',np.nan)
df_imm = df_imm.transform(lambda x: np.log2((2**x)+1))

# selecting uninflamed average expression data (log2 value, after z-score filtration) in the selected cell types in
# dictionary to store data frames and make it callable in a for cycle
#'celltype':df_imm[name_of_the_column]
imm_cell_types = {tup[0]: df_imm[tup[1]] for tup in imm_cells.items()}

##### Process #####

# Cell-cell interactions to process
# Epithelial cell DEGs - uninflamed CD immune cell pairs - change the list if you have different immune cells
cell_cell_interactions = [(x,y) for x in epi_degs.keys() for y in imm_cell_types.keys()]

#print(epi_degs['imm_enterocytes_2_degs'])#.keys())
cols = ["ligand","ligand_lfc","receptor","receptor_expression","ligand_cell_source","receptor_cell_target"]
epi_imm_interactions = pd.DataFrame(columns=cols)

# going through each cell pair and select interactions between them where one part is represented in one condition only
for i in cell_cell_interactions:

    # iterating through source cell type genes (here keys are all genes which have lfc values for the specified cell type in the epi_DEGs dictionary)
    for source in epi_degs[i[0]].keys():
        # selecting those ones which play role in intercellular communication as a transmitter
        if source in interactions.keys():
       # If gene is a ligand then get all the receptors
            for target in interactions[source]:
                # selecting those ones which play role in intercellular communication as a receiver
                if target in imm_cell_types[i[1]].keys():
                    #print(uninflamed_cell_types[i[1]][target])
                    if np.isnan(imm_cell_types[i[1]][target]) == False:
                        # Add to set - ligand name, ligand lfc, target name, target expression level, source cell type, target cell type
                        row = pd.Series([source, epi_degs[i[0]][source], target, imm_cell_types[i[1]][target], i[0], i[1]], index=cols)
                        epi_imm_interactions = epi_imm_interactions.append(row, ignore_index=True)
                            #add((source, epi_DEGs[i[0]][source], target, uninflamed_cell_types[i[1]][target], i[0],i[1]))

# Sort dataframe for cell type - remove replicates
epi_imm_interactions = epi_imm_interactions.sort_values(by=["ligand_cell_source","receptor_cell_target","ligand_lfc"])

# writing out all interactions in all conditions
epi_imm_interactions.to_csv(out_file1, sep='\t', index=False)
