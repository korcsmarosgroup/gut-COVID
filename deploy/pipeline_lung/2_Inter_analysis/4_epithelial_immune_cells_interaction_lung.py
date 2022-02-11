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
import argparse
import sys

def parse_args(args):
    help_text = \
        """
        === Epithelial immune cells interaction lung ===.
        """

    parser = argparse.ArgumentParser(description=help_text)

    parser.add_argument("-imm", "--input-file-imm",
                        help="<path to the input file> [mandatory]",
                        type=str,
                        dest="input_file_imm",
                        action="store",
                        required=True)

    parser.add_argument("-epi", "--input-file-epi",
                        help="<path to the input file> [mandatory]",
                        type=str,
                        dest="input_file_epi",
                        action="store",
                        required=True)

    parser.add_argument("-omni", "--omnipath-file",
                        help="<path to the input file> [mandatory]",
                        type=str,
                        dest="omnipath_file",
                        action="store",
                        required=True)

    parser.add_argument("-of", "--output-file",
                        help="<path to the output folder> [mandatory]",
                        type=str,
                        dest="output_file",
                        action="store",
                        required=True)

    results = parser.parse_args(args)
    return results.input_file_imm, results.input_file_epi, results.omnipath_file, results.output_file


def main():

    input_file_imm, input_file_epi, omnipath_file, output_file = parse_args(sys.argv[1:])

    # Importing expression information from healthy and diseased conditions in different cells as pandas dataframes
    df_imm = pd.read_csv(input_file_imm, index_col = 0) # contains all of the cells and conditions

    df_epi_deg = pd.read_csv(input_file_epi, index_col=0, sep="\t")  

    # specify ligand receptor interactions table
    lr_file_p = omnipath_file

    # specify output filepaths (1= all interactions)
    out_file1 =  output_file

    # Specify which comparisons we want from the epithelial degs table
    comparisons_degs = {'ciliated_cells': 'ciliated_sars_vs_healthy'}

    # Specify which cells we want from the immune expression table
    imm_cells = {'Bcell': 'B.cell_control',
                'Tcell_cyto': 'CTL_control',	
                'Mastcell': 'MC_control',
                'MD_macrophage': 'MoD.Ma_control',
                'DC_MD': 'moDC_control',	
                'NKcell': 'NKT_control',
                'NKcell_prolif': 'NKT.p_control',
                'Nonresident_macrophage': 'nrMa_control',	
                'DC_plasma': 'pDC_control',
                'Resident_macrophage': 'rMa_control',	
                'Treg': 'Treg_control'}

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
    # Epithelial cell DEGs - change the list if you have different immune cells
    cell_cell_interactions = [(x,y) for x in epi_degs.keys() for y in imm_cell_types.keys()]

    #print(epi_DEGs['imm_enterocytes_2_degs'])#.keys())
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
                        #if str(uninflamed_cell_types[i[1]][target]) != 'nan':
                        if np.isnan(imm_cell_types[i[1]][target]) == False:
                            #print(np.isnan(uninflamed_cell_types[i[1]][target]))
                            # Add to set - ligand name, ligand lfc, target name, target expression level, source cell type, target cell type
                            row = pd.Series([source, epi_degs[i[0]][source], target, imm_cell_types[i[1]][target], i[0], i[1]], index=cols)
                            epi_imm_interactions = epi_imm_interactions.append(row, ignore_index=True)
                                #add((source, epi_DEGs[i[0]][source], target, uninflamed_cell_types[i[1]][target], i[0],i[1]))

    # Sort dataframe for cell type - remove replicates
    epi_imm_interactions = epi_imm_interactions.sort_values(by=["ligand_cell_source","receptor_cell_target","ligand_lfc"])

    # writing out all interactions in all conditions
    epi_imm_interactions.to_csv(out_file1, sep='\t', index=False)

if __name__ == "__main__":
    main()
