# Main python script to run the whole Gut-COVID pipeline
from distutils.dir_util import copy_tree
from time import strftime
import subprocess
import sys
import os
from os import path
import time


scripts_folders = {
    "2_data_pre_processing": ["ileum_counts_table_from_10x_data.R",
                                "colon_counts_table_from_10x_data.R",
                                "log2_zscore_filter.py"],
    "3_intra_network_reconstruction": ["carnival/Classical_pipeline_viral_proteins.R",
                                        "carnival/Classical_pipeline_viral_miRNAs.R",
                                        "carnival/CARNIVAL.R",
                                        "virallink/filter_network_expressed_genes.R",
                                        "virallink/get_regulator_deg_network.R",
                                        "virallink/prepare_tiedie_input.R",
                                        "virallink/TieDie/tiedie.py",
                                        "virallink/combined_edge_node_tables.R"],
    "4_inter_network_reconstruction": ["download_lr_extra_omnipathR.R",
                                   "epithelial_organoid_immune_cells_interaction_colon.py",
                                   "epithelial_organoid_immune_cells_interaction_ileum.py",
                                   "ligand_receptor_stats.R",
                                   "ligand_imm_cell_network_creation.R",
                                   "ligand_receptor_family_network_creation.R"],
    "5_inter_network_analysis": ["get_num_ligands_per_network.R",
                                "get_num_ligands_per_network_up_down.R",
                                "barplot_num_ligands.R",
                                "cell_cell_heatplots_cellspecific.R",
                                "barplot_num_ints_per_epi_ligand.R",
                                "barplot_num_ints_per_receptor_num_imm_cells.R",
                                "ligand_receptor_heatplots.R",
                                "functional_analysis_ligand_receptors.R",
                                "ligand_imm_cell_heatplots_receptor_expression.R",
                                "receptor_imm_cell_heatplots.R"],
    "7_intra_network_analysis": ["Functional_analysis_merged_networks.R"]
}

scripts_parameters = {
    "ileum_counts_table_from_10x_data.R": [["R_object_input_data"]],
    "colon_counts_table_from_10x_data.R": [["R_object_input_data"]],
    "log2_zscore_filter.py": [["-i", "log2_zscore_filtering_cd_martin_input",
                                "-o", "/home/gutcovid/pipeline/1_input_data/intercellular_networks/filtered_genes_log2_CD_martin_et_al.csv"],
                                ["-i", "log2_zscore_filtering_uc_smillie_input",
                                "-o", "/home/gutcovid/pipeline/1_input_data/intercellular_networks/filtered_imm_log2_UC_smillie_et_al.csv"]],
    "Classical_pipeline_viral_miRNAs.R": [["classic_pipeline_colon",
                                            "dorothea_input_file_v2",
                                            "omnipath_input_file",
                                            "miRNA_input_file",
                                            "/home/gutcovid/pipeline/8_intra_output_data/carnival/colon_mirnas_results_top50tf_pleitropic_minsize15.rds",
                                            "/home/gutcovid/pipeline/8_intra_output_data/carnival/colon_mirnas_results_top50tf_pleitropic_minsize15"],
                                            ["classic_pipeline_ileum",
                                            "dorothea_input_file_v2",
                                            "omnipath_input_file",
                                            "miRNA_input_file",
                                            "/home/gutcovid/pipeline/8_intra_output_data/carnival/ileum_mirnas_results_top50tf_pleitropic_minsize15.rds",
                                            "/home/gutcovid/pipeline/8_intra_output_data/carnival/ileum_mirnas_results_top50tf_pleitropic_minsize15"],
                                            ],
    "Classical_pipeline_viral_proteins.R": [["classic_pipeline_colon",
                                            "dorothea_input_file_v2",
                                            "omnipath_input_file",
                                            "human_viral_input_file",
                                            "/home/gutcovid/pipeline/8_intra_output_data/carnival/colon_proteins_results_top50tf_pleitropic_minsize15.rds",
                                            "/home/gutcovid/pipeline/8_intra_output_data/carnival/colon_proteins_results_top50tf_pleitropic_minsize15"],
                                            ["classic_pipeline_ileum",
                                            "dorothea_input_file_v2",
                                            "omnipath_input_file",
                                            "human_viral_input_file",
                                            "/home/gutcovid/pipeline/8_intra_output_data/carnival/ileum_proteins_results_top50tf_pleitropic_minsize15.rds",
                                            "/home/gutcovid/pipeline/8_intra_output_data/carnival/ileum_proteins_results_top50tf_pleitropic_minsize15"]],
    "CARNIVAL.R": [["/home/gutcovid/pipeline/8_intra_output_data/carnival/colon_mirnas_results_top50tf_pleitropic_minsize15.rds"],
                    ["/home/gutcovid/pipeline/8_intra_output_data/carnival/ileum_mirnas_results_top50tf_pleitropic_minsize15.rds"],
                    ["/home/gutcovid/pipeline/8_intra_output_data/carnival/colon_proteins_results_top50tf_pleitropic_minsize15.rds"],
                    ["/home/gutcovid/pipeline/8_intra_output_data/carnival/ileum_proteins_results_top50tf_pleitropic_minsize15.rds"]],
    "filter_network_expressed_genes.R": [["ileum_expressed_genes_file",
                                            "dorothea_input_file",
                                            "omnipath_input_file",
                                            "id_type",
                                            "outdir_ileum_mirna"],
                                            ["ileum_expressed_genes_file",
                                            "dorothea_input_file",
                                            "omnipath_input_file",
                                            "id_type",
                                            "outdir_ileum_protein"],
                                            ["colon_expressed_genes_file",
                                            "dorothea_input_file",
                                            "omnipath_input_file",
                                            "id_type",
                                            "outdir_colon_mirna"],
                                            ["colon_expressed_genes_file",
                                            "dorothea_input_file",
                                            "omnipath_input_file",
                                            "id_type",
                                            "outdir_colon_protein"]],
    "get_regulator_deg_network.R": [["/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_mirna/2_process_a_priori_networks/dorothea_contextualised_network.txt",
                                    "/home/gutcovid/pipeline/1_input_data/intracellular_networks/expression_data/filt_degs_ileum_24h_infected_enterocytes_vs_mock_padj0.05_lfc0.5.csv",
                                    "id_type",
                                    "outdir_ileum_mirna"],
                                    ["/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_protein/2_process_a_priori_networks/dorothea_contextualised_network.txt",
                                    "/home/gutcovid/pipeline/1_input_data/intracellular_networks/expression_data/filt_degs_ileum_24h_infected_enterocytes_vs_mock_padj0.05_lfc0.5.csv",
                                    "id_type",
                                    "outdir_ileum_protein"],
                                    ["/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_mirna/2_process_a_priori_networks/dorothea_contextualised_network.txt",
                                    "/home/gutcovid/pipeline/1_input_data/intracellular_networks/expression_data/filt_degs_colon_24h_infected_enterocytes_vs_mock_padj0.05_lfc0.5.csv",
                                    "id_type",
                                    "outdir_colon_mirna"],
                                    ["/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_protein/2_process_a_priori_networks/dorothea_contextualised_network.txt",
                                    "/home/gutcovid/pipeline/1_input_data/intracellular_networks/expression_data/filt_degs_colon_24h_infected_enterocytes_vs_mock_padj0.05_lfc0.5.csv",
                                    "id_type",
                                    "outdir_colon_protein"]],
    "prepare_tiedie_input.R": [["/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_mirna/2_process_a_priori_networks/dorothea_contextualised_network.txt",
                                "/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_mirna/2_process_a_priori_networks/contextualised_regulator-deg_network.txt",
                                "/home/gutcovid/pipeline/1_input_data/intracellular_networks/expression_data/filt_degs_ileum_24h_infected_enterocytes_vs_mock_padj0.05_lfc0.5.csv",
                                "hbps_mirna",
                                "outdir_ileum_mirna"],
                                ["/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_protein/2_process_a_priori_networks/dorothea_contextualised_network.txt",
                                "/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_protein/2_process_a_priori_networks/contextualised_regulator-deg_network.txt",
                                "/home/gutcovid/pipeline/1_input_data/intracellular_networks/expression_data/filt_degs_ileum_24h_infected_enterocytes_vs_mock_padj0.05_lfc0.5.csv",
                                "hbps_protein",
                                "outdir_ileum_protein"],
                                ["/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_mirna/2_process_a_priori_networks/dorothea_contextualised_network.txt",
                                "/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_mirna/2_process_a_priori_networks/contextualised_regulator-deg_network.txt",
                                "/home/gutcovid/pipeline/1_input_data/intracellular_networks/expression_data/filt_degs_colon_24h_infected_enterocytes_vs_mock_padj0.05_lfc0.5.csv",
                                "hbps_mirna",
                                "outdir_colon_mirna"],
                                ["/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_protein/2_process_a_priori_networks/dorothea_contextualised_network.txt",
                                "/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_protein/2_process_a_priori_networks/contextualised_regulator-deg_network.txt",
                                "/home/gutcovid/pipeline/1_input_data/intracellular_networks/expression_data/filt_degs_colon_24h_infected_enterocytes_vs_mock_padj0.05_lfc0.5.csv",
                                "hbps_protein",
                                "outdir_colon_protein"]],
    "tiedie.py": [["-u", "/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_mirna/3_network_diffusion/input_files/upstream.input",
                    "-d", "/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_mirna/3_network_diffusion/input_files/downstream.input",
                    "-n", "/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_mirna/3_network_diffusion/input_files/pathway.sif",
                    "-o", "/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_mirna/3_network_diffusion/TieDIE"],
                    ["-u", "/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_mirna/3_network_diffusion/input_files/upstream.input",
                    "-d", "/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_mirna/3_network_diffusion/input_files/downstream.input",
                    "-n", "/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_mirna/3_network_diffusion/input_files/pathway.sif",
                    "-o", "/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_mirna/3_network_diffusion/TieDIE"],
                    ["-u", "/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_protein/3_network_diffusion/input_files/upstream.input",
                    "-d", "/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_protein/3_network_diffusion/input_files/downstream.input",
                    "-n", "/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_protein/3_network_diffusion/input_files/pathway.sif",
                    "-o", "/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_protein/3_network_diffusion/TieDIE"],
                    ["-u", "/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_protein/3_network_diffusion/input_files/upstream.input",
                    "-d", "/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_protein/3_network_diffusion/input_files/downstream.input",
                    "-n", "/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_protein/3_network_diffusion/input_files/pathway.sif",
                    "-o", "/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_protein/3_network_diffusion/TieDIE"]],
    "combined_edge_node_tables.R": [["/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_mirna/3_network_diffusion/TieDIE/tiedie.cn.sif",
                                    "/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_mirna/3_network_diffusion/TieDIE/heats.NA",
                                    "hbps_mirna",
                                    "sars",
                                    "/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_mirna/2_process_a_priori_networks/contextualised_regulator-deg_network.txt",
                                    "/home/gutcovid/pipeline/1_input_data/intracellular_networks/expression_data/unfilt_degs_ileum_24h_infected_enterocytes_vs_mock.csv",
                                    "id_type",
                                    "outdir_ileum_mirna"],
                                    ["/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_protein/3_network_diffusion/TieDIE/tiedie.cn.sif",
                                    "/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_protein/3_network_diffusion/TieDIE/heats.NA",
                                    "hbps_protein",
                                    "sars",
                                    "/home/gutcovid/pipeline/8_intra_output_data/virallink/ileum_protein/2_process_a_priori_networks/contextualised_regulator-deg_network.txt",
                                    "/home/gutcovid/pipeline/1_input_data/intracellular_networks/expression_data/unfilt_degs_ileum_24h_infected_enterocytes_vs_mock.csv",
                                    "id_type",
                                    "outdir_ileum_mirna"],
                                    ["/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_mirna/3_network_diffusion/TieDIE/tiedie.cn.sif",
                                    "/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_mirna/3_network_diffusion/TieDIE/heats.NA",
                                    "hbps_mirna",
                                    "sars",
                                    "/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_mirna/2_process_a_priori_networks/contextualised_regulator-deg_network.txt",
                                    "/home/gutcovid/pipeline/1_input_data/intracellular_networks/expression_data/unfilt_degs_colon_24h_infected_enterocytes_vs_mock.csv",
                                    "id_type",
                                    "outdir_ileum_mirna"],
                                    ["/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_protein/3_network_diffusion/TieDIE/tiedie.cn.sif",
                                    "/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_protein/3_network_diffusion/TieDIE/heats.NA",
                                    "hbps_protein",
                                    "sars",
                                    "/home/gutcovid/pipeline/8_intra_output_data/virallink/colon_protein/2_process_a_priori_networks/contextualised_regulator-deg_network.txt",
                                    "/home/gutcovid/pipeline/1_input_data/intracellular_networks/expression_data/unfilt_degs_colon_24h_infected_enterocytes_vs_mock.csv",
                                    "id_type",
                                    "outdir_ileum_mirna"]],
    "download_lr_extra_omnipathR.R": [["/home/gutcovid/pipeline/1_input_data/intercellular_networks/lr_Network_Omnipath.tsv",
                                        "/home/gutcovid/pipeline/1_input_data/intercellular_networks/lr_Network_Omnipath.rds"]],
    "epithelial_organoid_immune_cells_interaction_colon.py": [["-dfimm", "/home/gutcovid/pipeline/1_input_data/intercellular_networks/filtered_imm_log2_UC_smillie_et_al.csv",
                                                                "-dfepi", "/home/gutcovid/pipeline/6_inter_output_data/degs/filt_degs_colon_24h_bystander_v_mock_all_cells_padj0.05_lfc0.5.txt",
                                                                "-lr", "/home/gutcovid/pipeline/1_input_data/intercellular_networks/lr_Network_Omnipath.tsv",
                                                                "-o", "/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt"],
                                                                ["-dfimm", "/home/gutcovid/pipeline/1_input_data/intercellular_networks/filtered_imm_log2_UC_smillie_et_al.csv",
                                                                "-dfepi", "/home/gutcovid/pipeline/6_inter_output_data/degs/filt_degs_colon_24h_infected_v_mock_all_cells_padj0.05_lfc0.5.txt",
                                                                "-lr", "/home/gutcovid/pipeline/1_input_data/intercellular_networks/lr_Network_Omnipath.tsv",
                                                                "-o", "/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt"]],
    "epithelial_organoid_immune_cells_interaction_ileum.py": [["-dfimm", "/home/gutcovid/pipeline/1_input_data/intercellular_networks/filtered_genes_log2_CD_martin_et_al.csv",
                                                                "-dfepi", "/home/gutcovid/pipeline/6_inter_output_data/degs/filt_degs_ileum_24h_bystander_v_mock_all_cells_padj0.05_lfc0.5.txt",
                                                                "-lr", "/home/gutcovid/pipeline/1_input_data/intercellular_networks/lr_Network_Omnipath.tsv",
                                                                "-o", "/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt"],
                                                                ["-dfimm", "/home/gutcovid/pipeline/1_input_data/intercellular_networks/filtered_genes_log2_CD_martin_et_al.csv",
                                                                "-dfepi", "/home/gutcovid/pipeline/6_inter_output_data/degs/filt_degs_ileum_24h_infected_v_mock_all_cells_padj0.05_lfc0.5.txt",
                                                                "-lr", "/home/gutcovid/pipeline/1_input_data/intercellular_networks/lr_Network_Omnipath.tsv",
                                                                "-o", "/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt"]],
    "ligand_receptor_stats.R": [["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt",
                                "/home/gutcovid/pipeline/6_inter_output_data/expression/stats_tog_ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt",
                                "/home/gutcovid/pipeline/6_inter_output_data/expression/stats_sep_ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt",
                                "/home/gutcovid/pipeline/6_inter_output_data/expression/stats2_ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt"],
                                ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt",
                                "/home/gutcovid/pipeline/6_inter_output_data/expression/stats_tog_ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt",
                                "/home/gutcovid/pipeline/6_inter_output_data/expression/stats_sep_ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt",
                                "/home/gutcovid/pipeline/6_inter_output_data/expression/stats2_ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt"],
                                ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt",
                                "/home/gutcovid/pipeline/6_inter_output_data/expression/stats_tog_ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt",
                                "/home/gutcovid/pipeline/6_inter_output_data/expression/stats_sep_ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt",
                                "/home/gutcovid/pipeline/6_inter_output_data/expression/stats2_ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt"],
                                ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt",
                                "/home/gutcovid/pipeline/6_inter_output_data/expression/stats_tog_ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt",
                                "/home/gutcovid/pipeline/6_inter_output_data/expression/stats_sep_ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt",
                                "/home/gutcovid/pipeline/6_inter_output_data/expression/stats2_ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt"]],
    "ligand_imm_cell_network_creation.R": [["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt",
                                            "/home/gutcovid/pipeline/6_inter_output_data/colon_bystand_imm_entero2_ligand_imm_cell_network.txt"],
                                            ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt",
                                            "/home/gutcovid/pipeline/6_inter_output_data/colon_infected_imm_entero2_ligand_imm_cell_network.txt"],
                                            ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt",
                                            "/home/gutcovid/pipeline/6_inter_output_data/ileum_bystand_imm_entero2_ligand_imm_cell_network.txt"],
                                            ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt",
                                            "/home/gutcovid/pipeline/6_inter_output_data/ileum_infected_imm_entero2_ligand_imm_cell_network.txt"]],
    "ligand_receptor_family_network_creation.R": [["/home/gutcovid/pipeline/6_inter_output_data/expression/stats2_ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt",
                                                    "receptor_categories_colon",
                                                    "/home/gutcovid/pipeline/6_inter_output_data/colon_infected_ligand_receptor_category_network.txt"],
                                                    ["/home/gutcovid/pipeline/6_inter_output_data/expression/stats2_ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt",
                                                    "receptor_categories_ileum",
                                                    "/home/gutcovid/pipeline/6_inter_output_data/ileum_infected_ligand_receptor_category_network.txt"]],
    "get_num_ligands_per_network.R": [["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt",
                                        "/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt",
                                        "/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt",
                                        "/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt",
                                        "/home/gutcovid/pipeline/6_inter_output_data/number_ligands_per_network.txt"]],
    "get_num_ligands_per_network_up_down.R": [["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt",
                                                "/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt",
                                                "/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt",
                                                "/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt",
                                                "/home/gutcovid/pipeline/6_inter_output_data/number_ligands_per_network_up_down.txt"]],
    "barplot_num_ligands.R": [["/home/gutcovid/pipeline/6_inter_output_data/number_ligands_per_network_up_down.txt",
                                "/home/gutcovid/pipeline/6_inter_output_data/number_ligands_per_network_up_down_sep.png"]],
    "cell_cell_heatplots_cellspecific.R": [["/home/gutcovid/pipeline/6_inter_output_data/expression/stats_sep_ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt",
                                            "/home/gutcovid/pipeline/6_inter_output_data/expression/stats_sep_ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt",
                                            "/home/gutcovid/pipeline/6_inter_output_data/ileum_immentero2_heatplot_log2_num_interactions.png",
                                            "ileum"],
                                            ["/home/gutcovid/pipeline/6_inter_output_data/expression/stats_sep_ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt",
                                            "/home/gutcovid/pipeline/6_inter_output_data/expression/stats_sep_ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt",
                                            "/home/gutcovid/pipeline/6_inter_output_data/colon_immentero2_heatplot_log2_num_interactions.png",
                                            "colon"]],
    "barplot_num_ints_per_epi_ligand.R": [["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt",
                                            "/home/gutcovid/pipeline/6_inter_output_data/plots/ileum_infected_"],
                                            ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt",
                                            "/home/gutcovid/pipeline/6_inter_output_data/plots/ileum_bystander_"],
                                            ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt",
                                            "/home/gutcovid/pipeline/6_inter_output_data/plots/colon_infected_"],
                                            ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt",
                                            "/home/gutcovid/pipeline/6_inter_output_data/plots/colon_bystander_"]],
    "barplot_num_ints_per_receptor_num_imm_cells.R": [["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt",
                                                        "/home/gutcovid/pipeline/6_inter_output_data/plots/ileum_infected_"],
                                                        ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt",
                                                        "/home/gutcovid/pipeline/6_inter_output_data/plots/ileum_bystander_"],
                                                        ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt",
                                                        "/home/gutcovid/pipeline/6_inter_output_data/plots/colon_infected_"],
                                                        ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt",
                                                        "/home/gutcovid/pipeline/6_inter_output_data/plots/colon_bystander_"]],
    "ligand_receptor_heatplots.R": [["/home/gutcovid/pipeline/6_inter_output_data/expression/stats2_ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt",
                                    "/home/gutcovid/pipeline/6_inter_output_data/expression/stats2_ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt",
                                    "/home/gutcovid/pipeline/6_inter_output_data/plots/ileum_infected_",
                                    "/home/gutcovid/pipeline/6_inter_output_data/plots/ileum_bystand_"],
                                    ["/home/gutcovid/pipeline/6_inter_output_data/expression/stats2_ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt",
                                    "/home/gutcovid/pipeline/6_inter_output_data/expression/stats2_ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt",
                                    "/home/gutcovid/pipeline/6_inter_output_data/plots/colon_infected_",
                                    "/home/gutcovid/pipeline/6_inter_output_data/plots/colon_bystand_"]],
    "functional_analysis_ligand_receptors.R": [["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt",
                                                "/home/gutcovid/pipeline/1_input_data/intercellular_networks/lr_Network_Omnipath.tsv",
                                                "/home/gutcovid/pipeline/6_inter_output_data/functional_analysis"],
                                                ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt",
                                                "/home/gutcovid/pipeline/1_input_data/intercellular_networks/lr_Network_Omnipath.tsv",
                                                "/home/gutcovid/pipeline/6_inter_output_data/functional_analysis"],
                                                ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt",
                                                "/home/gutcovid/pipeline/1_input_data/intercellular_networks/lr_Network_Omnipath.tsv",
                                                "/home/gutcovid/pipeline/6_inter_output_data/functional_analysis"],
                                                ["/home/gutcovid/pipeline/6_inter_output_data/expression/ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt",
                                                "/home/gutcovid/pipeline/1_input_data/intercellular_networks/lr_Network_Omnipath.tsv",
                                                "/home/gutcovid/pipeline/6_inter_output_data/functional_analysis"]],
    "ligand_imm_cell_heatplots_receptor_expression.R": [["/home/gutcovid/pipeline/6_inter_output_data/ileum_infected_imm_entero2_ligand_imm_cell_network.txt",
                                                        "/home/gutcovid/pipeline/6_inter_output_data/plots/ileum_infected_"],
                                                        ["/home/gutcovid/pipeline/6_inter_output_data/colon_infected_imm_entero2_ligand_imm_cell_network.txt",
                                                        "/home/gutcovid/pipeline/6_inter_output_data/plots/colon_infected_"]],
    "receptor_imm_cell_heatplots.R": [["/home/gutcovid/pipeline/6_inter_output_data/expression/stats2_ligand_receptor_ints_triana_ileum_infected_martin_si_uninfl.txt",
                                        "/home/gutcovid/pipeline/6_inter_output_data/expression/stats2_ligand_receptor_ints_triana_ileum_bystand_martin_si_uninfl.txt",
                                        "/home/gutcovid/pipeline/6_inter_output_data/plots/ileum_infected_normal_",
                                        "/home/gutcovid/pipeline/6_inter_output_data/plots/ileum_bystand_normal_",
                                        "/home/gutcovid/pipeline/6_inter_output_data/plots/ileum_infected_smaller_",
                                        "/home/gutcovid/pipeline/6_inter_output_data/plots/ileum_bystand_smaller_"],
                                        ["/home/gutcovid/pipeline/6_inter_output_data/expression/stats2_ligand_receptor_ints_triana_colon_infected_smillie_colon_health.txt",
                                        "/home/gutcovid/pipeline/6_inter_output_data/expression/stats2_ligand_receptor_ints_triana_colon_bystand_smillie_colon_health.txt",
                                        "/home/gutcovid/pipeline/6_inter_output_data/plots/colon_infected_normal_",
                                        "/home/gutcovid/pipeline/6_inter_output_data/plots/colon_bystand_normal_",
                                        "/home/gutcovid/pipeline/6_inter_output_data/plots/colon_infected_smaller_",
                                        "/home/gutcovid/pipeline/6_inter_output_data/plots/colon_bystand_smaller_"]],
    "Functional_analysis_merged_networks.R": [["omnipath_input_file",
                                                "/home/gutcovid/pipeline/1_input_data/intracellular_networks/colon_viral_mirna_protein_merged_node_table.csv",
                                                "/home/gutcovid/pipeline/8_intra_output_data/colon/merged_shared_only_GO.txt",
                                                "/home/gutcovid/pipeline/8_intra_output_data/colon/merged_protein_only_REACT.txt"],
                                                ["omnipath_input_file",
                                                "/home/gutcovid/pipeline/1_input_data/intracellular_networks/ileum_viral_mirna_protein_merged_node_table.csv",
                                                "/home/gutcovid/pipeline/8_intra_output_data/ileum/merged_shared_only_GO.txt",
                                                "/home/gutcovid/pipeline/8_intra_output_data/ileum/merged_protein_only_REACT.txt"]]
}

def checking_input_parameters(script_parameters):
    """
    Checking the parameters
    """
    if not os.path.isfile("parameters.yml"):
        sys.stdout.write(f" WARNING: There is no appropriate parameter file! It should be parameters.yml\n")
        sys.exit(1)

    neccessary_parameters = ["R_object_input_data", "log2_zscore_filtering_cd_martin_input", "log2_zscore_filtering_uc_smillie_input", 
                            "miRNA_input_file", "human_viral_input_file", "dorothea_input_file", "dorothea_input_file_v2", "omnipath_input_file", "classic_pipeline_colon", 
                            "classic_pipeline_ileum", "receptor_categories_colon", "receptor_categories_ileum", "ileum_expressed_genes_file", 
                            "colon_expressed_genes_file", "id_type", "hbps_mirna", "hbps_protein", "sars", "outdir_ileum_mirna", "outdir_ileum_protein",
                            "outdir_colon_mirna", "outdir_colon_protein"]
    for nec_param in neccessary_parameters:
        if nec_param not in script_parameters:
            sys.stdout.write(f" WARNING: A parameter is missing from the parameters.yml file: {nec_param}\n\n")
            sys.exit(2)


def get_parameters():
    """
    Get all of the parameters from the parameter.yml file
    """
    parameter_file = "parameters.yml"
    parameters = {}

    with open(parameter_file, 'r') as param:

        for line in param:
            line = line.strip().split(": ")

            parameter_name = line[0]
            parameter = line[1]

            if parameter_name not in parameters:
                parameters[parameter_name] = parameter

    return parameters


def checking_parameters_of_the_scripts(call_command):
    """
    Checking that the parameter files for the given script exists or not
    """
    for param in call_command:

        if "/" in param and not param.endswith(".R") and not param.endswith(".py") and "." in param:
            if not os.path.isfile(param):
                sys.stdout.write(f" WARNING: One of the parameters of the script does not exist: "
                                 f"{param}\n\n")
                sys.exit(3)


def checking_errors():
    """
    Checking that the scripts have errors or not
    """
    with open("gutcovid.out", 'r') as log:

        for line in log:
            line = line.strip()

            if "Execution halted" in line:
                sys.stdout.write(f" WARNING: There was an error during the running! Please check the 'gutcovid.out' "
                                 f"file for more details.\n\n")
                sys.exit(4)


def run_script(script_parameters, parameters_for_the_script, script, step, call_script, command_log_file, script_number):
    """
    Function to run the given script
    """
    diffpath_scripts_carnival = ["Classical_pipeline_viral_proteins.R",
                        "Classical_pipeline_viral_miRNAs.R",
                        "CARNIVAL.R"]

    diffpath_scripts_virallink = ["filter_network_expressed_genes.R",
                                    "get_regulator_deg_network.R",
                                    "prepare_tiedie_input.R",
                                    "combined_edge_node_tables.R"]

    diffpath_scripts_tiedie = ["tiedie.py"]

    if script in diffpath_scripts_carnival:
        script = f"carnival/{script}"

    if script in diffpath_scripts_virallink:
        script = f"virallink/{script}" 

    if script in diffpath_scripts_tiedie:
        script = f"virallink/TieDie/{script}"
        if path.exists('/home/gutcovid/pipeline/tiedie_kernel.pkl'):
            os.remove('/home/gutcovid/pipeline/tiedie_kernel.pkl')
            print("The kernel file removed!")

    call_command = [call_script, f"{step}/{script}"]

    for p in parameters_for_the_script:

        if p in script_parameters:
            call_command.append(script_parameters[p])

        if p not in script_parameters:
            call_command.append(p)

        if "/" in p and "." in p:

            if not os.path.isfile(p):
                f = open(p, "a")
                f.close()

    checking_parameters_of_the_scripts(call_command)
    print(f'*** Command: {" ".join(call_command)} ***')
    with open(command_log_file, 'a') as command_log:
        command_log.write(f"#{script_number} " + " ".join(call_command) + '\n\n')
    time.sleep(1)
    run = subprocess.Popen(call_command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    run.communicate()

    checking_errors()

    print(f'*** [{strftime("%H:%M:%S")}] {script} finished successfully... ***\n')


def main():
    """
    Main function of the script
    """
    print(f'\n*** Gut-COVID pipeline ***')

    if not os.path.isfile("gutcovid.out"):
        open("gutcovid.out", 'a')

    if os.path.isfile("gutcovid.out"):
        os.remove("gutcovid.out")
        open("gutcovid.out", 'a')

    if not os.path.isfile("pipeline_command.log"):
        open("pipeline_command.log", 'a')

    if os.path.isfile("pipeline_command.log"):
        os.remove("pipeline_command.log")
        open("pipeline_command.log", 'a')

    script_parameters = get_parameters()
    print(f'\n*** Checking input parameters... ***')
    checking_input_parameters(script_parameters)
    print(f'*** Input parameters are fine, starting... ***\n')

    script_number = 1
    step_number = 1
    for step in scripts_folders:

        step_name_array = step.split("_")[1:]
        step_name = " ".join(step_name_array).upper()

        print(f"\n*** Step {step_number}/5 - {step_name} ***\n")
        for script in scripts_parameters:

            if step_number == 1: # Change the number to 2, when the 2nd folder is active

                if script == "tiedie.py":

                    if f"virallink/TieDie/{script}" not in scripts_folders[step]:
                        continue

                else:

                    if f"carnival/{script}" not in scripts_folders[step]:

                        if f"virallink/{script}" not in scripts_folders[step]:
                            continue

            if step_number != 1: # Change the number to 2, when the 2nd folder is active

                if script not in scripts_folders[step]:
                    continue

            for param in scripts_parameters[script]:

                print(f'*** [{strftime("%H:%M:%S")}] {script} ***')
                parameters_for_the_script = param

                if script.endswith(".R"):
                    call_script = "Rscript"
                    run_script(script_parameters, parameters_for_the_script, script, step, call_script, "pipeline_command.log", script_number)

                elif script.endswith(".py"):
                    call_script = "python3"
                    run_script(script_parameters, parameters_for_the_script, script, step, call_script, "pipeline_command.log", script_number)

                script_number += 1

        step_number += 1

    print(f'\n*** Gut-COVID pipeline was successfully finished! ***\n')


if __name__ == '__main__':
    main()
