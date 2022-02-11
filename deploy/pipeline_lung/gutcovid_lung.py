# Main python script to run the whole Gut-COVID pipeline on lung data
from distutils.dir_util import copy_tree
from time import strftime
from os import path
import subprocess
import sys
import os
import time


scripts_folders = {
    "1_data_preparation": ["1_creating_lung_data_tables.R",
                            "2_extract_matrices.R",
                            "3_initial_data_mining.R"],
    "2_Inter_analysis": ["2_filter_expression_gaussian.py",
                            "3_log2_zscore_filter.py",
                            "4_epithelial_immune_cells_interaction_lung.py",
                            "5_ligand_imm_cell_network_creation.R",
                            "6_ligand_receptor_stats.R",
                            "7_barplot_num_ints_per_epi_ligand.R",
                            "8_barplot_num_ints_per_receptor_num_imm_cells.R",
                            "11_cell_cell_heatplots_cellspecific.R",
                            "12_ligand_imm_cell_heatplots_receptor_expression.R",
                            "14_ligand_receptor_heatplots.R",
                            "15_functional_analysis_ligand_receptors.R",
                            "16_ligand_receptor_family_network_creation.R",
                            "17_filter_degs_for_ligands.R"],
    "3_Intra_analysis_ViralLink": ["filter_network_expressed_genes.R",
                                    "get_regulator_deg_network.R",
                                    "prepare_tiedie_input.R",
                                    "TieDie/tiedie.py",
                                    "combined_edge_node_tables_for_miRNA.R",
                                    "combined_edge_node_tables.R"],
    "5_Intra_functional_analysis": ["Functional_analysis_from_Virallink_mirna.R",
                                    "Functional_analysis_from_Virallink_protein.R",
                                    "Functional_analysis_from_Virallink_shared.R"]
}

scripts_parameters = {
    "1_creating_lung_data_tables.R": [["R_object_input_data",
                                        "/home/gutcovid/pipeline_lung/1_data_preparation/metadata_covid_nb_main.txt",
                                        "/home/gutcovid/pipeline_lung/1_data_preparation/deg_lung_squamousous_alveolar_cells_control_vs_critcal.txt",
                                        "/home/gutcovid/pipeline_lung/1_data_preparation/deg_lung_secreatory_alveolar_cells_control_vs_critcal.txt",
                                        "/home/gutcovid/pipeline_lung/1_data_preparation/deg_lung_ciliated_alveolar_cells_control_vs_critcal.txt",
                                        "/home/gutcovid/pipeline_lung/1_data_preparation/SEURAT_object_normed.rds"]],
    "2_extract_matrices.R": [["/home/gutcovid/pipeline_lung/1_data_preparation/SEURAT_object_normed.rds",
                                "/home/gutcovid/pipeline_lung/1_data_preparation/metadata_covid_nb_main.txt"]],
    "3_initial_data_mining.R": [["/home/gutcovid/pipeline_lung/1_data_preparation/COVID_lung_dataset_nb_mainv2_average_epxression.tsv",
                                "/home/gutcovid/pipeline_lung/1_data_preparation/average_immune_expression.csv",
                                "/home/gutcovid/pipeline_lung/1_data_preparation/average_epithelial_ciliated_expression.csv"]],
    "2_filter_expression_gaussian": [["-i", "average_epithelial_ciliated_expression.txt",
                                        "-o", "/home/gutcovid/pipeline_lung/2_Inter_analysis/2_filter_expression_gaussian"]],
    "3_log2_zscore_filter.py": [["-i", "/home/gutcovid/pipeline_lung/1_data_preparation/average_immune_expression.csv",
                                "-o", "/home/gutcovid/pipeline_lung/2_Inter_analysis/3_log2_zscore_filter/immune_healthy_expressed_genes.csv"]],
    "4_epithelial_immune_cells_interaction_lung.py": [["-imm", "/home/gutcovid/pipeline_lung/2_Inter_analysis/3_log2_zscore_filter/immune_healthy_expressed_genes.csv",
                                                        "-epi", "/home/gutcovid/pipeline_lung/1_data_preparation/degs_filt_lung_ciliated_cells_adjpval0.05_lfc0.5_inter.txt",
                                                        "--omni" "omnipath_file",
                                                        "-o", "/home/gutcovid/pipeline_lung/2_Inter_analysis/4_epithelial_immune_cells_interaction_lung/ligand_receptor_ints_lung_epithelial_immune_cells.txt"]],
    "5_ligand_imm_cell_network_creation.R": [["/home/gutcovid/pipeline_lung/2_Inter_analysis/4_epithelial_immune_cells_interaction_lung/ligand_receptor_ints_lung_epithelial_immune_cells.txt",
                                                "/home/gutcovid/pipeline_lung/2_Inter_analysis/5_ligand_imm_cell_network_creation/lung_ciliated_ligand_imm_cell_network.txt"]],
    "6_ligand_receptor_stats.R": [["/home/gutcovid/pipeline_lung/2_Inter_analysis/4_epithelial_immune_cells_interaction_lung/ligand_receptor_ints_lung_epithelial_immune_cells.txt",
                                    "/home/gutcovid/pipeline_lung/2_Inter_analysis/6_ligand_receptor_stats/stats_tog_ligand_receptor_ints_lung_ciliated_immune_cells.txt",
                                    "/home/gutcovid/pipeline_lung/2_Inter_analysis/6_ligand_receptor_stats/stats_sep_ligand_receptor_ints_lung_ciliated_immune_cells.txt",
                                    "/home/gutcovid/pipeline_lung/2_Inter_analysis/6_ligand_receptor_stats/stats2_ligand_receptor_ints_lung_ciliated_immune_cells.txt"]],
    "7_barplot_num_ints_per_epi_ligand.R": [["/home/gutcovid/pipeline_lung/2_Inter_analysis/4_epithelial_immune_cells_interaction_lung/ligand_receptor_ints_lung_epithelial_immune_cells.txt",
                                            "/home/gutcovid/pipeline_lung/2_Inter_analysis/7_barplot_num_ints_per_epi_ligand"]],
    "8_barplot_num_ints_per_receptor_num_imm_cells.R": [["/home/gutcovid/pipeline_lung/2_Inter_analysis/4_epithelial_immune_cells_interaction_lung/ligand_receptor_ints_lung_epithelial_immune_cells.txt",
                                            "/home/gutcovid/pipeline_lung/2_Inter_analysis/7_barplot_num_ints_per_epi_ligand"]],
    "11_cell_cell_heatplots_cellspecific.R": [["/home/gutcovid/pipeline_lung/2_Inter_analysis/6_ligand_receptor_stats/stats_sep_ligand_receptor_ints_lung_ciliated_immune_cells.txt",
                                                "/home/gutcovid/pipeline_lung/2_Inter_analysis/11_cell_cell_heatplots_cellspecific/ciliated_immune_cell_cell_heatplot_log2_num_interactions.png"]],
    "12_ligand_imm_cell_heatplots_receptor_expression.R": [["/home/gutcovid/pipeline_lung/2_Inter_analysis/5_ligand_imm_cell_network_creation/lung_ciliated_ligand_imm_cell_network.txt",
                                                            "/home/gutcovid/pipeline_lung/2_Inter_analysis/12_ligand_imm_cell_heatplots_receptor_expression"]],
    "14_ligand_receptor_heatplots.R": [["/home/gutcovid/pipeline_lung/2_Inter_analysis/6_ligand_receptor_stats/stats2_ligand_receptor_ints_lung_ciliated_immune_cells.txt",
                                        "/home/gutcovid/pipeline_lung/2_Inter_analysis/14_ligand_receptor_heatplots"]],
    "15_functional_analysis_ligand_receptors.R": [["/home/gutcovid/pipeline_lung/2_Inter_analysis/4_epithelial_immune_cells_interaction_lung/ligand_receptor_ints_lung_epithelial_immune_cells.txt",
                                                    "omnipath_file",
                                                    "/home/gutcovid/pipeline_lung/2_Inter_analysis/15_functional_analysis_ligand_receptors"]],
    "16_ligand_receptor_family_network_creation.R": [["/home/gutcovid/pipeline_lung/2_Inter_analysis/6_ligand_receptor_stats/stats2_ligand_receptor_ints_lung_ciliated_immune_cells.txt",
                                                        "receptor_categories",
                                                        "/home/gutcovid/pipeline_lung/2_Inter_analysis/16_ligand_receptor_family_network_creation/ligand_receptor_category_network.txt"]],
    "17_filter_degs_for_ligands.R": [["all_ligands",
                                        "/home/gutcovid/pipeline_lung/1_data_preparation/degs_filt_lung_ciliated_cells_adjpval0.05_lfc0.5_inter.txt",
                                        "/home/gutcovid/pipeline_lung/2_Inter_analysis/17_filter_degs_for_ligands/ligands_lung_ciliary_cells_adjpval0.05_lfc0.5.txt",
                                        "/home/gutcovid/pipeline_lung/2_Inter_analysis/17_filter_degs_for_ligands/ligands_lung_ciliary_cells_adjpval0.05_lfc0.5.csv"]],
    "filter_network_expressed_genes.R": [["/home/gutcovid/pipeline_lung/2_Inter_analysis/2_filter_expression_gaussian/epithelial_ciliated_moderate_expressed_genes.txt",
                                            "dorothea_file",
                                            "omnipath_file",
                                            "id_type",
                                            "outdir_miRNA"],
                                            ["/home/gutcovid/pipeline_lung/2_Inter_analysis/2_filter_expression_gaussian/epithelial_ciliated_moderate_expressed_genes.txt",
                                            "dorothea_file",
                                            "omnipath_file",
                                            "id_type",
                                            "outdir_proteins"]],
    "get_regulator_deg_network.R": [["/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_rna/dorothea_contextualised_network.txt",
                                    "/home/gutcovid/pipeline_lung/2_Inter_analysis/17_filter_degs_for_ligands/ligands_lung_ciliary_cells_adjpval0.05_lfc0.5.csv",
                                    "id_type",
                                    "outdir_miRNA"],
                                    ["/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_proteins/dorothea_contextualised_network.txt",
                                    "/home/gutcovid/pipeline_lung/2_Inter_analysis/17_filter_degs_for_ligands/ligands_lung_ciliary_cells_adjpval0.05_lfc0.5.csv",
                                    "id_type",
                                    "outdir_proteins"]],
    "prepare_tiedie_input.R": [["/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_rna/omnipath_contextualised_network.txt",
                                "/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_rna/contextualised_regulator-deg_network.txt",
                                "/home/gutcovid/pipeline_lung/2_Inter_analysis/17_filter_degs_for_ligands/ligands_lung_ciliary_cells_adjpval0.05_lfc0.5.csv",
                                "hbps_mirna",
                                "outdir_miRNA"],
                                ["/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_proteins/omnipath_contextualised_network.txt",
                                "/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_proteins/contextualised_regulator-deg_network.txt",
                                "/home/gutcovid/pipeline_lung/2_Inter_analysis/17_filter_degs_for_ligands/ligands_lung_ciliary_cells_adjpval0.05_lfc0.5.csv",
                                "hbps_protein",
                                "outdir_proteins"]],
    "tiedie.py": [["-u", "/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_rna/input_files/upstream.input",
                    "-d", "/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_rna/input_files/downstream.input",
                    "-n", "/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_rna/input_files/pathway.sif",
                    "-o", "/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_rna/TieDIE"],
                    ["-u", "/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_proteins/input_files/upstream.input",
                    "-d", "/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_proteins/input_files/downstream.input",
                    "-n", "/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_proteins/input_files/pathway.sif",
                    "-o", "/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_proteins/TieDIE"]],
    "combined_edge_node_tables_for_miRNA.R": [["/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_rna/TieDIE/tiedie.cn.sif",
                                    "/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_rna/TieDIE/heats.NA",
                                    "hbps_mirna",
                                    "sars",
                                    "/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_rna/contextualised_regulator-deg_network.txt",
                                    "/home/gutcovid/pipeline_lung/1_data_preparation/degs_unfilt_lung_ciliated_cells.csv",
                                    "id_type",
                                    "outdir_miRNA"]],
    "combined_edge_node_tables.R": [["/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_proteins/TieDIE/tiedie.cn.sif",
                                    "/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_proteins/TieDIE/heats.NA",
                                    "hbps_protein",
                                    "sars",
                                    "/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_proteins/contextualised_regulator-deg_network.txt",
                                    "/home/gutcovid/pipeline_lung/1_data_preparation/degs_unfilt_lung_ciliated_cells.csv",
                                    "id_type",
                                    "outdir_proteins"]],
    "Functional_analysis_from_Virallink_mirna.R": [["/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_proteins/omnipath_contextualised_network.txt",
                                                    "/home/gutcovid/pipeline_lung/4_Intra_create_merged_network/ciliated_cells_moderate_intra_merged_network_edited.csv"]],
    "Functional_analysis_from_Virallink_protein.R": [["/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_proteins/omnipath_contextualised_network.txt",
                                                    "/home/gutcovid/pipeline_lung/4_Intra_create_merged_network/ciliated_cells_moderate_intra_merged_network_edited.csv"]],
    "Functional_analysis_from_Virallink_shared.R": [["/home/gutcovid/pipeline_lung/3_Intra_analysis_ViralLink/output_directory_viral_proteins/omnipath_contextualised_network.txt",
                                                    "/home/gutcovid/pipeline_lung/4_Intra_create_merged_network/ciliated_cells_moderate_intra_merged_network_edited.csv"]]
}

def checking_input_parameters(script_parameters):
    """
    Checking the parameters
    """
    if not os.path.isfile("parameters_lung.yml"):
        sys.stdout.write(f" WARNING: There is no appropriate parameter file! It should be parameters_lung.yml\n")
        sys.exit(1)

    neccessary_parameters = ["R_object_input_data"]
    for nec_param in neccessary_parameters:
        if nec_param not in script_parameters:
            sys.stdout.write(f" WARNING: A parameter is missing from the parameters_lung.yml file: {nec_param}\n\n")
            sys.exit(2)


def get_parameters():
    """
    Get all of the parameters from the parameter_lung.yml file
    """
    parameter_file = "parameters_lung.yml"
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
    with open("gutcovid_lung.out", 'r') as log:

        for line in log:
            line = line.strip()

            if "Execution halted" in line:
                sys.stdout.write(f" WARNING: There was an error during the running! Please check the 'gutcovid_lung.out' "
                                 f"file for more details.\n\n")
                sys.exit(4)


def run_script(script_parameters, parameters_for_the_script, script, step, call_script, command_log_file, script_number):
    """
    Function to run the given script
    """
    diffpath_scripts_virallink = ["filter_network_expressed_genes.R",
                                    "get_regulator_deg_network.R",
                                    "prepare_tiedie_input.R",
                                    "combined_edge_node_tables_for_miRNA.R",
                                    "combined_edge_node_tables.R"]

    diffpath_scripts_tiedie = ["tiedie.py"]

    if script in diffpath_scripts_virallink:
        script = f"virallink/{script}" 

    if script in diffpath_scripts_tiedie:
        script = f"virallink/TieDie/{script}"
        if path.exists('/home/gutcovid/pipeline_lung/tiedie_kernel.pkl'):
            os.remove('/home/gutcovid/pipeline_lung/tiedie_kernel.pkl')
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
    print(f'\n*** Gut-COVID (lung) pipeline ***')

    if not os.path.isfile("gutcovid_lung.out"):
        open("gutcovid_lung.out", 'a')

    if os.path.isfile("gutcovid_lung.out"):
        os.remove("gutcovid_lung.out")
        open("gutcovid_lung.out", 'a')

    if not os.path.isfile("pipeline_command_lung.log"):
        open("pipeline_command_lung.log", 'a')

    if os.path.isfile("pipeline_command_lung.log"):
        os.remove("pipeline_command_lung.log")
        open("pipeline_command_lung.log", 'a')

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

            if step_number == 3:

                if script == "tiedie.py":

                    if f"virallink/TieDie/{script}" not in scripts_folders[step]:
                        continue

                else:

                    if f"carnival/{script}" not in scripts_folders[step]:

                        if f"virallink/{script}" not in scripts_folders[step]:
                            continue

            if step_number != 3:

                if script not in scripts_folders[step]:
                    continue

            for param in scripts_parameters[script]:

                print(f'*** [{strftime("%H:%M:%S")}] {script} ***')
                parameters_for_the_script = param

                if script.endswith(".R"):
                    call_script = "Rscript"
                    run_script(script_parameters, parameters_for_the_script, script, step, call_script, "pipeline_command_lung.log", script_number)

                elif script.endswith(".py"):
                    call_script = "python3"
                    run_script(script_parameters, parameters_for_the_script, script, step, call_script, "pipeline_command_lung.log", script_number)

                script_number += 1

        step_number += 1

    print(f'\n*** Gut-COVID (lung) pipeline was successfully finished! ***\n')


if __name__ == '__main__':
    main()
