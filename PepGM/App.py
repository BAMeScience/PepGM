#!/usr/bin/env python

from PySimpleGUI import Col, FolderBrowse, FileBrowse, Slider, InputText, Text, Combo, Radio, Frame, Input, Button
import os
import subprocess
import sys
import webbrowser
from pathlib import Path

import PySimpleGUI as sg
import yaml
from yaml.loader import SafeLoader

config_keys = ['ExperimentName', 'SampleName', 'HostName', 'ReferenceDBName', 'ScientificHostName', 'FilterSpectra',
               'AddHostandCrapToDB', 'SamplePath', 'ParametersFile', 'DataDir', 'DatabaseDir', 'PeptideShaker',
               'SearchGUI', 'ResourcesDir', 'ResultsDir', 'TaxidMapping', 'searchengines', 'psmFDR', 'peptideFDR',
               'proteinFDR', 'Alpha', 'Beta', 'prior', 'TaxaInPlot', 'TaxaInProteinCount', 'sourceDB', 'APImail',
               'APIkey', 'config_file_name']


def build_frames(experiment, sample, host, host_scientific, ref_database, path_sample, path_parameter,
                 path_data, path_database, peptide_shaker, search_gui, alpha, beta, prior, psm_fdr, peptide_fdr,
                 protein_fdr, n_taxa, mail, key):
    # CONFIG FRAME
    config_frame = Col([[Text("API mail", size=(12, 1)), InputText(key="APImail", default_text=mail, text_color="grey"),
                         Text('API key', size=(12, 1)), InputText(key="APIkey", default_text=key, text_color="grey")]])
    # RUN FRAME
    run_frame = Col([
        [InputText(key="ExperimentName", default_text=experiment, expand_x=True, enable_events=True)],
        [Text("Sample", size=(12, 1)),
         Input(key="SampleName", default_text=sample, expand_x=True),
         Text("Host", size=(12, 1)),
         Input(key="HostName", default_text=host, expand_x=True)],
        [Text("Reference", size=(12, 1)),
         Input(key="ReferenceDBName", default_text=ref_database, expand_x=True),
         Text("Scientific host", size=(12, 1)),
         Input(key="ScientificHostName", default_text=host_scientific, expand_x=True)],
        [Radio('Filter spectra', default=False, key="FilterSpectra", group_id=1),
         Radio('Add host and crap databases', default=True, key="AddHostandCrapToDB", group_id=1)]])
    # INPUT FRAMES
    input_file_frame = Col([
        [Text("Sample spectra", size=(20, 1)), InputText(key="SamplePath", default_text=path_sample, size=(90, 1)),
         FileBrowse()],
        [Text("Parameter", size=(20, 1), tooltip="Specify absolute path to parameter file (.par)"),
         Input(key="ParametersFile", default_text=path_parameter, size=(90, 1)), FileBrowse()]])
    input_dir_frame = Col(
        [[Text("Sample data", size=(20, 1)), Input(key="DataDir", default_text=path_data, size=(90, 1)),
          FolderBrowse()],
         [Text("Database", size=(20, 1)), Input(key="DatabaseDir", default_text=path_database, size=(90, 1)),
          FolderBrowse()],
         [Text("Peptide Shaker", size=(20, 1)),
          Input(key="PeptideShaker", default_text=peptide_shaker, size=(90, 1)), FolderBrowse()],
         [Text("Search GUI", size=(20, 1)), Input(key="SearchGUI", default_text=search_gui, size=(90, 1)),
          FolderBrowse()],
         [Text("Resources", size=(20, 1)), Input(key="ResourcesDir", default_text="resources/", size=(90, 1))],
         [Text("Results", size=(20, 1)), Input(key="ResultsDir", default_text="results/", size=(90, 1))],
         [Text("TaxID mapping", size=(20, 1)), Input(key="TaxidMapping", default_text="taxidMapping/", size=(90, 1))]])
    # PEPGM FRAME
    gridsearch_frame = Col([[Text("Alpha", size=(12, 2)), InputText(default_text=alpha, key="Alpha")],
                            [Text("Beta", size=(12, 2)), Input(default_text=beta, key="Beta")],
                            [Text("Prior", size=(12, 2)), Input(default_text=prior, key="prior")]])
    pepgm_frame = Col([[Text("# Taxa", size=(12, 1)), Input(default_text=n_taxa, key="TaxaInPlot")],
                       [Frame("Grid search", [[gridsearch_frame]], expand_y=True)]], expand_x=True, expand_y=True)
    fdr_frame = Col([[Text("PSM", size=(10, 1)),
                      Slider(range=(0, 20), default_value=psm_fdr, resolution=1, orientation='horizontal', key="psmFDR",
                             size=(40, 10))],
                     [Text("Peptide", size=(10, 1)),
                      Slider(range=(0, 20), default_value=peptide_fdr, resolution=1, orientation='horizontal',
                             key="peptideFDR", size=(40, 10))],
                     [Text("Protein", size=(10, 1)),
                      Slider(range=(0, 20), default_value=protein_fdr, resolution=1, orientation='horizontal',
                             key="proteinFDR", size=(40, 10))]])
    searchgui_frame = Col([[Text("Search engine", size=(15, 1)), Combo(list(
        ["-xtandem", "-myrimatch", "-ms_amanda", "-msgf", "-omssa", "-comet", "-tide", "-andromeda", "-meta_morpheus",
         "-novor", "-directtag"]), readonly=True, key="searchengines", default_value="-xtandem")],
                           [Frame("FDR calculation", [[fdr_frame]], expand_x=True, expand_y=True)]], expand_x=True,
                          expand_y=True)

    return config_frame, run_frame, input_file_frame, input_dir_frame, pepgm_frame, \
           searchgui_frame


def setup_layout(experiment="", sample="", host="", host_scientific="", ref_database="", path_sample="",
                 path_parameter="", path_data="", path_database="", peptide_shaker="", search_gui="",
                 alpha=[0.01, 0.05, 0.1, 0.2, 0.4, 0.6], beta=[0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 0.7],
                 prior=[0.1, 0.3, 0.5], psm_fdr=1, peptide_fdr=1, protein_fdr=1, n_taxa=15, key="", mail=""):
    config_frame, run_frame, input_file_frame, input_dir_frame, pepgm_frame, searchgui_frame = build_frames(
        experiment, sample, host, host_scientific, ref_database, path_sample, path_parameter, path_data, path_database,
        peptide_shaker, search_gui, list(alpha), list(beta), list(prior), psm_fdr, peptide_fdr, protein_fdr, n_taxa,
        key, mail)
    scaffold = [[Frame("Configuration details", [[config_frame]], expand_x=True)],
                [Frame("Run details", [[run_frame]], tooltip="Run configuration", expand_x=True)],
                [Frame("Input paths", [[input_file_frame], [input_dir_frame]], expand_x=True)],
                [Frame("Search parameter", [[searchgui_frame]], expand_x=True, size=(200, 185)),
                 Frame("PepGM settings", [[pepgm_frame]], expand_x=True, size=(200, 185))],
                [Text("Cores"), Input(default_text=30, key="core_number"), Button("Write", button_color="grey"),
                 Button('Dry run'), Button('Run', button_color="LightGreen"), Button('Exit', button_color="red"),
                 Button("Help", button_color="orange")]]

    return scaffold


def parse_config(configurations, path_config):
    """
    Write configuration from GUI into config file.
    :param configurations: list, configuration as retrieved from GUI
    :param path_config: str, path to config file
    """
    with open(path_config, "w") as config_file_out:
        # iterate over all dictionary keys and entries
        for param, value in configurations.items():
            # integer
            if param in ["TaxaInPlot", "TaxaInProteinCount"]:
                config_file_out.write("%s: %i\n" % (param, int(value)))
            # bool
            elif param in ["AddHostandCrapToDB", "FilterSpectra"]:
                config_file_out.write("%s: %r\n" % (param, bool(value)))
            # lists needs different formatting
            elif param in ["Alpha", "Beta", "prior"]:
                # changed values
                if "(" not in configurations[param]:
                    values = [float(x) for x in configurations[param].split()]
                # unchanged values
                elif "(" in configurations[param]:
                    values = configurations[param].replace("(", "").replace(")", "").split(",")
                # gridsearch parameter have to be saved in a list with the following format
                # [0.01, 0.05, 0.10, 0.20, 0.40]
                # opening bracket
                config_file_out.write("%s: [" % param)
                for idx, val in enumerate(values):
                    if idx is not len(values) - 1:
                        # max 2 digits
                        config_file_out.write("%.2f," % float(val))
                    else:
                        config_file_out.write("%.2f" % float(val))
                # closing bracket
                config_file_out.write("]\n")
            elif param == "ScientificHostName":
                # PepGM format requirement: " 'homo sapiens' "
                config_file_out.write("%s: \"\'%s\'\" \n" % (param, value))
            else:
                config_file_out.write("%s: '%s'\n" % (param, value))


def main():
    sg.theme("SystemDefaultForReal")
    # wd = os.path.dirname(os.path.realpath(workflow.snakefile))
    # get config file path
    config_file = Path(os.path.normpath(os.path.join(os.path.dirname(__file__), '../config/config.yaml')))

    if not config_file.exists():
        scaffold = setup_layout()
    else:
        # load previous configurations
        prev_configs = yaml.load(config_file.read_text(), Loader=SafeLoader)
        # do not display single quotes
        prev_configs["ScientificHostName"] = prev_configs["ScientificHostName"].replace("'", "")
        # auto fill input from config file
        scaffold = setup_layout(prev_configs["ExperimentName"], prev_configs["SampleName"], prev_configs["HostName"],
                                prev_configs["ScientificHostName"], prev_configs["ReferenceDBName"],
                                prev_configs["SamplePath"],
                                prev_configs["ParametersFile"], prev_configs["DataDir"], prev_configs["DatabaseDir"],
                                prev_configs["PeptideShaker"], prev_configs["SearchGUI"], prev_configs["Alpha"],
                                prev_configs["Beta"],
                                prev_configs["prior"], prev_configs["psmFDR"], prev_configs["peptideFDR"],
                                prev_configs["proteinFDR"],
                                prev_configs["TaxaInPlot"], prev_configs["APImail"], prev_configs["APIkey"])
    # initialize window
    window = sg.Window(title="Run PepGM", layout=scaffold, resizable=True, scaling=True, grab_anywhere=True,
                       auto_size_text=True, auto_size_buttons=True, finalize=True)
    # event loop
    while True:
        event, values = window.Read()
        input_key_list = [key for key, value in window.key_dict.items()]
        # run button
        if event == 'Run':
            cores = int(values["core_number"])

            # print("Your snakemake command: ")
            # run_command(snakemake_cmd)
            snakemake_cmd = 'cd ..; snakemake --use-conda --conda-frontend conda --cores 30'
            # subprocess.Popen([f'{snakemake_cmd}'], shell=True)
            # subprocess.Popen(snakemake_cmd, shell=True)
            subprocess.Popen([f'{snakemake_cmd}'], shell=True)
        # dry run button
        if event == 'Dry run':
            snakemake_cmd = f"cd ..; snakemake -np"
            run_command(cmd=snakemake_cmd, window=window)
        # update config file
        if event == "Write":
            print("Write configs successfull.")
            configs = {key: val for key, val in values.items() if key in config_keys}
            parse_config(configs, config_file)
        # help pages are on GitHub
        if event == 'Help':
            print("You will be redirected to the GitHub help pages.")
            webbrowser.open("https://github.com/BAMeScience/PepGM/blob/master/readme.md")
        # exit
        if event == 'Exit':
            break
        # close window event (X)
        if event == sg.WINDOW_CLOSED:
            break


if __name__ == "__main__":
    main()
