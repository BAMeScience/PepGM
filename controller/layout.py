""" Combine frames."""
from PySimpleGUI import Frame, Input, Button, Text, theme, Multiline

import frames

theme("SystemDefaultForReal")


def setup(ExperimentName="", SampleName="", HostName="", ScientificHostName="", ReferenceDBName="",
          SamplePath="", ParametersFile="", DataDir="", DatabaseDir="", PeptideShakerDir="", SearchGUIDir="",
          Alpha=[0.01, 0.05, 0.1, 0.2, 0.4, 0.6], Beta=[0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 0.7],
          prior=[0.1, 0.3, 0.5], psmFDR=1, peptideFDR=1, proteinFDR=1, TaxaInPlot=15, TaxaInProteinCount=15,
          sourceDB="all[FILT]", APIkey="", APImail=""):
    config_frame, run_frame, input_file_frame, input_dir_frame, gridsearch_frame, plotting_frame, pepgm_frame, \
    fdr_frame, searchgui_frame = frames.build("config.yaml", ExperimentName, SampleName, HostName,
                                              ScientificHostName, ReferenceDBName, SamplePath, ParametersFile,
                                              DataDir, DatabaseDir, PeptideShakerDir, SearchGUIDir, list(Alpha), list(Beta),
                                              list(prior), psmFDR, peptideFDR, proteinFDR, TaxaInPlot,
                                              TaxaInProteinCount, sourceDB, APIkey, APImail)
    scaffold = [[Frame("Run details", [[run_frame]], tooltip="Run configuration", expand_x=True)],
                [Frame("Input paths", [[input_file_frame], [input_dir_frame]], expand_x=True),
                 Frame("Search parameter", [[searchgui_frame]], expand_x=True)],
                [Frame("PepGM parameter", [[pepgm_frame]], expand_x=True)],
                [Frame("Configuration details", [[config_frame]], expand_x=True)],
                [[Multiline(expand_x=True, expand_y=True, reroute_stdout=True, reroute_stderr=True,
                            echo_stdout_stderr=True, autoscroll=True, size=12)]],
                [Text("Cores"), Input(default_text=30, key="core_number"), Button("Read", button_color="grey"),
                 Button('Dry run'), Button('Run', button_color="LightGreen"), Button('Stop', button_color="red"),
                 Button("Help", button_color="orange")]]

    return scaffold
