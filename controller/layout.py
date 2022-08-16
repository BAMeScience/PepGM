""" Combine frames."""
from PySimpleGUI import Frame, Input, Button, Output, Text, theme

import frames

theme("SystemDefaultForReal")


def setup(ExperimentName="", SampleName="", HostName="", ScientificHostName="", ReferenceDBName="",
          SamplePath="", ParametersFile="", DataDir="", DatabaseDir="", PeptideShakerDir="", SearchGUIDir="",
          Alpha="[0.01, 0.05, 0.1, 0.2, 0.4, 0.6]", Beta="[0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 0.7]",
          prior="[0.1, 0.3, 0.5]", psmFDR=1, peptideFDR=1, proteinFDR=1, TaxaInPlot=15, TaxaInProteinCount=15,
          sourceDB="all[FILT]"):
    config_frame, run_frame, input_file_frame, input_dir_frame, gridsearch_frame, plotting_frame, pepgm_frame, \
    fdr_frame, searchgui_frame = frames.build("config.yaml", ExperimentName, SampleName, HostName,
                                              ScientificHostName, ReferenceDBName, SamplePath, ParametersFile,
                                              DataDir, DatabaseDir, PeptideShakerDir, SearchGUIDir, Alpha, Beta,
                                              prior, psmFDR, peptideFDR, proteinFDR, TaxaInPlot,
                                              TaxaInProteinCount, sourceDB)

    scaffold = [[Frame("Config file", [[config_frame]])],
                [Frame("Run", [[run_frame]])],
                [Frame("Input files", [[input_file_frame]])],
                [Frame("Input directories", [[input_dir_frame]])],
                [Frame("SearchGUI parameter", [[searchgui_frame]])],
                [Frame("PepGM parameter", [[pepgm_frame]])],
                [[Output(size=(73, 18))]],
                [Text("Cores"), Input(default_text=30, key="core_number", size=(10, 30)), Button('Dry run'),
                 Button('Run', button_color="LightGreen"), Button('Exit', button_color="red"),
                 Button("Help", button_color="orange")]]

    return scaffold
