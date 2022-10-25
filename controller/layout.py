""" Combine frames."""
from PySimpleGUI import Frame, Input, Button, Text, Multiline

import frames


def setup(ExperimentName="", SampleName="", HostName="", ScientificHostName="", ReferenceDBName="",
          SamplePath="", ParametersFile="", DataDir="", DatabaseDir="", PeptideShakerDir="", SearchGUIDir="",
          Alpha="[0.01, 0.05, 0.1, 0.2, 0.4, 0.6]", Beta=[0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 0.7],
          prior=[0.1, 0.3, 0.5], psmFDR=1, peptideFDR=1, proteinFDR=1, TaxaInPlot=15,
          APIkey="", APImail=""):
    config_frame, run_frame, input_file_frame, input_dir_frame, pepgm_frame, searchgui_frame = frames.build(
        "config.yaml", ExperimentName, SampleName, HostName,
        ScientificHostName, ReferenceDBName, SamplePath, ParametersFile,
        DataDir, DatabaseDir, PeptideShakerDir, SearchGUIDir, list(Alpha),
        list(Beta), list(prior), psmFDR, peptideFDR, proteinFDR, TaxaInPlot, APIkey, APImail)
    scaffold = [[Frame("Configuration details", [[config_frame]], expand_x=True)],
                [Frame("Run details", [[run_frame]], tooltip="Run configuration", expand_x=True)],
                [Frame("Input paths", [[input_file_frame], [input_dir_frame]], expand_x=True)],
                [Frame("Search parameter", [[searchgui_frame]], expand_x=True, size=(200, 185)),
                 Frame("PepGM settings", [[pepgm_frame]], expand_x=True, size=(200, 185))],
                [[Multiline(expand_x=True, expand_y=True, reroute_stdout=True, reroute_stderr=True,
                            echo_stdout_stderr=True, autoscroll=True, size=(50, 30))]],
                [Text("Cores"), Input(default_text=30, key="core_number"), Button("Write", button_color="grey"),
                 Button('Dry run'), Button('Run', button_color="LightGreen"), Button('Exit', button_color="red"),
                 Button("Help", button_color="orange")]]

    return scaffold
