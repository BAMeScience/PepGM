""" Combine frames."""
from PySimpleGUI import Frame, Input, Button, Output, Text, theme

import frames

theme("SystemDefaultForReal")

scaffold = [[Frame("Config file", [[frames.config_frame]])],
            [Frame("Run", [[frames.run_frame]])],
            [Frame("Input files", [[frames.input_file_frame]])],
            [Frame("Input directories", [[frames.input_dir_frame]])],
            [Frame("SearchGUI parameter", [[frames.searchgui_frame]])],
            [Frame("PepGM parameter", [[frames.pepgm_frame]])],
            [[Output(size=(73, 18))]],
            [Text("Cores"), Input(default_text=30, key="core_number", size=(10, 30)), Button('Dry run'),
             Button('Run', button_color="LightGreen"), Button('Exit', button_color="red"),
             Button("Help", button_color="orange")]]
