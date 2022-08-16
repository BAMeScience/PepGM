""""Main for GUI."""

import json
import os
import re
import subprocess
import sys
import webbrowser
from os.path import exists as file_exists

import PySimpleGUI as sg
import yaml

import layout


def run_command(cmd, timeout=None, window=None):
    """
    Display snakemake run details in a windows inside the GUI.
    """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = ""
    for line in p.stdout:
        line = line.decode(errors='replace' if sys.version_info < (3, 5) else 'backslashreplace').rstrip()
        output += line
        print(line)
        window.Refresh() if window else None
    retval = p.wait(timeout)
    return retval, output


def parse_config(settings):
    """
    Parse user input into a list that creates config file.
    """
    config_file_name = settings["config_file_name"]

    # remove old config file
    if os.path.exists(f"../config/{config_file_name}"):
        os.remove(f"../config/{config_file_name}")

    # remove unnecessary entries
    configs = {k: settings[k] for k in settings if not re.match('^Browse', k)}

    # parser
    for key, _ in configs.items():
        if key in ('Alpha', 'Beta', 'prior'):
            configs[key] = json.loads(configs[key])
        elif key in ("AddHostandCrapToDB", "FilterSpectra"):
            configs[key] = bool(configs[key])
        elif key in ("peptideFDR", "proteinFDR", "psmFDR"):
            configs[key] = str(int(configs[key]))
        elif key in ("TaxaInPlot", "TaxaInProteinCount"):
            configs[key] = int(configs[key])
        elif key in ("searchengines", "ResourcesDir", "ResultsDir", "TaxidMapping"):
            configs[key] = str(configs[key])

    # save
    with open(f"../config/{config_file_name}", "w") as outfile:
        yaml.safe_dump(configs, outfile, default_flow_style=None)


if __name__ == "__main__":
    sg.theme("SystemDefaultForReal")

    # load config settings if already present
    if file_exists("../config/config.yaml"):
        # load yaml
        with open("../config/config.yaml", 'r') as config_file:
            configs = yaml.load(config_file, Loader=yaml.FullLoader)
        scaffold = layout.setup(configs["ExperimentName"], configs["SampleName"], configs["HostName"],
                                configs["ScientificHostName"], configs["ReferenceDBName"], configs["SamplePath"],
                                configs["ParametersFile"], configs["DataDir"], configs["DatabaseDir"],
                                configs["PeptideShakerDir"], configs["SearchGUIDir"], configs["Alpha"], configs["Beta"],
                                configs["prior"], configs["psmFDR"], configs["peptideFDR"], configs["proteinFDR"],
                                configs["TaxaInPlot"], configs["TaxaInProteinCount"], configs["sourceDB"])
    # set up fresh GUI else
    else:
        scaffold = layout.setup()

# initialize window
window = sg.Window(title="Run PepGM", layout=scaffold, resizable=True)
while True:
    event, values = window.Read()

    # run button
    if event == 'Run':
        parse_config(values)
        cores = int(values["core_number"])
        snakemake_cmd = f"cd ..; snakemake --cores {cores} -p"
        print(snakemake_cmd)
        run_command(cmd=snakemake_cmd, window=window)

    # dry run button
    if event == 'Dry run':
        parse_config(values)
        snakemake_cmd = f"cd ..; snakemake -np"
        run_command(cmd=snakemake_cmd, window=window)

    # help pages are on GitHub
    if event == 'Help':
        webbrowser.open("https://github.com/BAMeScience/PepGM/blob/master/readme.md")

    # exit
    if event in (None, 'Exit'):
        break
