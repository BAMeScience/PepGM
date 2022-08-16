import json
import os
import re
import subprocess
import sys

import PySimpleGUI as sg
import yaml

from layout import scaffold
import webbrowser

sg.theme("SystemDefaultForReal")

def run_command(cmd, timeout=None, window=None):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = ""
    for line in p.stdout:
        line = line.decode(errors='replace' if sys.version_info < (3, 5) else 'backslashreplace').rstrip()
        output += line
        print(line)
        window.Refresh() if window else None
    retval = p.wait(timeout)
    return retval, output


def parse_config(values):
    config_file_name = values["config_file_name"]

    # remove old config file
    if os.path.exists(f"../config/{config_file_name}"):
        os.remove(f"../config/{config_file_name}")

    # remove unnecessary entries from gui
    configs = {k: values[k] for k in values if not re.match('^Browse', k)}

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
        elif key == "ScientificHostName":
            configs[key] = '"' + configs[key] + '"'
        elif key in ("PeptideShakerDir", "SearchGUIDir", "DataDir", "DatabaseDir"):
            configs[key] = str(configs[key] + '/')
        elif key in ("searchengines", "ResourcesDir", "ResultsDir", "TaxidMapping"):
            configs[key] = str(configs[key])

    with open(f"../config/{config_file_name}", "w") as outfile:
        yaml.safe_dump(configs, outfile, default_flow_style=None)


# initialize window
window = sg.Window(title="Run PepGM", layout=scaffold, resizable=True)

# main
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

    if event == 'Help':
        #TODO: replace private GitLab by public GitHub link
        webbrowser.open("https://git.bam.de/tholstei/pepgm/-/blob/master/readme.md")

    if event in (None, 'Exit'):
        break
