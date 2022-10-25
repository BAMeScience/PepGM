""""Main for GUI."""
import os
import subprocess
import sys
import webbrowser
from pathlib import Path

import PySimpleGUI as sg
import yaml
from yaml.loader import SafeLoader

import const
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


def parse_config(configurations, config_file_path, prev_configs):
    """
    Write configuration from GUI into config file.
    :param configurations: list, configuration as retrieved from GUI
    :param config_file_path: str, path to config file
    """
    with open(config_file_path, "w") as config_file_out:
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


if __name__ == "__main__":
    sg.theme("SystemDefaultForReal")
    # get config file path
    config_file = Path(os.path.normpath(os.path.join(os.path.dirname(__file__), '../config/config.yaml')))

    if not config_file.exists():
        scaffold = layout.setup()
    else:
        # load previous configurations
        prev_configs = yaml.load(config_file.read_text(), Loader=SafeLoader)
        # do not display single quotes
        prev_configs["ScientificHostName"] = prev_configs["ScientificHostName"].replace("'", "")
        # auto fill input from config file
        scaffold = layout.setup(prev_configs["ExperimentName"], prev_configs["SampleName"], prev_configs["HostName"],
                                prev_configs["ScientificHostName"], prev_configs["ReferenceDBName"],
                                prev_configs["SamplePath"],
                                prev_configs["ParametersFile"], prev_configs["DataDir"], prev_configs["DatabaseDir"],
                                prev_configs["PeptideShakerDir"], prev_configs["SearchGUIDir"], prev_configs["Alpha"],
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
            snakemake_cmd = f"cd ..; snakemake --cores {cores} -p"
            print("Your snakemake command: ")
            print(snakemake_cmd)
            run_command(cmd=snakemake_cmd, window=window)
        # dry run button
        if event == 'Dry run':
            snakemake_cmd = f"cd ..; snakemake -np"
            run_command(cmd=snakemake_cmd, window=window)
        # update config file
        if event == "Write":
            configs = {key: val for key, val in values.items() if key in const.keys_to_keep}
            parse_config(configs, config_file, prev_configs)
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
