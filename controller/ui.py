""""Main for GUI."""
import subprocess
import sys
import webbrowser
from os.path import exists as file_exists

import PySimpleGUI as sg
import yaml
from yaml.loader import SafeLoader

import layout
import const


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


def parse_config(configs, config_file_path="../config/config.yaml"):
    """
    Parse configuration from GUI into config file.
    :param configs: list, configuration as retrieved from GUI
    :param config_file_path: str, path to config file
    """
    with open(config_file_path, 'w') as config_file_out:
        # iterate over all dictionary keys and entries that are retrieved from the GUI
        for key, value in configs.items():
            # integer
            if key in ["TaxaInPlot", "TaxaInProteinCount"]:
                config_file_out.write("%s: %i\n" % (key, int(value)))
            # bool
            elif key in ["AddHostandCrapToDB", "FilterSpectra"]:
                config_file_out.write("%s: %r\n" % (key, bool(value)))
            # lists needs different formatting
            elif key in ["Alpha", "Beta", "prior"]:
                # gridsearch parameter have to be saved in a list with the following format
                # [0.010000, 0.050000, 0.100000, 0.200000, 0.400000]
                config_file_out.write("%s: [" % key)
                # init
                temp = []
                digit = ""
                # merge digits into a float
                for element in configs[key]:
                    if element == '(' or element == ')':
                        continue
                    if element.isdigit() or element == ".":
                        digit = digit + element
                    if element == ",":
                        temp.append(float(digit))
                        digit = ""
                # write merged digits
                for element in temp:
                    config_file_out.write("%f," % element)
                # write end of array
                config_file_out.write("]\n")
            # use bool format
            elif key == "ScientificHostName":
                config_file_out.write('%s: '"%s"'\n' % (key, value))
            # use string format
            else:
                config_file_out.write("%s: '%s'\n" % (key, value))


if __name__ == "__main__":
    sg.theme("SystemDefaultForReal")
    # load settings if config file is present
    # file has to have name 'config.yaml'
    if file_exists("../config/config.yaml"):
        # load yaml
        with open("../config/config.yaml") as f:
            prev_configs = yaml.load(f, Loader=SafeLoader)
        # auto fill input
        scaffold = layout.setup(prev_configs["ExperimentName"], prev_configs["SampleName"], prev_configs["HostName"],
                                prev_configs["ScientificHostName"], prev_configs["ReferenceDBName"],
                                prev_configs["SamplePath"],
                                prev_configs["ParametersFile"], prev_configs["DataDir"], prev_configs["DatabaseDir"],
                                prev_configs["PeptideShakerDir"], prev_configs["SearchGUIDir"], prev_configs["Alpha"],
                                prev_configs["Beta"],
                                prev_configs["prior"], prev_configs["psmFDR"], prev_configs["peptideFDR"],
                                prev_configs["proteinFDR"],
                                prev_configs["TaxaInPlot"], prev_configs["TaxaInProteinCount"],
                                prev_configs["sourceDB"], prev_configs["APImail"],
                                prev_configs["APIkey"])
    # set up fresh GUI else
    else:
        scaffold = layout.setup()
    # initialize window
    window = sg.Window(title="Run PepGM", layout=scaffold, resizable=True, scaling=True, grab_anywhere=True,
                       auto_size_text=True, auto_size_buttons=True)
    # event loop
    while True:
        event, values = window.Read()
        input_key_list = [key for key, value in window.key_dict.items()]
        # run button
        if event == 'Run':
            cores = int(values["core_number"])
            snakemake_cmd = f"cd ..; snakemake --cores {cores} -p"
            print(snakemake_cmd)
            run_command(cmd=snakemake_cmd, window=window)
        # dry run button
        if event == 'Dry run':
            snakemake_cmd = f"cd ..; snakemake -np"
            run_command(cmd=snakemake_cmd, window=window)
        # update config file
        if event == "Read":
            configs = {key: val for key, val in values.items() if key in const.keys_to_keep}
            parse_config(configs)
        # help pages are on GitHub
        if event == 'Help':
            webbrowser.open("https://github.com/BAMeScience/PepGM/blob/master/readme.md")
        # exit
        if event == 'Exit':
            break
