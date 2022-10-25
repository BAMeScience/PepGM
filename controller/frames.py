""" Build GUI layout. """
from PySimpleGUI import Col, Input, Frame, FolderBrowse, FileBrowse, Slider, InputText, Text, Combo, Radio


def build(ConfigName, ExperimentName, SampleName, HostName, ScientificHostName, ReferenceDBName, SamplePath,
          ParametersFile, DataDir, DatabaseDir, PeptideShakerDir, SearchGUIDir, Alpha, Beta, prior, psmFDR,
          peptideFDR, proteinFDR, TaxaInPlot, APImail, APIkey):
    ## CONFIGURATION FRAME
    # API credentials
    config_frame = Col(
        [[Text("API mail", size=(12, 1)), InputText(key="APImail", default_text=APImail, text_color="grey"),
          Text('API key', size=(12, 1)), InputText(key="APIkey", default_text=APIkey, text_color="grey")]])
    # RUN FRAME
    # configures folder structure
    run_frame = Col([
        [InputText(key="ExperimentName", default_text=ExperimentName, expand_x=True, enable_events=True)],
        [Text("Sample", size=(12, 1)),
         Input(key="SampleName", default_text=SampleName, expand_x=True),
         Text("Host", size=(12, 1)),
         Input(key="HostName", default_text=HostName, expand_x=True)],
        [Text("Reference", size=(12, 1)),
         Input(key="ReferenceDBName", default_text=ReferenceDBName, expand_x=True),
         Text("Scientific host", size=(12, 1)),
         Input(key="ScientificHostName", default_text=ScientificHostName, expand_x=True)],
        [Radio('Filter spectra', default=False, key="FilterSpectra", group_id=1),
         Radio('Add host and crap databases', default=True, key="AddHostandCrapToDB", group_id=1)]])
    # INPUT FRAMES
    # locates input files
    input_file_frame = Col([
        [Text("Sample spectra", size=(20, 1)), InputText(key="SamplePath", default_text=SamplePath, size=(90, 1)),
         FileBrowse()],
        [Text("Parameter", size=(20, 1), tooltip="Specify absolute path to parameter file (.par)"),
         Input(key="ParametersFile", default_text=ParametersFile, size=(90, 1)), FileBrowse()]])
    input_dir_frame = Col(
        [[Text("Sample data", size=(20, 1)), Input(key="DataDir", default_text=DataDir, size=(90, 1)), FolderBrowse()],
         [Text("Database", size=(20, 1)), Input(key="DatabaseDir", default_text=DatabaseDir, size=(90, 1)),
          FolderBrowse()],
         [Text("Peptide Shaker", size=(20, 1)),
          Input(key="PeptideShakerDir", default_text=PeptideShakerDir, size=(90, 1)), FolderBrowse()],
         [Text("Search GUI", size=(20, 1)), Input(key="SearchGUIDir", default_text=SearchGUIDir, size=(90, 1)),
          FolderBrowse()],
         [Text("Resources", size=(20, 1)), Input(key="ResourcesDir", default_text="resources/", size=(90, 1))],
         [Text("Results", size=(20, 1)), Input(key="ResultsDir", default_text="results/", size=(90, 1))],
         [Text("TaxID mapping", size=(20, 1)), Input(key="TaxidMapping", default_text="taxidMapping/", size=(90, 1))]])
    # PEPGM FRAME
    # configure your PepGM run
    gridsearch_frame = Col([[Text("Alpha", size=(12, 2)), InputText(default_text=Alpha, key="Alpha")],
                            [Text("Beta", size=(12, 2)), Input(default_text=Beta, key="Beta")],
                            [Text("Prior", size=(12, 2)), Input(default_text=prior, key="prior")]])
    pepgm_frame = Col([[Text("# Taxa", size=(12, 1)), Input(default_text=TaxaInPlot, key="TaxaInPlot")],
                       [Frame("Grid search", [[gridsearch_frame]], expand_y=True)]], expand_x=True, expand_y=True)
    # configure FDR settings
    fdr_frame = Col([[Text("PSM", size=(10, 1)),
                      Slider(range=(0, 20), default_value=psmFDR, resolution=1, orientation='horizontal', key="psmFDR",
                             size=(40, 10))],
                     [Text("Peptide", size=(10, 1)),
                      Slider(range=(0, 20), default_value=peptideFDR, resolution=1, orientation='horizontal',
                             key="peptideFDR", size=(40, 10))],
                     [Text("Protein", size=(10, 1)),
                      Slider(range=(0, 20), default_value=proteinFDR, resolution=1, orientation='horizontal',
                             key="proteinFDR", size=(40, 10))]])
    searchgui_frame = Col([[Text("Search engine", size=(15, 1)), Combo(list(
        ["-xtandem", "-myrimatch", "-ms_amanda", "-msgf", "-omssa", "-comet", "-tide", "-andromeda", "-meta_morpheus",
         "-novor", "-directtag"]), readonly=True, key="searchengines", default_value="-xtandem")],
                           [Frame("FDR calculation", [[fdr_frame]], expand_x=True, expand_y=True)]], expand_x=True,
                          expand_y=True)

    return config_frame, run_frame, input_file_frame, input_dir_frame, pepgm_frame, \
           searchgui_frame
