""" Build GUI layout. """
from PySimpleGUI import Col, Input, Frame, FolderBrowse, FileBrowse, Slider, InputText, Text, Combo, Radio, theme

def build(ConfigName, ExperimentName, SampleName, HostName, ScientificHostName, ReferenceDBName, SamplePath,
          ParametersFile, DataDir, DatabaseDir, PeptideShakerDir, SearchGUIDir, Alpha, Beta, prior, psmFDR,
          peptideFDR, proteinFDR, TaxaInPlot, TaxaInProteinCount, sourceDB, APImail, APIkey):
    # configuration frame
    config_frame = Col([[Text("API mail"),
                         InputText(key="APImail", default_text=APImail, text_color="grey"), Text('API key'),
                         InputText(key="APIkey", default_text=APIkey, text_color="grey")]])
    # run frame
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
    # input frame
    input_file_frame = Col([
        [Text("Sample spectra", size=(12, 1), tooltip="Specify absolute path to sample spectra file (.mgf)"),
         Input(key="SamplePath", default_text=SamplePath), FileBrowse()],
        [Text("Parameter", size=(12, 1), tooltip="Specify absolute path to parameter file (.par)"),
         Input(key="ParametersFile", default_text=ParametersFile), FileBrowse()]])
    input_dir_frame = Col(
        [[Text("Sample data", size=(12, 1)), Input(key="DataDir", default_text=DataDir), FolderBrowse()],
         [Text("Database", size=(12, 1)), Input(key="DatabaseDir", default_text=DatabaseDir), FolderBrowse()],
         [Text("Peptide Shaker", size=(12, 1)), Input(key="PeptideShakerDir", default_text=PeptideShakerDir),
          FolderBrowse()],
         [Text("Search GUI", size=(12, 1)), Input(key="SearchGUIDir", default_text=SearchGUIDir), FolderBrowse()],
         [Text("Resources", size=(12, 1)), Input(key="ResourcesDir", default_text="resources/")],
         [Text("Results", size=(12, 1)), Input(key="ResultsDir", default_text="results/")],
         [Text("TaxID mapping", size=(12, 1)),
          Input(key="TaxidMapping", default_text="taxidMapping/")]])
    # pepgm gridsearch subframe
    gridsearch_frame = Col([[Text("alpha", size=(5, 1)),
                             Input(default_text=Alpha, key="Alpha")],
                            [Text("beta", size=(5, 1)),
                             Input(default_text=Beta, key="Beta")],
                            [Text("prior", size=(5, 1)),
                             Input(default_text=prior, key="prior")]], expand_x=True, expand_y=True)
    # pepgm plotting subframe
    plotting_frame = Col([[Text("Number of taxa", size=(15, 1)),
                           Input(default_text=TaxaInPlot, key="TaxaInPlot")],
                          [Text("Number of taxa \n in protein count", size=(15, 1)),
                           Input(default_text=TaxaInProteinCount, key="TaxaInProteinCount")],
                          [Text("Source database", size=(15, 1)),
                           Input(default_text=sourceDB, key="sourceDB")]], expand_x=True, expand_y=True)
    # pepgm frame
    pepgm_frame = Col([[Frame("Grid search", [[gridsearch_frame]]),
                        Frame("Results plotting", [[plotting_frame]], expand_x=True, expand_y=True)]])
    # search gui fdr subframe
    fdr_frame = Col([[Text("PSM", size=(10, 1)),
                      Slider(range=(0, 20), default_value=psmFDR, resolution=1, orientation='horizontal',
                             key="psmFDR")],
                     [Text("Peptide", size=(10, 1)),
                      Slider(range=(0, 20), default_value=peptideFDR, resolution=1, orientation='horizontal',
                             key="peptideFDR")],
                     [Text("Protein", size=(10, 1)),
                      Slider(range=(0, 20), default_value=proteinFDR, resolution=1, orientation='horizontal',
                             key="proteinFDR")]], expand_x=True, expand_y=True)
    # searchgui frame
    searchgui_frame = Col([[Text("Search engine", size=(15, 1)),
                            Combo(list(
                                ["-xtandem", "-myrimatch", "-ms_amanda", "-msgf", "-omssa", "-comet", "-tide",
                                 "-andromeda", "-meta_morpheus", "-novor", "-directtag"]), readonly=True,
                                key="searchengines", default_value="-xtandem")],
                           [Frame("FDR calculation", [[fdr_frame]])]], expand_x=True, expand_y=True)

    return config_frame, run_frame, input_file_frame, input_dir_frame, gridsearch_frame, plotting_frame, pepgm_frame, \
           fdr_frame, searchgui_frame
