""" Build GUI layout. """
from PySimpleGUI import Col, Input, Frame, FolderBrowse, FileBrowse, Slider, InputText, Text, Combo, Radio, theme

theme("SystemDefaultForReal")


def build(ConfigName, ExperimentName, SampleName, HostName, ScientificHostName, ReferenceDBName, SamplePath,
          ParametersFile, DataDir, DatabaseDir, PeptideShakerDir, SearchGUIDir, Alpha, Beta, prior, psmFDR,
          peptideFDR, proteinFDR, TaxaInPlot, TaxaInProteinCount, sourceDB):
    config_frame = Col([[InputText(size=(70, 5), key="config_file_name", default_text=ConfigName)]], size=(600, 30),
                       pad=(0, 0))

    run_frame = Col([[InputText(size=(70, 5), key="ExperimentName", default_text=ExperimentName)],
                     [Text("Sample", size=(15, 1)), Input(size=(43, 5), key="SampleName", default_text=SampleName)],
                     [Text("Host", size=(15, 1)), Input(size=(43, 5), key="HostName", default_text=HostName)],
                     [Text("Scientific host", size=(15, 1)),
                      Input(size=(43, 5), key="ScientificHostName", default_text=ScientificHostName)],
                     [Text("Reference", size=(15, 1)),
                      Input(size=(43, 5), key="ReferenceDBName", default_text=ReferenceDBName)],
                     [Radio('Filter spectra', default=False, key="FilterSpectra", group_id=1),
                      Radio('Add host and crap databases', default=True, key="AddHostandCrapToDB", group_id=1)]],
                    size=(600, 150), pad=(0, 0))

    input_file_frame = Col(
        [[Text("Sample", size=(15, 1)), Input(key="SamplePath", default_text=SamplePath), FileBrowse()],
         [Text("Parameter", size=(15, 1)), Input(key="ParametersFile", default_text=ParametersFile), FileBrowse()]],
        size=(600, 70),
        pad=(0, 0))

    input_dir_frame = Col(
        [[Text("Sample data", size=(15, 1)), Input(key="DataDir", default_text=DataDir), FolderBrowse()],
         [Text("Database", size=(15, 1)), Input(key="DatabaseDir", default_text=DatabaseDir), FolderBrowse()],
         [Text("Peptide Shaker", size=(15, 1)), Input(key="PeptideShakerDir", default_text=PeptideShakerDir),
          FolderBrowse()],
         [Text("Search GUI", size=(15, 1)), Input(key="SearchGUIDir", default_text=SearchGUIDir), FolderBrowse()],
         [Text("Resources", size=(15, 1)), Input(key="ResourcesDir", default_text="resources/")],
         [Text("Results", size=(15, 1)), Input(key="ResultsDir", default_text="results/")],
         [Text("TaxID mapping", size=(15, 1)),
          Input(key="TaxidMapping", default_text="taxidMapping/")]],
        size=(600, 220),
        pad=(0, 0))

    gridsearch_frame = Col([[Text("alpha", size=(5, 1)),
                             Input(default_text=Alpha, key="Alpha", size=(30, 5))],
                            [Text("beta", size=(5, 1)),
                             Input(default_text=Beta, key="Beta", size=(30, 5))],
                            [Text("prior", size=(5, 1)),
                             Input(default_text=prior, key="prior", size=(30, 5))]])

    plotting_frame = Col([[Text("Number of taxa", size=(15, 1)),
                           Input(default_text=TaxaInPlot, key="TaxaInPlot", size=(7, 1))],
                          [Text("Number of taxa in protein count", size=(15, 1)),
                           Input(default_text=TaxaInProteinCount, key="TaxaInProteinCount", size=(7, 1))],
                          [Text("Source database", size=(15, 1)),
                           Input(default_text=sourceDB, key="sourceDB", size=(7, 1))]], size=(230, 70))

    pepgm_frame = Col([[Frame("Grid search", [[gridsearch_frame]]),
                        Frame("Results plotting", [[plotting_frame]])]],
                      size=(600, 120), pad=(0, 0))

    fdr_frame = Col([[Text("PSM"), Slider(range=(0, 20), default_value=psmFDR, resolution=1,
                                          orientation='horizontal', key="psmFDR", size=(15, 10)),
                      Text("Peptide"),
                      Slider(range=(0, 20), default_value=peptideFDR, resolution=1,
                             orientation='horizontal', key="peptideFDR", size=(15, 10)),
                      Text("Protein"),
                      Slider(range=(0, 20), default_value=proteinFDR, resolution=1,
                             orientation='horizontal', key="proteinFDR", size=(15, 10))]], size=(580, 50), pad=(0, 0))

    searchgui_frame = Col([[Text("Search engine", size=(15, 1)),
                            Combo(list(
                                ["-xtandem", "-myrimatch", "-ms_amanda", "-msgf", "-omssa", "-comet", "-tide",
                                 "-andromeda", "-meta_morpheus", "-novor", "directtag"]), readonly=True,
                                key="searchengines", default_value="-xtandem")],
                           [Frame("FDR settings", [[fdr_frame]])]], size=(600, 100), pad=(0, 0))

    return config_frame, run_frame, input_file_frame, input_dir_frame, gridsearch_frame, plotting_frame, pepgm_frame, \
           fdr_frame, searchgui_frame
