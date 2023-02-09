configfile: 'config/config.yaml'

DataDirectory = config['DataDir']
DatabaseDirectory = config ['DatabaseDir']
ResultsDir = config['ResultsDir'] + config['ExperimentName'] +'/'+config['SampleName']+'/'
XTandemDir = config['ResultsDir'] + config['ExperimentName'] +'/'+config['SampleName']+'/XTandem/'
MS2RescoreDir = config['ResultsDir'] + config['ExperimentName'] +'/'+config['SampleName']+'/MS2Rescore/'
ResourcesDir = config['ResourcesDir']
TaxidMapping = config['TaxidMapping']
ResultsDirStrain = config['ResultsDir']

FilterSpectra = config['FilterSpectra']
AddHostandCrapToDB = config['AddHostandCrapToDB']
ExperimentName = config['ExperimentName']
SpectraFileType = config['SpectraFileType']
SampleName = config['SampleName']
HostName = config['HostName']
ScientificHostName = config['ScientificHostName']
ReferenceDBName = config['ReferenceDBName']


TaxaInPlot = config['TaxaInPlot']
#TaxaInProteinCount = config['TaxaInProteinCount']
#sourceDB = config['sourceDB']

AlphaRange = config['Alpha']
BetaRange = config['Beta']
prior = config['prior']

# X!Tandem parameters
XTANDEM_DEFAULT = config['xtandem_default']
XTANDEM_FMME = config['xtandem_fmme']
XTANDEM_FMMEU = config['xtandem_fmmeu']
XTANDEM_PMMEP = config['xtandem_pmmep']
XTANDEM_PMMEM = config['xtandem_pmmem']
XTANDEM_PMMEU = config['xtandem_pmmeu']
# X!Tandem PTMs (comma separated if more than one)
XTANDEM_MODS_FIX = config['xtandem_mods_fixed'] 
XTANDEM_MODS_VAR = config['xtandem_mods_variable'] 
XTANDEM_MODS_VAR_NTERM = config['xtandem_mods_variable_nterm'] 
# additional parameters (dict: param name -> value)
XTANDEM_PARAS = config['xtandem_add_params']

#MS2rescore parameters
RescorePipeline = config['RescorePipeline']
RescoreFeatures = config['RescoreFeatures']
RunPercolator = config['RunPercolator']
FragModel = config['FragModel']
mods = config['Mods']

#UnipeptQueryParameter
TaxaNumber = config['TaxaNumber']
targetTaxa = config['targetTaxa']
FDR = config['FDR']

#when using searchGUI/PepShaker
PeptideShaker = config['PeptideShaker']
SearchGUI = config['SearchGUI']

searchengines = config['searchengines']
peptideFDR = config['peptideFDR']
proteinFDR = config['proteinFDR']
psmFDR = config['psmFDR']

SamplePath = config['SamplePath']
ParametersFile = config['ParametersFile']