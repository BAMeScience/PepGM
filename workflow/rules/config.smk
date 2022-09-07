configfile: 'config/config.yaml'


EntrezMail = config['APImail']
EntrezKey = config['APIkey']

DataDirectory = config['DataDir']
DatabaseDirectory = config ['DatabaseDir']
ResultsDir = config['ResultsDir'] + config['ExperimentName'] +'/'
ResourcesDir = config['ResourcesDir']
TaxidMapping = config['TaxidMapping']
ResultsDirStrain = config['ResultsDir']
SamplePath = config['SamplePath']
ParametersFile = config['ParametersFile']


PeptideShakerDir = config['PeptideShakerDir']
SearchGUIDir = config['SearchGUIDir']

searchengines = config['searchengines']
peptideFDR = config['peptideFDR']
proteinFDR = config['proteinFDR']
psmFDR = config['psmFDR']


FilterSpectra = config['FilterSpectra']
AddHostandCrapToDB = config['AddHostandCrapToDB']
ExperimentName = config['ExperimentName']
SampleName = config['SampleName']
HostName = config['HostName']
ScientificHostName = config['ScientificHostName']
ReferenceDBName = config['ReferenceDBName']

TaxaInPlot = config['TaxaInPlot']
TaxaInProteinCount = config['TaxaInProteinCount']
sourceDB = config['sourceDB']

AlphaRange = config['Alpha']
BetaRange = config['Beta']
prior = config['prior']
