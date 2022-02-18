
import argparse
from os.path import exists
from FactorGraphGeneration import *



parser = argparse.ArgumentParser(description = 'Run the PepGM algorithm from command line')

parser.add_argument('--targetTaxa', nargs = '*', help ='enter a list of taxa to include in your graphical model')
parser.add_argument('--PSM_Report', type =str, required =True, help = 'path to your PSM report txt file (output from peptideshaker)')
parser.add_argument('--PeptideMapPath',type=str, required =True, help = 'path to where you want to save you taxon-peptide map .json file')#make it so that this works as argument for both functions
parser.add_argument('--out', type = str, required = True, help = 'path to where you want to save the GraphML file of the factorgraph')
parser.add_argument('--sourceDB',type = str, nargs ='?', const ='', help = 'name of the DB queried through Entrez')

args = parser.parse_args()

#targetTaxa = ['coronaviridae','chlorocebus']
#PeptideMapPath = '/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD025130_Sars_CoV_2/coronaviridae.json'
#PSM_Report = '/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD025130_Sars_CoV_2/chlorocebus_refseq_Default_PSM_Report.txt'
#out = '/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD025130_Sars_CoV_2/chlorocebus_refseq_PepGM_graph.graphml'

#targetTaxa = ['herpesviridae','homo sapiens']
#PeptideMapPath = '/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD005104_Herpessimplex_1/herpesviridae.json'
#PSM_Report = '/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD005104_Herpessimplex_1/human_refseq_Default_PSM_Report.txt'
#out = '/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD005104_Herpessimplex_1/human_refseq_PepGM_graph.graphml'



Taxongraph = TaxonGraph()
Taxongraph.GetAllLeafTaxa(args.targetTaxa)
#only fetch TaxonData if the PeptideMap File doesn't exist yet:
if not exists(args.PeptideMapPath):
    print(exists(args.PeptideMapPath))
    Taxongraph.FetchTaxonData(args.PeptideMapPath, args.sourceDB)
Taxongraph.CreateTaxonPeptidegraphFromPSMresults(args.PSM_Report,args.PeptideMapPath,0.001)
#Taxongraph.CreateExample()
Factorgraph = FactorGraph()
Factorgraph.ConstructFromTaxonGraph(Taxongraph)
CTFactorgraph = GenerateCTFactorGraphs(Factorgraph)
CTFactorgraph.SaveToGraphML(args.out)
