
import argparse
from FactorGraphGeneration import *



parser = argparse.ArgumentParser(description = 'Run the PepGM algorithm from command line')

parser.add_argument('--targetTaxa', nargs = '*', help ='enter a list of taxa to include in your graphical model')
parser.add_argument('--PSM_Report', type =str, required =True, help = 'path to your PSM report txt file (output from peptideshaker)')
parser.add_argument('--PeptideMapPath',type=str, required =True, help = 'path to where you want to save you taxon-peptide map .json file')#make it so that this works as argument for both functions
parser.add_argument('--out', type = str, required = True, help = 'path to where you want to save the GraphML file of the factorgraph')

args = parser.parse_args()

Taxongraph = TaxonGraph()
Taxongraph.GetAllLeafTaxa(args.targetTaxa)
#only fetch TaxonData if the PeptideMap File doesn't exist yet:
if not exists(args.PeptideMapPath):
    Taxongraph.FetchTaxonData(args.PeptideMapPath)
Taxongraph.CreateTaxonPeptidegraphFromMzID(args.PSM_Report,args.PeptideMapPath,0.001)
#Taxongraph.CreateExample()
Factorgraph = FactorGraph()
Factorgraph.ConstructFromTaxonGraph(Taxongraph)
CTFactorgraph = GenerateCTFactorGraphs(Factorgraph)
CTFactorgraph.SaveToGraphML(args.out)
