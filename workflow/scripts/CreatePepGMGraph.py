"""
This module initializes the graphical models and adds taxons, peptides and factors.
"""
import argparse
from os.path import exists
from FactorGraphGeneration import *


# argparse preliminaries
parser = argparse.ArgumentParser(description = 'Run the PepGM algorithm from command line')
parser.add_argument('--targetTaxa', type = str, nargs = 1, help ='path to the .txt file containing a list of taids to include in the graph, one taxid per line')
parser.add_argument('--PSM_Report', type =str, required =True, help = 'path to your PSM report txt file (output from peptideshaker)')
parser.add_argument('--PeptideMapPath',type=str, required =True, help = 'path to where you want to save you taxon-peptide map .json file')#make it so that this works as argument for both functions
parser.add_argument('--out', type = str, required = True, help = 'path to where you want to save the GraphML file of the factorgraph')
parser.add_argument('--sourceDB',type = str, nargs ='?', const ='', help = 'name of the DB queried through Entrez')
args = parser.parse_args()

# init networkx graph object
Taxongraph = TaxonGraph()
Taxongraph.GetAllLeafTaxaFromTaxids(args.targetTaxa[0])
# only fetch taxon data if the taxon peptide map file doesn't exist yet
if not exists(args.PeptideMapPath):
    print(exists(args.PeptideMapPath))
    Taxongraph.FetchTaxonData(args.PeptideMapPath, args.sourceDB)
# add peptides
Taxongraph.CreateTaxonPeptidegraphFromPSMresults(args.PeptideMapPath, args.PSM_Report, 0.001)

# add factors
Factorgraph = FactorGraph()
Factorgraph.ConstructFromTaxonGraph(Taxongraph)
CTFactorgraph = GenerateCTFactorGraphs(Factorgraph)
CTFactorgraph.SaveToGraphML(args.out)
