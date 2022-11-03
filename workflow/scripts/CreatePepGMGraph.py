
import argparse
from os.path import exists
from FactorGraphGeneration import *



parser = argparse.ArgumentParser(description = 'Run the PepGM algorithm from command line')

parser.add_argument('--UnipeptCSV', type = str, required = True, help = 'path to where you want to save the GraphML file of the factorgraph')
parser.add_argument('--out', type = str, required = True, help = 'path to output file where graphml will be saved')

args = parser.parse_args()


Taxongraph = TaxonGraph()
Taxongraph.CreateFromUnipeptResponseCSV(args.UnipeptCSV)
Factorgraph = FactorGraph()
Factorgraph.ConstructFromTaxonGraph(Taxongraph)
CTFactorgraph = GenerateCTFactorGraphs(Factorgraph)
CTFactorgraph.SaveToGraphML(args.out)
