
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

#if __name__ =='__main__':
#    UnipeptCSV = '/home/tholstei/repos/PepGM_all/PepGM/results/SIHUMIx_no_subspecies/CAMPI_SIHUMIx/differentGraphs300.csv'
#    Taxongraph = TaxonGraph()
#    Taxongraph.CreateFromUnipeptResponseCSV(UnipeptCSV)
#    Factorgraph = FactorGraph()
#    Factorgraph.ConstructFromTaxonGraph(Taxongraph)
#    CTFactorgraph = GenerateCTFactorGraphs(Factorgraph)
#    CTFactorgraph.SaveToGraphML('/home/tholstei/repos/PepGM_all/PepGM/results/SIHUMIx_no_subspecies/CAMPI_SIHUMIx/differentGraphs300.graphml')