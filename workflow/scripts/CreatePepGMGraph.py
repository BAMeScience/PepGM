"""
This module initializes a taxon graph and fills in peptides and factors.
"""
import argparse
from os.path import exists
from FactorGraphGeneration import *

# argparser preliminaries
parser = argparse.ArgumentParser(description = 'Initialize graphical model with taxons, peptides and factors.')
parser.add_argument('--targetTaxa', type = str, nargs = 1,
                    help ='Path to mapped taxids text file. Format: One taxid per line.')
parser.add_argument('--PSM_Report', type =str, required =True,
                    help = 'Path to psm report (peptide shaker ouput).')
parser.add_argument('--PeptideMapPath',type=str, required =True,
                    help = 'Path to taxon peptide map. Format: json file')
parser.add_argument('--out', type = str, required = True,
                    help = "Output directory and filename.")
args = parser.parse_args()

# init taxon networkx graph object
Taxongraph = TaxonGraph()
Taxongraph.GetAllLeafTaxaFromTaxids(args.targetTaxa[0])

# fetch taxon data if the peptide map file does not exist yet
if not exists(args.PeptideMapPath):
    print(exists(args.PeptideMapPath))
    Taxongraph.FetchTaxonData(args.PeptideMapPath, args.sourceDB)

# add peptides from peptide shaker output
Taxongraph.CreateTaxonPeptidegraphFromPSMresults(args.PeptideMapPath, args.PSM_Report, min_score = 0.001)

# add factors
Factorgraph = FactorGraph()
Factorgraph.ConstructFromTaxonGraph(Taxongraph)
CTFactorgraph = GenerateCTFactorGraphs(Factorgraph)
# save
CTFactorgraph.SaveToGraphML(args.out)
