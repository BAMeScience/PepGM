
import argparse
from belief_propagation import *
from os.path import exists


parser = argparse.ArgumentParser(description = 'Run the PepGM algorithm from command line')

parser.add_argument('--targetTaxa', nargs = '*', help ='enter a list of taxa to include in your graphical model')
parser.add_argument('--PeptideMapPath',type=str, required =True, help = 'path to where you want to save you taxon-peptide map .json file')#make it so that this works as argument for both functions
parser.add_argument('--PSM_Report', type =str, required =True, help = 'path to your PSM report txt file (output from peptideshaker)')
parser.add_argument('--max_iter', nargs = '?' ,type = int , default = 10000, help ='max. number of iterations the belief propagation algo will go through')
parser.add_argument('--tol', nargs = '?' ,type = float, default = 0.006, help = 'residual error allowed for the BP algorithm')
parser.add_argument('--out', type = str, help = 'path to the file you want to save your results as')

args = parser.parse_args()


#targetTaxa = ['igacovirus','gallus gallus']
#PeptideMapPath = '/home/tholstei/repos/VirusGraph/Data/SampleData/PXD002936_avian_bronchitis/igacovirus.json'
#PSM_Report = '/home/tholstei/repos/VirusGraph/Data/SampleData/PXD002936_avian_bronchitis/chicken_refseq_Default_PSM_Report.txt'
#out = '/home/tholstei/repos/VirusGraph/Data/SampleData/PXD002936_avian_bronchitis/chicken_refseq_PepGMresults.csv'
#max_iter=1000
#tol=0.006


Taxongraph = TaxonGraph()
Taxongraph.GetAllLeafTaxa(args.targetTaxa)
#only fetch TaxonData if the PeptideMap File doesn't exist yet:
if not exists(args.PeptideMapPath):
    Taxongraph.FetchTaxonData(args.PeptideMapPath)
Taxongraph.CreateTaxonPeptidegraphFromMzID(args.PSM_Report,args.PeptideMapPath,0.001)
#Taxongraph.CreateExample()
Factorgraph = FactorGraph()
Factorgraph.ConstructFromTaxonGraph(Taxongraph)


#To save a gml of your graph, uncomment this
#GmlCTFactorgraphs = CTFactorGraph(Factorgraph)
#GmlCTFactorgraphs.SaveToGml('TaxonpeptideGraph_herpesviridae_pxd005104')


Factorgraphs = Factorgraph.SeparateSubgraphs()

CTFactorgraphs = GenerateCTFactorGraphs(Factorgraphs)

#TODO introduce the grid search here
Resultlist,Resultsdict,Nodetypes = CalibrateAllSubgraphs(CTFactorgraphs,args.max_iter,args.tol)
save = SaveResultsToCsv(Resultsdict,Nodetypes,args.out)
#VisualizeResults(Nodetypes,Resultsdict,'taxon',Taxongraph)