import argparse
from belief_propagation import *
from FactorGraphGeneration import * 
from os.path import exists


parser = argparse.ArgumentParser(description = 'Run the PepGM algorithm from command line')

parser.add_argument('--GraphMLPath', type = str, required = True, help = 'path to where you want to save the GraphML file of the factorgraph')
parser.add_argument('--max_iter', nargs = '?' ,type = int , default = 10000, help ='max. number of iterations the belief propagation algo will go through')
parser.add_argument('--tol', nargs = '?' ,type = float, default = 0.006, help = 'residual error allowed for the BP algorithm')
parser.add_argument('--out', type = str, required= True, help = 'path to the file you want to save your results as')
parser.add_argument('--alpha', type = float, required = True, help ='detection probability of a peptide for the noisy-OR model')
parser.add_argument('--beta',type = float, required = True, help = 'probability of wrong detection')

args = parser.parse_args()


#targetTaxa = ['igacovirus','gallus gallus']
#PeptideMapPath = '/home/tholstei/repos/VirusGraph/Data/SampleData/PXD002936_avian_bronchitis/igacovirus.json'
#PSM_Report = '/home/tholstei/repos/VirusGraph/Data/SampleData/PXD002936_avian_bronchitis/chicken_refseq_Default_PSM_Report.txt'
#out = '/home/tholstei/repos/VirusGraph/Data/SampleData/PXD002936_avian_bronchitis/chicken_refseq_PepGMresults.csv'
#max_iter=1000
#tol=0.006



CTFactorgraph = nx.read_graphml(args.GraphMLPath)
CTFactorgraph.FillInFactors(args.alpha,args.beta)

CTFactorgraphs = SeparateSubgraphs(CTFactorGraph)

#TODO introduce the grid search here
Resultlist,Resultsdict,Nodetypes = CalibrateAllSubgraphs(CTFactorgraphs,args.max_iter,args.tol)
save = SaveResultsToCsv(Resultsdict,Nodetypes,args.out)
#VisualizeResults(Nodetypes,Resultsdict,'taxon',Taxongraph)