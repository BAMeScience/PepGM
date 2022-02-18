import argparse
from belief_propagation_log import *
from FactorGraphGeneration import * 
from os.path import exists


parser = argparse.ArgumentParser(description = 'Run the PepGM algorithm from command line')

parser.add_argument('--GraphMLPath', type = str, required = True, help = 'path to where you want to save the GraphML file of the factorgraph')
parser.add_argument('--max_iter', nargs = '?' ,type = int , default = 10000, help ='max. number of iterations the belief propagation algo will go through')
parser.add_argument('--tol', nargs = '?' ,type = float, default = 0.006, help = 'residual error allowed for the BP algorithm')
parser.add_argument('--out', type = str, required= True, help = 'path to the file you want to save your results as')
parser.add_argument('--alpha', type = float, required = True, help ='detection probability of a peptide for the noisy-OR model')
parser.add_argument('--beta',type = float, required = True, help = 'probability of wrong detection')
parser.add_argument('--prior', type = float, required = True, help = 'prior assigned to all taxa')

args = parser.parse_args()

#graphMLPath = '/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD005104_Herpessimplex_1/human_refseq_PepGM_graph.graphml'
#alpha = 0.9
#beta = 0.3
#out =  'resources/SampleData/PXD005104_Herpessimplex_1/human_refseq_PepGM_Results_0.8_0.4.csv'
#max_iter = 1000
#tol = 0.003



CTFactorgraph = CTFactorGraph(args.GraphMLPath)
CTFactorgraph.FillInFactors(args.alpha,args.beta)
CTFactorgraph.FillInPriors(args.prior)



CTFactorgraphs = [SeparateSubgraphs(CTFactorgraph,filternodes) for filternodes in nx.connected_components(CTFactorgraph)]

Resultlist,Resultsdict,Nodetypes = CalibrateAllSubgraphs(CTFactorgraphs,args.max_iter,args.tol)
save = SaveResultsToCsv(Resultsdict,Nodetypes,args.out)
#VisualizeResults(Nodetypes,Resultsdict,'taxon',Taxongraph)