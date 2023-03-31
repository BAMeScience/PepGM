import argparse
from zero_lookahead_belief_propagation import *
from FactorGraphGeneration import * 



#parser = argparse.ArgumentParser(description = 'Run the PepGM algorithm from command line')

#parser.add_argument('--GraphMLPath', type = str, required = True, help = 'path to where you want to save the GraphML file of the factorgraph')
#parser.add_argument('--max_iter', nargs = '?' ,type = int , default = 10000, help ='max. number of iterations the belief propagation algo will go through')
#parser.add_argument('--tol', nargs = '?' ,type = float, default = 0.006, help = 'residual error allowed for the BP algorithm')
#parser.add_argument('--out', type = str, required= True, help = 'path to the file you want to save your results as')
#parser.add_argument('--alpha', type = float, required = True, help ='detection probability of a peptide for the noisy-OR model')
#parser.add_argument('--beta',type = float, required = True, help = 'probability of wrong detection')
#parser.add_argument('--prior', type = float, required = True, help = 'prior assigned to all taxa')

#args = parser.parse_args()

#CTFactorgraph = CTFactorGraph(args.GraphMLPath)
#CTFactorgraph.FillInFactors(args.alpha,args.beta)
#CTFactorgraph.FillInPriors(args.prior)
#CTFactorgraph.AddCTNodes()



#CTFactorgraphs = [SeparateSubgraphs(CTFactorgraph,filternodes) for filternodes in nx.connected_components(CTFactorgraph)]

#Resultlist,Resultsdict,Nodetypes = CalibrateAllSubgraphs(CTFactorgraphs,args.max_iter,args.tol)
#save = SaveResultsToCsv(Resultsdict,Nodetypes,args.out)
#VisualizeResults(Nodetypes,Resultsdict,'taxon',Taxongraph)

if __name__=='__main__':

    GraphMLPath = '/home/tholstei/repos/PepGM_all/PepGM/results/SIHUMIx_CT_Underflow/CAMPI_SIHUMIx/SIHUMI_DB1UNIPROT_PepGM_graph.graphml'
    alpha = 0.8
    beta = 0.4
    prior = 0.5
    max_iter = 10000
    tol = 0.05
    CTFactorgraph = CTFactorGraph(GraphMLPath)
    CTFactorgraph.FillInFactors(alpha,beta)
    CTFactorgraph.FillInPriors(prior)
    CTFactorgraph.AddCTNodes()
    out = 'test.csv'



    CTFactorgraphs = [SeparateSubgraphs(CTFactorgraph,filternodes) for filternodes in nx.connected_components(CTFactorgraph)]

    Resultlist,Resultsdict,Nodetypes = CalibrateAllSubgraphs(CTFactorgraphs,max_iter,tol)
    save = SaveResultsToCsv(Resultsdict,Nodetypes,out)