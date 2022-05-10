from GridSearchAnalysis import *
from PhyloTreeView import *
import pandas as pd
import os
import argparse 
from ete3 import NCBITaxa

parser = argparse.ArgumentParser(description = 'Downstream analysis of the grid search for the PepGM algorithm')
parser.add_argument('--resultsfolder', required = True, help = 'folder with the results from the beliefpropagation in csv format')
parser.add_argument('--out', required=True, help ='output png file')
parser.add_argument('--host', required =True, help = 'name of the host, to be excluded from parameter checked taxa')
parser.add_argument('--reference_db',required = True, help = 'name of the reference datasbe used to include in the filepath string')

args = parser.parse_args()

ncbi = NCBITaxa()


def SaveReducedCSV(filepath,host,output):
    '''
    Saves CSV of the 20 best results to a file in the main results folder
    :param filepath: str, path of cvs to be reduced
    :param host: str,if included in the search results, the host will be excluded from the results list
    :param output: str,output path of csv
    '''
        
    HostTaxid = ncbi.get_name_translator([host])[host][0]
    HostTaxidList = [str(i) for i in ncbi.get_descendant_taxa(HostTaxid)]+[str(HostTaxid)]
            
    Results = pd.read_csv(filepath,names = ['ID','score','type'])
    TaxIDS = Results.loc[Results['type']=='taxon']
    TaxIDS.loc[:,'score'] = pd.to_numeric(TaxIDS['score'],downcast = 'float')
    TaxIDS = TaxIDS.sort_values('score', ascending = False)
    TaxIDS = TaxIDS[TaxIDS.ID.isin(HostTaxidList)==False]
    TaxIDS.drop(TaxIDS[TaxIDS.ID.isin(HostTaxidList)].index, inplace=True)
    
    TaxIDS.to_csv(output)
           


def MoveBestResultsPlot(filepath,out):
    '''
    Copies barplot of best identifiedparameter set to main resultsfolder
    :param filepath: list, the three grid search parameters identified as best for the sample at hand
    :param out: output path where the file will be copied

    '''
    os.system('cp '+filepath +' '+out)


#analyse the PepGM grid search with an empirical metric, retrieve best parameters
Parameters = ComputeMetric(args.resultsfolder, args.host, args.out)
Resultsfile = args.resultsfolder+'/Prior'+str(Parameters[2]) +'/'+ args.reference_db + '_PepGM_Results_a'+str(Parameters[0])+'_b'+str(Parameters[1])+'_p'+str(Parameters[2])

#save the reduced results csv
SaveReducedCSV(Resultsfile+'.csv', args.host, args.resultsfolder +'/PepGm_Results.csv')
MoveBestResultsPlot(Resultsfile+'.png',args.resultsfolder+'/PepGM_ResultsPlot.png')

#save a phylogenetic tree view of the PepGm results
CreatePhyloTreeView(Resultsfile+'.csv',args.host, args.resultsfolder + '/PhyloTreeView.png')