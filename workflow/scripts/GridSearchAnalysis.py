import numpy as np
import pandas as pd
import os
import re 
import matplotlib
from matplotlib import pyplot as plt
from ete3 import NCBITaxa
from  scipy.stats import entropy


#package s ettings
matplotlib.use('Agg')
ncbi = NCBITaxa()



def ComputeMetric(resultsfolder, host, output, weightsfile):

    """
    Compute a "goodness" metric for the parameters used in the grid search. 
    Generates a barplot of the metric and returns the best identified parameters.
    Return a list of the best parameter set in order [alpha,beta,gamma].

    :param resultsfolder: str, path to PepGM resultsfolder
    :param host: str, host to be excluded from the parameter checked taxa
    :param output: str, name of the output .png file 
    :param weightsfile: str, path to the .csv file with Unipept output
    
    """
    
    #predefine necessary lists
    Params = []
    Metrics = []
    SumProportions = []
    Entropies = []
    TaxDistances = []
    WeightCoeffs = []

    #file with weights of taxids
    Weights = pd.read_csv(weightsfile,usecols=['taxa','weight'])
    Weights = Weights.groupby(['taxa']).sum().reset_index()
    #drop host taxids
    HostTaxid = ncbi.get_name_translator([host])[host][0]
    HostTaxidList = [str(i) for i in ncbi.get_descendant_taxa(HostTaxid)]+[str(HostTaxid)]
    Weights = Weights[Weights.taxa.isin(HostTaxidList)==False]
    Weights.drop(Weights[Weights.taxa.isin(HostTaxidList)].index, inplace=True)
    Maxweight = Weights.max()['weight']
    AllTaxidsToAdd = []
    #add descendant taxa into weight dataframe
    for Taxid in Weights['taxa']:
        TaxidsToAdd = ncbi.get_descendant_taxa(Taxid)
        TaxidsToAdd = [[txd,float((Weights.loc[Weights['taxa']==Taxid]['weight']).to_string(index=False))]for txd in TaxidsToAdd]
        AllTaxidsToAdd.append(TaxidsToAdd[:])
       
    AllTaxidsToAdd = [txd_1 for txd in AllTaxidsToAdd for txd_1 in txd]
    AllWeights = pd.concat([Weights, pd.DataFrame(AllTaxidsToAdd,columns = ['taxa','weight'])],axis =0)

        

    
    for folders in os.listdir(resultsfolder):
        if os.path.isdir(resultsfolder+ '/' + folders):
            for file in os.listdir(resultsfolder+ '/' + folders):
    
                if file.endswith('.csv'):
            
                    Results = pd.read_csv(resultsfolder+ '/' + folders +'/' + file,names = ['ID','score','type'])
                    TaxIDS = Results.loc[Results['type']=='taxon']
                    TaxIDS.loc[:,'score'] = pd.to_numeric(TaxIDS['score'],downcast = 'float')
                    TaxIDS = TaxIDS.sort_values('score', ascending = False)
                    TaxIDS = TaxIDS[TaxIDS.ID.isin(HostTaxidList)==False]
                    TaxIDS.drop(TaxIDS[TaxIDS.ID.isin(HostTaxidList)].index, inplace=True)
                    
                    #compute the metric
                    
                    #what's the weight of the highest scoring taxid?
                    Weight = AllWeights.loc[AllWeights['taxa']==int(TaxIDS.ID.head(1).item())]['weight'].head(1).item()
                    WeightCoeff = Weight/Maxweight
                    WeightCoeffs.append(WeightCoeff)
    
                    #compare the posterior probbilities
                    FirstSum = np.sum(TaxIDS[:3]['score'])
                    NextSum = np.sum(TaxIDS[3:8]['score'])
                    SumProportion = FirstSum/NextSum
                    SumProportions.append(SumProportion)
                    Entropy = entropy(TaxIDS[:10]['score'])

                  
    
                    #compute the pairwise taxonomic distance of the first 4 Taxa
                    FourTaxa = np.array(TaxIDS[:4]['ID'])
                    tree = ncbi.get_topology(FourTaxa)
                    DistanceSum = (tree.get_distance(FourTaxa[0],FourTaxa[1]))**2+tree.get_distance(FourTaxa[1],FourTaxa[2])+tree.get_distance(FourTaxa[2],FourTaxa[3])
                    TaxDistances.append(DistanceSum)
                    Matching = (1/(Entropy))*SumProportion*(1/DistanceSum)*WeightCoeff
                    checkpoint= 3
                    
                    if not np.isnan(Matching) and Matching < 3000:


                        Metrics.append(Matching)
    
                        #get the corresponding parameters
                        params = re.findall('\d+\.\d+',file)
                        Params.append([params[0],params[1],params[2]])

     
    zipLists = zip(Metrics,Params)
    SortedPairs = sorted(zipLists,reverse= True)

    
    tuples = zip(*SortedPairs)
    Metrics,Params = [list(tuple) for tuple in tuples]
    
    figure1 = plt.figure(figsize=(30,15), tight_layout=True)
    
    plt.barh(list(range(20)),Metrics[0:20],color = 'mediumvioletred')
    plt.xticks(fontsize =35)
    plt.yticks(list(range(20)),Params[0:20],fontsize=35)
    
    print(output)
    plt.savefig(output)
    plt.close()

    return Params[0]



if __name__=='__main__':
    resultsfolder = '/home/tholstei/repos/PepGM_all/PepGM/results/SearchGUI_full_Uniprot/PXD003013_Cowpox_BR'
    host = 'homo sapiens'
    output = '/home/tholstei/repos/PepGM_all/PepGM/results/SearchGUI_full_Uniprot/PXD003013_Cowpox_BR/paramchecktest.png'
    weightsfile = '/home/tholstei/repos/PepGM_all/PepGM/results/SearchGUI_full_Uniprot/PXD003013_Cowpox_BR/GraphDataframe.csv'
    ComputeMetric(resultsfolder, host, output, weightsfile)
    




