import numpy as np
import pandas as pd
import os
import argparse 
import re 
import matplotlib
from matplotlib import pyplot as plt
from ete3 import NCBITaxa
from ete3 import Tree

matplotlib.use('Agg')
ncbi = NCBITaxa()
# script to compute a "goodness" metric for the parameters used in the grid search. the metric is the sum of the 3 best scores divided by the sum of the next 3 best scores

parser = argparse.ArgumentParser(description = 'Run the PepGM algorithm from command line')
parser.add_argument('--resultsfolder', required = True, help = 'folder with the results from the beliefpropagation in cvs format')
parser.add_argument('--out', required=True, help ='output png file')

args = parser.parse_args()

#resultsfolder = '/home/tholstei/repos/PepGM_all/PepGM/results/PXD002936_avian_bronchitis'
#paramcheck = '/home/tholstei/repos/PepGM_all/PepGM/results/PXD002936_avian_bronchitis/paramcheck.png'

Params = []

Metrics = []

for folders in os.listdir(args.resultsfolder):
    if os.path.isdir(args.resultsfolder+ '/' + folders):
        for file in os.listdir(args.resultsfolder+ '/' + folders):

            if file.endswith('.csv'):
        
                Results = pd.read_csv(args.resultsfolder+ '/' + folders +'/' + file,names = ['ID','score','type'])
                TaxIDS = Results.loc[Results['type']=='taxon']
                TaxIDS.loc[:,'score'] = pd.to_numeric(TaxIDS['score'],downcast = 'float')
                TaxIDS = TaxIDS.sort_values('score', ascending = False)
        
                #compute the metric

                #compare the posterior probbilities
                FirstSum = np.sum(TaxIDS[:3]['score'])
                NextSum = np.sum(TaxIDS[3:6]['score'])
                SumProportion = FirstSum/NextSum
                print(SumProportion,'sumprop')

                #compute the pairwise taxonomic distance of the first 4 Taxa
                FourTaxa = np.array(TaxIDS[:4]['ID'])
                tree = ncbi.get_topology(FourTaxa)
                DistanceSum = (tree.get_distance(FourTaxa[0],FourTaxa[1]))**9+tree.get_distance(FourTaxa[1],FourTaxa[2])+tree.get_distance(FourTaxa[2],FourTaxa[3])
                print(DistanceSum)
                Matching = SumProportion/DistanceSum
                print(Matching,'match')

                Metrics.append(Matching)

                #get the corresponding parameters
                params = re.findall('\d+\.\d+',file)
               
                
        
                Params.append([params[0],params[1],params[2]])

zipLists = zip(Metrics,Params)
SortedPairs = sorted(zipLists,reverse= True)

tuples = zip(*SortedPairs)
Metrics,Params = [list(tuple) for tuple in tuples]

figure1 = plt.figure(figsize=(30,15), tight_layout=True)
test = Metrics[0:18]
test1 = list(range(20))
print(test)

plt.barh(list(range(20)),Metrics[0:20],color = 'mediumvioletred')
plt.xticks(fontsize =35)
plt.yticks(list(range(20)),Params[0:20],fontsize=35)
#plt.set_ticklabels(fontsize =35)

plt.savefig(args.out)




    




