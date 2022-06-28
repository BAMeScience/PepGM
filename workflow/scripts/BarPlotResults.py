import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
from ete3 import NCBITaxa
import argparse

matplotlib.use('Agg')

parser = argparse.ArgumentParser(description = 'generate BarPlot of PepGM results')

parser.add_argument('--ResultsFile', type = str, help = 'path(s) to your PepGM results CSV')
parser.add_argument('--NumberofResults', type = int, default = 12, help = 'how many taxa you want to show up on the results plot')
parser.add_argument('--out',type =str, help = 'path(s) to your results file')

args = parser.parse_args()

'''
Script that takes PepGM .csv output, translates taxIDS to scientific names, and barplots the *number of results* highest scoring taxa
'''


ncbi = NCBITaxa()

#read csv using pandas  
IDs = pd.read_csv(args.ResultsFile, names = ['ID','score','type'])
TaxIDS = IDs.loc[IDs['type']=='taxon']
TaxIDS.loc[:,'score'] = pd.to_numeric(TaxIDS['score'],downcast = 'float')
TaxIDS = TaxIDS.sort_values('score')
TaxaCheck = TaxIDS.ID.tolist()

#translate taxids to scientific names
TaxaNameDict = ncbi.get_taxid_translator(TaxIDS['ID'])
TaxaNames = [TaxaNameDict[int(tax)] for tax in TaxaCheck]
Scores = TaxIDS['score']



#make the barplot
fig, ax = plt.subplots()
fig.set_size_inches(30,15)
ax.barh(range(len(TaxaNames[-args.NumberofResults:])),Scores[-args.NumberofResults:], color='royalblue')

ax.set_yticks(range(len(TaxaNames[-args.NumberofResults:])))
ax.set_yticklabels(TaxaNames[-args.NumberofResults:], fontsize = 35)
plt.xlim((0,1))
plt.xlabel('Posterior probability',fontsize=35)
ax.xaxis.set_ticks(np.arange(0,1.2,0.2))
ax.xaxis.set_ticklabels([0,0.2,0.4,0.6,0.8,1.0], fontsize =35)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

fig.tight_layout()
        
plt.savefig(args.ResultsFile.replace('.csv','.png'))
plt.close()
    
 

