import json
import matplotlib.pyplot as plt
from ete3 import NCBITaxa
import argparse
import numpy as np


'''
Script that plots the number of proteins per taxon in the dict used to build the factor graph

'''

parser = argparse.ArgumentParser(description = 'count the number of proteins per taxon in json dictionary')
parser.add_argument('--ResultsFile', type = str, help = 'path(s) to your PepGM results CSV')
parser.add_argument('--NumberofResults', type = int, default = 12, help = 'how many taxa you want to show up on the results plot')
parser.add_argument('--out',type =str, help = 'path(s) to your results file')

args = parser.parse_args()

ncbi = NCBITaxa()

with open (args.ResultsFile) as file:
    TaxonPeptideDict = json.load(file)

TaxIDS = list(TaxonPeptideDict.keys())
Proteins = TaxonPeptideDict.values()
ProteinNumber = [len(protlist) for protlist in Proteins]
ProteinNumber_sorted = sorted(ProteinNumber)
print(ProteinNumber_sorted)

SortingIndex = np.argsort(ProteinNumber)
TaxIDS[:] = [TaxIDS[i] for i in SortingIndex]



TaxaNameDict = ncbi.get_taxid_translator(TaxonPeptideDict.keys())
TaxaNames = [TaxaNameDict[int(tax)] for tax in TaxIDS]


fig, ax = plt.subplots()
fig.set_size_inches(30,15)
ax.barh(range(len(TaxaNames[-args.NumberofResults:])),ProteinNumber_sorted[-args.NumberofResults:], color='red')

ax.set_yticks(range(len(TaxaNames[-args.NumberofResults:])))
ax.set_yticklabels(TaxaNames[-args.NumberofResults:], fontsize = 20)
plt.xlabel('Posterior probability',fontsize=35)
#ax.set_xticklabels(ax.get_xticklabels, fontsize =35)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

fig.tight_layout()
        
plt.savefig(args.out)
plt.close()
