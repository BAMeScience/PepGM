import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
from ete3 import NCBITaxa

ncbi = NCBITaxa()
ResultsFile = '/home/tholstei/repos/PepGM_all/PepGM/results/Cowpox_cowpox_removed/PXD014913_cowpox_BR/Prior0.5/refseqViralwocpxv_PepGM_Results_a0.01_b0.7_p0.5.png'
NumberofResults = 15
out = '/home/tholstei/repos/PepGM_all/PepGM/results/Cowpox_cowpox_removed/PXD014913_cowpox_BR/PrettyPlot.png'



TaxIDS = pd.read_csv(ResultsFile)
TaxIDS.loc[:,'score'] = pd.to_numeric(TaxIDS['score'],downcast = 'float')
TaxIDS = TaxIDS.head(15)
print(TaxIDS)
TaxaCheck = TaxIDS.ID.tolist()

#translate taxids to scientific names
TaxaNameDict = ncbi.get_taxid_translator(TaxIDS['ID'])
TaxaNames = [TaxaNameDict[int(tax)] for tax in TaxaCheck]
Scores = TaxIDS['score']



#make the barplot
fig, ax = plt.subplots()
fig.set_size_inches(30,15)
bars = ax.barh(range(len(TaxaNames)),Scores[::-1], color='royalblue')

ax.set_yticks(range(len(TaxaNames)))
ax.set_yticklabels(TaxaNames[::-1], fontsize = 35)
plt.xlim((0,1))
plt.xlabel('Posterior probability',fontsize=35)
ax.xaxis.set_ticks(np.arange(0,1.2,0.2))
ax.xaxis.set_ticklabels([0,0.2,0.4,0.6,0.8,1.0], fontsize =35)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

print(bars)

bar_color = bars[0].get_facecolor()

for bar in bars:
  ax.text(
      bar.get_width() + 0.05,
      bar.get_y() + bar.get_height() / 5,
      round(bar.get_width(), 2), fontsize = 35,
      horizontalalignment='center',
      color=bar_color,
      weight='bold'
  )

fig.tight_layout()
        
plt.savefig(out)
plt.close()