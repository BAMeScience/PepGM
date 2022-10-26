import numpy as np
from ete3 import NCBITaxa
from csv import reader
import argparse
import pandas as pd
import json
import UnipeptGetTaxonomyfromPout as Unipept

parser = argparse.ArgumentParser()

parser.add_argument('--UnipeptResponseFile', type = str, required = True, help = 'path to Unipept response .json file')
parser.add_argument('--TaxonomyQuery', required = True, help = 'taxa to query in Unipept. If querying all taxa, put [1]') 
parser.add_argument('--NumberOfTaxa', type = int, required = True, help = 'number of taxa to include in the output' )
parser.add_argument('--MinScore', type = float, required = True, help = 'min peptide score for the peptide to be included in the search')
parser.add_argument('--PoutFile', type = str, required = True, help = 'path to percolator(ms2rescore) Pout file')
parser.add_argument('--out', type = str, required = True, help = 'path to csv out file')


args = parser.parse_args()

def WeightAllTaxaFromJson(JsonPath,PeptScoreDict,MaxTax):
    '''
    Takes Unipept response Json and returns a dataframe with peptide sequence, corresponding taxon, corresponding FAs, number of psms and peptide score filtered by PSM-weighted taxa
    :param JsonPath: str, path toUnipept response JSON file
    :param PeptScoreDict: dict, dictionary {peptide:[peptide_score,#psms]}
    :param MaxTax: int, maximum number of taxa to include
    :return: pandas dataframe including all petide-taxa links where the taxon weight was above the median taxon weight or where max. MaxTax taxa were present
    '''

    with open(JsonPath,'r') as file:
        UnipeptDict = json.load(file)
    
    UnipeptFrame = pd.json_normalize(UnipeptDict,record_path = ['peptides'])
    UnipeptFrame['psms_score']= UnipeptFrame['sequence'].map(PeptScoreDict)
    print(UnipeptFrame['psms_score'])
    UnipeptFrame = pd.concat([UnipeptFrame.drop('psms_score',axis=1),pd.json_normalize(UnipeptFrame['psms_score'])],axis = 1)
    UnipeptFrame['weight']= UnipeptFrame['psms'].div([len(element) for element in UnipeptFrame['taxa']])
    UnipeptFrame = UnipeptFrame.explode('taxa',ignore_index = True)

    UnipeptFrameTaxaWeights = UnipeptFrame.groupby('taxa')['weight'].sum().reset_index()
    UnipeptFrameTaxaWeights = UnipeptFrameTaxaWeights.sort_values(by=['weight'], ascending=False)
    threshold = UnipeptFrameTaxaWeights.weight.median()
    TopTaxa = UnipeptFrameTaxaWeights.loc[UnipeptFrameTaxaWeights["weight"] >= threshold]

    if len(TopTaxa.taxa)<50:
        return UnipeptFrame[UnipeptFrame['taxa'].isin(TopTaxa.taxa)]
    else:
        return UnipeptFrame[UnipeptFrame['taxa'].isin(TopTaxa.taxa[0:MaxTax])]
   

#get peptides, spsms and score from pout file into dicitonary shape
pep_score_psm = Unipept.Poutparser(args.PoutFile,args.MinScore,'')
UnipeptPeptides = dict()
for peptide in pep_score_psm.keys():
    FullyTrypticPeptides = Unipept.PepListNoMissedCleavages(peptide)
    for pep in FullyTrypticPeptides:
        UnipeptPeptides[pep] ={'score':pep_score_psm[peptide][0], 'psms':pep_score_psm[peptide][1]} 

#get and save Info from Unipept
request = Unipept.generatePostRequest(list(UnipeptPeptides.keys()),[int(item) for item in args.TaxonomyQuery.split(',')])    
save = Unipept.PostInfoFromUnipept(request,args.UnipeptResponseFile)


#format and return dataframe with weighted taxa
DF = WeightAllTaxaFromJson(args.UnipeptResponseFile,UnipeptPeptides,args.NumberOfTaxa)
DF.to_csv(args.out)







if __name__=='__main__':

    pout_file = '/home/tholstei/repos/PepGM_all/PepGM/results/Xtandem_rescore_test/PXD018594_Sars_CoV_2/MS2Rescore/rescored_searchengine_ms2pip_rt_features.pout'
    pep_score_psm = Unipept.Poutparser(pout_file,0.05,'')
    
    print(pd.__version__)
    UnipeptPeptides = dict()
    for peptide in pep_score_psm.keys():
        FullyTrypticPeptides = Unipept.PepListNoMissedCleavages(peptide)
        for pep in FullyTrypticPeptides:
            UnipeptPeptides[pep] ={'score':pep_score_psm[peptide][0], 'psms':pep_score_psm[peptide][1]} 
    
    print(UnipeptPeptides)
    request = Unipept.generatePostRequest(list(UnipeptPeptides.keys()),[11118])
    
    
    out = Unipept.PostInfoFromUnipept(request,'test_unipept.json')
    InputFile = '/home/tholstei/repos/PepGM_all/PepGM/test_unipept.json'
    DF = WeightAllTaxaFromJson(InputFile,UnipeptPeptides,30)
    print(DF)


