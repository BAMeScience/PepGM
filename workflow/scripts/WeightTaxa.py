import numpy as np
from ete3 import NCBITaxa
from csv import reader
import argparse
import pandas as pd
import json
import os.path

parser = argparse.ArgumentParser()

parser.add_argument('--UnipeptResponseFile', type = str, required = True, help = 'path to Unipept response .json file')
parser.add_argument('--NumberOfTaxa', type = int, required = True, help = 'number of taxa to include in the output' )
parser.add_argument('--out', type = str, required = True, help = 'path to csv out file')
parser.add_argument('--UnipeptPeptides', type = str, required = True, help = 'path to Unipept response .json file')


args = parser.parse_args()

def WeightAllTaxaFromJson(JsonPath,PeptScoreDict,MaxTax,chunks = True):
    '''
    Takes Unipept response Json and returns a dataframe with peptide sequence, corresponding taxon, corresponding FAs, number of psms and peptide score filtered by PSM-weighted taxa
    :param JsonPath: str, path toUnipept response JSON file
    :param PeptScoreDict: dict, dictionary {peptide:[peptide_score,#psms]}
    :param MaxTax: int, maximum number of taxa to include
    :return: pandas dataframe including all petide-taxa links where the taxon weight was above the median taxon weight or where max. MaxTax taxa were present
    '''

    with open(PeptScoreDict,'r') as file:
        PeptScoreDictload = json.load(file)

    if chunks:
        with open(JsonPath,'r') as file:
            UnipeptDict = {"peptides":[]}
            for line in file:
                try:
                    print()
                    UnipeptDict["peptides"].extend(json.loads(line)["peptides"])
                except KeyError:
                    UnipeptDict["peptides"]=[json.loads(line)["peptides"]]

    
    else:
        with open(JsonPath,'r') as file:
            UnipeptDict = json.load(file)


    UnipeptFrame = pd.json_normalize(UnipeptDict,record_path = ['peptides'])
    UnipeptFrame['psms_score']= UnipeptFrame['sequence'].map(PeptScoreDictload)
    UnipeptFrame = pd.concat([UnipeptFrame.drop('psms_score',axis=1),pd.json_normalize(UnipeptFrame['psms_score'])],axis = 1)
    UnipeptFrame['weight']= UnipeptFrame['psms'].div([len(element) for element in UnipeptFrame['taxa']])
    UnipeptFrame = UnipeptFrame.explode('taxa',ignore_index = True)

    UnipeptFrameTaxaWeights = UnipeptFrame.groupby('taxa')['weight'].sum().reset_index()
    UnipeptFrameTaxaWeights = UnipeptFrameTaxaWeights.sort_values(by=['weight'], ascending=False)
    UnipeptFrameTaxaWeights.to_csv('/home/tholstei/repos/PepGM_all/PepGM/results/Xtandem_rescore_decoy_test_full_Uniprot/CAMPI_SIHUMIx/fullweightframe.csv')
    threshold = UnipeptFrameTaxaWeights.weight.median()
    TopTaxa = UnipeptFrameTaxaWeights.loc[UnipeptFrameTaxaWeights["weight"] >= threshold]

    if len(TopTaxa.taxa)<50:
        return UnipeptFrame[UnipeptFrame['taxa'].isin(TopTaxa.taxa)]
    else:
        return UnipeptFrame[UnipeptFrame['taxa'].isin(TopTaxa.taxa[0:MaxTax])]
   




#get peptides, spsms and score from pout file into dicitonary shape
#careful: the ms2rescore 'score' is the posterior error probability!!



#format and return dataframe with weighted taxa
DF = WeightAllTaxaFromJson(args.UnipeptResponseFile,args.UnipeptPeptides,args.NumberOfTaxa)

DF.to_csv(args.out)







if __name__=='__main__':

    InputFile = '/home/tholstei/repos/PepGM_all/PepGM/results/SIHUMIx_no_subspecies/CAMPI_SIHUMIx/UnipeptResponse.json'
    DF = WeightAllTaxaFromJson(InputFile,'/home/tholstei/repos/PepGM_all/PepGM/results/SIHUMIx_no_subspecies/CAMPI_SIHUMIx/UnipeptPeptides.json',300)
    DF.to_csv('/home/tholstei/repos/PepGM_all/PepGM/results/SIHUMIx_no_subspecies/CAMPI_SIHUMIx/differentGraphs300.csv')


