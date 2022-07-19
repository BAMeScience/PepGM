import numpy as np
from ete3 import NCBITaxa
from csv import reader
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('UnipeptResponseFile', type = str, required = True, help = 'path to Unipept response .csv file')
#parser.add_argument('taxonomic Resolution') include this file later
parser.add_argument('NumberOfTaxa', type = int, required = True, help = 'number of taxa to include in the output' )

def GetAllTaxa(InputFile,TaxResolution):
    '''
    Takes .csv file with Unipept returns for a list of peptides,
    queries taxa at the TaxResolution level for all peptides, taxon weight = #psms
    :input InputFile: str, path to Unipept responses in .csv file
    :input TaxResolution: str,Taxonomic resolution of the returned taxa per peptide
    :return PepToTaxaDict: dict, nested dictionary {peptide:{taxa,weight,score}}
    '''

    ncbi = NCBITaxa()
    PepToTaxaDict = {}
    with open(InputFile,'r') as File:
        csv_reader = reader(File)
        header = next(csv_reader)
        for row in csv_reader:
            check = row
            if row[1] == '1':
                continue
            if row[5] == 'None':
                if row[4] == 'None':
                    if row[3] == 'None':
                        if row[2] == 'None':
                            continue
                            MappedTaxa = ncbi.get_descendant_taxa(row[1],collapse_subspecies=True)
                            PepToTaxaDict.update({row[0]:{'Taxa': MappedTaxa , 'weight':float(row[8])/len(MappedTaxa),'score':float(row[7])}})
                        else:
                            continue
                            MappedTaxa = ncbi.get_descendant_taxa(row[2],collapse_subspecies=True)
                            print(MappedTaxa)
                            PepToTaxaDict.update({row[0]:{'Taxa': MappedTaxa, 'weight':float(row[8])/len(MappedTaxa),'score':float(row[7])}})
                    else:
                        MappedTaxa = [tax for tax in ncbi.get_descendant_taxa(row[3],collapse_subspecies=True)]
                        PepToTaxaDict.update({row[0]:{'Taxa': MappedTaxa, 'weight':float(row[8])/len(MappedTaxa),'score':float(row[7])}})
                else:
                    MappedTaxa = [tax for tax in ncbi.get_descendant_taxa(row[4],collapse_subspecies=True)]
                    PepToTaxaDict.update({row[0]:{'Taxa': MappedTaxa, 'weight':float(row[8])/len(MappedTaxa),'score':float(row[7])}})
            else:
                PepToTaxaDict.update({row[0]:{'Taxa':[row[5]], 'weight':float(row[8]),'score':float(row[7])}})

    return PepToTaxaDict


        
def WeightTaxa(Pep2TaxDict,TaxNum):
    '''
    Takes dictionary {peptide:{taxa,weight,score}} as input and rearanges it to {taxa:{weight,{peptides:{peptide:score}}} 
    :input Pep2TaxDict: dict, mapping peptides to taxa, weight and score
    :input TaxNum: int, # of taxa to be included in the returned list
    :output WeightSortedTaxa: list of tuples, (taxon,attribute dict) '''

    TaxaDict = {}
    for peptide in Pep2TaxDict.keys():
        Taxa = Pep2TaxDict[peptide]['Taxa']
        for tax in Taxa:
            try:
                TaxaDict[tax]['weight']= TaxaDict[tax]['weight']+Pep2TaxDict[peptide]['weight']
                TaxaDict[tax]['peptides'][peptide]= Pep2TaxDict[peptide]['score']
            except :
                TaxaDict.update({tax : {'weight': Pep2TaxDict[peptide]['weight'], 'peptides': {peptide: Pep2TaxDict[peptide]['score']}}})

    WeightSortedTaxa = sorted(TaxaDict.items(), key=lambda x:x[1]['weight'],reverse = True)
    return WeightSortedTaxa[0:TaxNum]







if __name__=='__main__':
    InputFile = '/home/tholstei/repos/PepGM_all/PepGM/test_unipept.csv'
    Dict1 = GetAllTaxa(InputFile,0.5)
    Dict2 = WeightTaxa(Dict1,100)
    print(Dict2)
    #ncbi = NCBITaxa()
    #test = ncbi.get_descendant_taxa(7711,collapse_subspecies = True)
    #print(len([x for x in test]))
   
