from LoadSimplePSMResults import loadSimplePepScore


import requests
import json
import re
import re
from ete3 import NCBITaxa
import argparse

ncbi = NCBITaxa()

parser = argparse.ArgumentParser()

parser.add_argument('--UnipeptResponseFile', type = str, required = True, help = 'path to Unipept response .json file')
parser.add_argument('--TaxonomyQuery', required = True, help = 'taxa to query in Unipept. If querying all taxa, put [1]') 
parser.add_argument('--NumberOfTaxa', type = int, required = True, help = 'number of taxa to include in the output' )
parser.add_argument('--PSMResultsFile', type = str, required = True, help = 'path to percolator(ms2rescore) Pout file')
parser.add_argument('--pep_out', type = str, required = True, help = 'path to csv out file')


args = parser.parse_args()

def PepListNoMissedCleavages(Pepnames):
    """
    Takes a peptide and cleaves it into Unipept format (0 missed cleavages,
    cleaves after K or R except followed by P)
    :param peptides_in: list of peptides
    :return: list of peptides in Unipept format
    """
    
    peptides = list()
    trypsin = lambda peptide: re.sub(r'(?<=[RK])(?=[^P])', '\n', peptide, re.DOTALL).split()
    peptides += trypsin(peptide)

    return peptides



def generatePostRequestChunks(peptides,TargetTaxa,chunksize=20):
    '''
    Generates POST requests (json) queryong a chunk of peptides from petide list and target taxon
    :param peptides: list of peptides to query in Unipept
    :param TargetTaxa: list of one or more taxa to include in the Unipept query
    :param chunksize: number of peptides to be requested from Unipept
    '''
    
    AllTargetTaxa = []
    for Taxon in TargetTaxa:
        AllTargetTaxa.append(Taxon)
        AllTargetTaxa.extend(ncbi.get_descendant_taxa(Taxon))
    
    
    Listofpeptides = [peptides[i:i + chunksize] for i in range(0, len(peptides), chunksize)]
    Listofrequests = [{"peptides":chunk, "taxa":AllTargetTaxa} for chunk in Listofpeptides]

    return Listofrequests


def PostInfoFromUnipeptChunks(request_json, out_file):
    """
    Send all requests, get for each peptide the phylum, family, genus and collection of EC-numbers
    :param request_list: list of Get Requests
    :param result_file: csv file with Unipept info (phylum, family, genus and collection of EC-numbers)
    :return: None
    """
    
    url = "http://api.unipept.ugent.be/mpa/pept2filtered.json"
    print('now querying Unipept')

    for chunk in request_json:
        request = requests.post(url,json.dumps(chunk),headers={'content-type':'application/json'}, timeout = None)    
        with open(out_file, 'a') as f_out:
            print(request.text,file=f_out)


Pepnames,Pepscores = loadSimplePepScore(args.PSMResultsFile)
PepScoreDict = dict(zip(Pepnames,Pepscores))

PSMnumber = dict()
for i in Pepnames:
    PSMnumber[i] = PSMnumber.get(i, 0) + 1

UnipeptPeptides = dict()
for peptide in PepScoreDict.keys():
    FullyTrypticPeptides = PepListNoMissedCleavages(peptide)
    for pep in FullyTrypticPeptides:
        UnipeptPeptides[pep] = {'score': PepScoreDict[peptide], 'psms': PSMnumber[peptide]} 
with open(args.pep_out, 'a') as f_out:
    f_out.write( json.dumps(UnipeptPeptides))


#get and save Info from Unipept if the response file doesn't exist yet
request = generatePostRequestChunks(list(PepScoreDict.keys()),[int(item) for item in args.TaxonomyQuery.split(',')])    
save = PostInfoFromUnipeptChunks(request,args.UnipeptResponseFile)

#if __name__=='__main__':
#    PSMResultsFile = '/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD002936_avian_bronchitis/chicken_refseq_Default_PSM_Report.txt'
   
#    Pepnames,Pepscores = loadSimplePepScore(PSMResultsFile)

#    PSMnumber = dict()
#    for i in Pepnames:
#        PSMnumber[i] = PSMnumber.get(i, 0) + 1
    
#    PepScoreDict = dict(zip(Pepnames,Pepscores))

#    UnipeptPeptides = dict()
#    for peptide in PepScoreDict.keys():
#        FullyTrypticPeptides = PepListNoMissedCleavages(peptide)
#        for pep in FullyTrypticPeptides:
#            UnipeptPeptides[pep] = {'score': PepScoreDict[peptide], 'psms': PSMnumber[peptide]}
#    with open('unipept_peptides_test.json', 'a') as f_out:
#        f_out.write( json.dumps(UnipeptPeptides))

#    request = generatePostRequestChunks(list(UnipeptPeptides.keys()),[11118])
#    out = PostInfoFromUnipeptChunks(request,'test_unipept.json')