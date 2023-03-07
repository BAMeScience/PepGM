#from https://gitlab.com/rki_bioinformatics/gnomo/-/blob/master/scripts/unipept-get-peptinfo.py


import requests
import json
import re
import re
import os.path
from ete3 import NCBITaxa
import argparse

ncbi = NCBITaxa()

parser = argparse.ArgumentParser()

parser.add_argument('--UnipeptResponseFile', type = str, required = True, help = 'path to Unipept response .json file')
parser.add_argument('--TaxonomyQuery', required = True, help = 'taxa to query in Unipept. If querying all taxa, put [1]') 
parser.add_argument('--NumberOfTaxa', type = int, required = True, help = 'number of taxa to include in the output' )
parser.add_argument('--FDR', type = float, required = True, help = 'min peptide score for the peptide to be included in the search')
parser.add_argument('--PoutFile', type = str, required = True, help = 'path to percolator(ms2rescore) Pout file')
parser.add_argument('--pep_out', type = str, required = True, help = 'path to csv out file')


args = parser.parse_args()

#code adapted from pout to prot
def Poutparser(pout_file, fdr_threshold, decoy_flag):
    '''
    Parses the ms2rescore pout file for peptides, psm numbers and peptide scores
    :param pout_file: str, path to pout file
    :param fdr_threshold: float, fdr threshold below which psms are kept
    :param decoy_flag: str, can be emtpy string, decoy flag in pout file
    :return: dict, peptides:[score,#psms]
    '''
    
    pep_score = dict()
    pep_psm = dict()
    pep_score_psm = dict()

    assert os.path.exists(pout_file), "input file or folder does not exist"

    with open(pout_file, "r") as f:
        next(f)  # skip header
        for line in f:
            # skip empty lines
            if line.rstrip() == "":
                continue
            splitted_line = line.rstrip().split("\t", maxsplit=5)
            assert len(splitted_line) >= 6, "Input file is wrongly formatted. Make sure that the input is a valid .pout file."
            psm_id, _, q, pep, peptide,_ = splitted_line
            if float(q) < fdr_threshold:
                peptide = re.sub("\[.*?\]", "", peptide)
                peptide = peptide.split(".")[1]
                # update pep_psm
                if peptide not in pep_psm.keys():
                    pep_psm[peptide] = set()
                    pep_psm[peptide].add(psm_id)
                else:
                    pep_psm[peptide].add(psm_id)
                # update pep_score
                if peptide not in pep_score.keys():
                    if float(pep) <0.001:
                        pep_score[peptide] = '0.001'
                    else:
                        pep_score[peptide] = pep          #adjustement necessary to not have 0 and 1 fuck up probability calculations
                else:
                    if float(pep) <0.001:
                        pep_score[peptide] = '0.001'
                    else:
                        pep_score[peptide] = min(pep,pep_score[peptide])
                pep_score_psm[peptide] = [pep_score[peptide],len(pep_psm[peptide])]
            
    return pep_score_psm

def PepListNoMissedCleavages(peptide):
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
    print('querying taxa ', TargetTaxa)
    AllTargetTaxa = []
    for Taxon in TargetTaxa:
        AllTargetTaxa.append(Taxon)
        AllTargetTaxa.extend(ncbi.get_descendant_taxa(Taxon, collapse_subspecies=True))
    
    
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

def generatePostRequest(peptides,TargetTaxa):
    '''
    Generates POST request (json) from petide and target taxon
    :param peptides: list of peptides to query in Unipept
    :param TargetTaxa: list of one or more taxa to include in the Unipept query
    '''

    print('number of peptides to be queried: ', len(peptides))
    
    AllTargetTaxa = []
    for Taxon in TargetTaxa:
        AllTargetTaxa.append(Taxon)
        AllTargetTaxa.extend(ncbi.get_descendant_taxa(Taxon))
    
    request = {"peptides":peptides, "taxa":AllTargetTaxa}

    return request



def PostInfoFromUnipept(request_json, out_file):
    """
    Send all requests, get for each peptide the phylum, family, genus and collection of EC-numbers
    :param request_list: list of Get Requests
    :param result_file: csv file with Unipept info (phylum, family, genus and collection of EC-numbers)
    :return: None
    """
    
    url = "http://api.unipept.ugent.be/mpa/pept2filtered.json"
    print('now querying Unipept')
    request = requests.post(url,json.dumps(request_json),headers={'content-type':'application/json'},timeout=None )    
    
    with open(out_file, 'w') as f_out:
        print(request.text,file=f_out)



pep_score_psm = Poutparser(args.PoutFile,args.FDR,'')
UnipeptPeptides = dict()
for peptide in pep_score_psm.keys():
    FullyTrypticPeptides = PepListNoMissedCleavages(peptide)
    for pep in FullyTrypticPeptides:
        UnipeptPeptides[pep] ={'score':pep_score_psm[peptide][0], 'psms':pep_score_psm[peptide][1]} 
with open(args.pep_out, 'a') as f_out:
    f_out.write( json.dumps(UnipeptPeptides))

#get and save Info from Unipept if the response file doesn't exist yet
request = generatePostRequestChunks(list(UnipeptPeptides.keys()),[int(item) for item in args.TaxonomyQuery.split(',')])    
save = PostInfoFromUnipeptChunks(request,args.UnipeptResponseFile)


#if __name__=='__main__':
#    pout_file = '/home/tholstei/repos/PepGM_all/PepGM/results/Xtandem_rescore_test/PXD018594_Sars_CoV_2/MS2Rescore/rescored_searchengine_ms2pip_rt_features.pout'
#    pep_score_psm = Poutparser(pout_file,0.05,'')


#    UnipeptPeptides = dict()
#    for peptide in pep_score_psm.keys():
#        FullyTrypticPeptides = PepListNoMissedCleavages(peptide)
#        for pep in FullyTrypticPeptides:
#            UnipeptPeptides[pep] = pep_score_psm[peptide]

    #out = getInfoFromUnipept(UnipeptPeptides,'test_unipept.csv' )
#    print(list(UnipeptPeptides.keys()))
#    request = generatePostRequest(list(UnipeptPeptides.keys()),[11118])


#    out = PostInfoFromUnipept(request,'test_unipept.json')