#from https://gitlab.com/rki_bioinformatics/gnomo/-/blob/master/scripts/unipept-get-peptinfo.py


import requests
import json
import re
import re
import os.path
from ete3 import NCBITaxa

ncbi = NCBITaxa()

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
                    pep_score[peptide] = pep
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


def generatePostRequest(peptides,TargetTaxa):
    '''
    Generates POST request (json) from petide and target taxon
    :param peptides: list of peptides to query in Unipept
    :param TargetTaxa: list of one or more taxa to include in the Unipept query
    '''
    
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
    
    url = "https://sherlock.ugent.be/mpa/pept2filtered.json"
    request = requests.post(url,json.dumps(request_json),headers={'content-type':'application/json'})

    with open(out_file, 'w') as f_out:
        print(request.text,file=f_out)








if __name__=='__main__':
    pout_file = '/home/tholstei/repos/PepGM_all/PepGM/results/Xtandem_rescore_test/PXD018594_Sars_CoV_2/MS2Rescore/rescore.pin_searchengine_ms2pip_rt_features.pout'
    pep_score_psm = Poutparser(pout_file,0.05,'')


    UnipeptPeptides = dict()
    for peptide in pep_score_psm.keys():
        FullyTrypticPeptides = PepListNoMissedCleavages(peptide)
        for pep in FullyTrypticPeptides:
            UnipeptPeptides[pep] = pep_score_psm[peptide]

    #out = getInfoFromUnipept(UnipeptPeptides,'test_unipept.csv' )
    
    request = generatePostRequest(list(UnipeptPeptides.keys()),[11118])


    out = PostInfoFromUnipept(request,'test_unipept.json')