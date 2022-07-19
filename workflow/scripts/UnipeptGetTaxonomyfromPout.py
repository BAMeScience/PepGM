#from https://gitlab.com/rki_bioinformatics/gnomo/-/blob/master/scripts/unipept-get-peptinfo.py

import urllib.request
import json
import argparse
import re
import re
import os.path

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




def generateGetRequest(peptides):
    """
    Given a list of Unipept peptides, it generates a list of Unipept Get Requests
    :param peptides: Unipept peptides
    :return: list of get requests for Unipept
    """
    request_list = list()
    char_count = 74  # number of chars minimum per request
    request = "http://api.unipept.ugent.be/api/v1/peptinfo.json?"
    for peptide in peptides:
        if (len(peptide) + 9 + char_count) < 2048:
            request += "input[]={}&".format(peptide)
            char_count += 9 + len(peptide)
        else:
            request += "equate_il=true&extra=true&names=true"
            request_list.append(request)
            request = "http://api.unipept.ugent.be/api/v1/peptinfo.json?" + "input[]={}&".format(peptide)
            char_count = 74 + 9 + len(peptide)

    if len(request) != 49:  # first part of request
        request += "equate_il=true&extra=true"
        request_list.append(request)


    return request_list



def getInfoFromUnipept(request_dict, result_file):
    """
    Send all requests, get for each peptide the phylum, family, genus and collection of EC-numbers
    :param request_list: list of Get Requests
    :param result_file: csv file with Unipept info (phylum, family, genus and collection of EC-numbers)
    :return: None
    """

    unipept_list = list()
    requests = []
    
    requests = generateGetRequest(request_dict.keys())
    

    for element in requests:
        unipept_json = urllib.request.urlopen(element).read()
        unipept_list.append(json.loads(unipept_json.decode('utf-8')))

    
    tmp = open(result_file, "w")
    tmp.close()


    with open(result_file, "a") as f_out:
        print('peptide,','lca_id,','phylum_id,','family_id,','genus_id,','species_id,','ec,','score,','number of psms',file=f_out)
        for response in unipept_list:
            for element in response:
                lca_id = element['taxon_id']
                phylum_id = element["phylum_id"]
                family_id = element["family_id"]
                genus_id = element["genus_id"]
                species_id = element["species_id"]
                try:
                    ec_list = list()
                    for i in element['ec']:
                        ec_list.append(i['ec_number'].strip())
                    ec = ";".join(ec_list)
                except:
                    ec = ""
                print("{},{},{},{},{},{},{},{},{}".format(element["peptide"].strip(), lca_id, phylum_id, family_id, genus_id, species_id, ec, 
                        1-float(request_dict[element["peptide"].strip()][0]),request_dict[element["peptide"].strip()][1]), file=f_out)


if __name__=='__main__':
    pout_file = '/home/tholstei/repos/PepGM_all/PepGM/results/Xtandem_rescore_test/PXD018594_Sars_CoV_2/MS2Rescore/rescore.pin_searchengine_ms2pip_rt_features.pout'
    pep_score_psm = Poutparser(pout_file,0.05,'')


    UnipeptPeptides = dict()
    for peptide in pep_score_psm.keys():
        FullyTrypticPeptides = PepListNoMissedCleavages(peptide)
        for pep in FullyTrypticPeptides:
            UnipeptPeptides[pep] = pep_score_psm[peptide]

    out = getInfoFromUnipept(UnipeptPeptides,'test_unipept.csv' )