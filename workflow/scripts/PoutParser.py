import re
import os.path

#code adapted from pout to prot
def parser(pout_file, fdr_threshold, decoy_flag):
    ### Initiate dictionaries ###
    # 1) PSM -> experiment dict(PSM : exp (string of file name), ...)
    pep_score = dict()

    # 2) Peptide -> Set<PSM>
    pep_psm = dict()

    ### Read pout file ###
    assert os.path.exists(pout_file), "input file or folder does not exist"

    with open(pout_file, "r") as f:
        next(f)  # skip header
        for line in f:
            # skip empty lines
            if line.rstrip() == "":
                continue
            splitted_line = line.rstrip().split("\t", maxsplit=5)
            assert len(splitted_line) >= 6, "Input file is wrongly formatted. Make sure that the input is a valid .pout file."
            psm_id, _, q, pep, peptide, proteins = splitted_line
            if float(q) < fdr_threshold:
                # retrieving all information
                peptide = re.sub("\[.*?\]", "", peptide)
                peptide = peptide.split(".")[1]
                proteins = proteins.split("\t")
                for i, protein in enumerate(proteins):
                    # SwissProt proteins: # >sp|<accession>|<description>
                    try:
                        protein = protein.split("|", maxsplit=2)[1]
                        if protein != "" and protein != "sp" and protein != "tr":  # could be the case in typo,
                            # e.g. >generic||<accession>|<description>; sometimes empty with sp or tr
                            proteins[i] = protein
                    except IndexError:
                        pass  # leave protein name just as complete header
                # Normally, only decoy proteins are returned IF there's other proteins linked to a PSM as well
                assert decoy_flag == "" or not all(decoy_flag in protein for protein in proteins), \
                    "{} is only contained in decoy proteins".format(psm_id)
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

    return pep_score, pep_psm


if __name__ == "__main__":
    pep_score,pep_psm = parser('/home/tholstei/repos/PepGM_all/PepGM/results/Xtandem_rescore_test/PXD018594_Sars_CoV_2/MS2Rescore/rescore.pin_searchengine_ms2pip_rt_features.pout',0.05,'')
    pep_psm = dict([(peptide,len(pep_psm[peptide])) for peptide in pep_psm.keys()])
    print(pep_psm,pep_score)

