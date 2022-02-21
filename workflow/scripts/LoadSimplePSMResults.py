import pandas as pd


def loadfoundProteins(filepath):
    '''take PeptideShaker default PSM results and retrieve all identified protein accessions into list'''
    pepIDs= pd.read_csv (filepath, sep = '\t', error_bad_lines=False) #error bad lines should be remove when development is done
    proteins_all = pepIDs['Protein(s)'].tolist()
    proteins = []
    for protein in proteins_all:
        protein.split(',')
        proteins.append(protein)
    proteins.extend()

    return proteins
    

def loadSimplePepScore(filepath):
    '''returns confidence value and peptide names from PeptideShaker default PSm results'''
    pepIDs= pd.read_csv (filepath, sep = '\t', error_bad_lines=False) #error bad lines should be remove when development is done
    pepnames = pepIDs['Sequence'].tolist()
    pepscores = pepIDs['Confidence [%]'].tolist()

    #directly remove IDS with zero confidence
    pepnames = [pepnames[i] for i in range(len(pepscores)) if pepscores[i] != 0]
    pepscores = [pepscore for pepscore in pepscores if pepscore!=0]
    pepscores = [99.99 if pepscore ==100 else pepscore for pepscore in pepscores]


    return pepnames,pepscores



