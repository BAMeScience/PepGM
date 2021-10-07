import re
import linecache
import pandas as pd


def loadPepsMzID(filepath):

    pepshaker_substr1 = r'PeptideShaker PSM confidence" value='        # Substring to search for any pepshaker match
    pepshaker_peptidesubstr = r'SpectrumIdentificationItem passThreshold=' # Substring to search for


    pepShakerName = []
    pepShakerScore =[] 

    with open (filepath, 'rt') as pepshaker_lines:
        for line_num, line in enumerate(pepshaker_lines):
            if line.find(pepshaker_peptidesubstr) != -1:    # if case-sensitive match
                pepshaker_peptideline = line_num+1
                #print(peptideline)
        if line.find(pepshaker_substr1) != -1:    # if case-sensitive match,
                print('match2')
                pepshaker_start_pepname = linecache.getline(pepshaker_searchfile, pepshaker_peptideline).index('peptide_ref=')+13                           # find the index of the substring befor the Pepname
                pepshaker_end_pepname = linecache.getline(pepshaker_searchfile, pepshaker_peptideline).index('calculatedMassToCharge=')-2
                pepShakerName.append(linecache.getline(pepshaker_searchfile, pepshaker_peptideline)[pepshaker_start_pepname:pepshaker_end_pepname])    # and use it to cut everything from the line except the pepname
    
                pepshaker_start_confidence = line.index('value=') + 7 # find the index of the substring befor the Pepname
                pepshaker_end_confidence = line.index('/>') - 1
                pepShakerScore.append(line[pepshaker_start_confidence:pepshaker_end_confidence])  # and use it to cut everything from the line except

    return pepShakerName,pepShakerScore



def loadSimplePepScore(filepath):
    pepIDs= pd.read_csv (filepath, sep = '\t', error_bad_lines=False) #error bad lines should be remove when development is done
    pepnames = pepIDs['Sequence'].tolist()
    pepscores = pepIDs['Confidence [%]'].tolist()

    #directly remove IDS with zero confidence
    pepnames = [pepnames[i] for i in range(len(pepscores)) if pepscores[i] != 0]
    pepscores = [pepscore for pepscore in pepscores if pepscore!=0]
    pepscores = [99.99 if pepscore ==100 else pepscore for pepscore in pepscores]


    return pepnames,pepscores



