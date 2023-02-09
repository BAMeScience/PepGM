import pandas as pd 
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = 'Run the PepGM algorithm from command line')
parser.add_argument('--proteinfile', type = str, nargs = 1,
                    help ='Path to Peptideshaker protein report')
parser.add_argument('--out', type = str, nargs = 1,
                    help ='Path to output txt file of accessions')
parser.add_argument('--out_score', type = str, nargs = 1,
                    help ='Path to save output file with accessions and score')
args = parser.parse_args()



def ParseProteinReport(ProteinReportFile,AccessionsOutput,ScoreOutput):


    proteinReport = pd.read_table(ProteinReportFile[0])
    proteinReport['Secondary Accessions'] = proteinReport['Secondary Accessions'].apply(lambda x: '' if pd.isna(x) else x)
    proteinReport['All Accessions'] = proteinReport['Main Accession']+ ',' + proteinReport['Secondary Accessions']
    proteinReport['All Accessions'] = proteinReport['All Accessions'].apply(lambda x: x.split(','))
    proteinReport = proteinReport.explode('All Accessions').reset_index()
    proteinReport['All Accessions'].replace('',np.nan, inplace = True)
    proteinReport.dropna(subset = ['All Accessions'], inplace = True)
    
    proteinReport.to_csv(AccessionsOutput[0],columns = ['All Accessions'], header = False, index = False)
    proteinReport.to_csv(ScoreOutput[0], columns = ['All Accessions','Confidence [%]'], header = True, index = False)

ParseProteinReport(args.proteinfile, args.out, args.out_score)
#ParseProteinReport('/home/tholstei/repos/PepGM_all/PepGM/results/VMBenchmark/sars2/refseqViral_Default_Protein_Report.txt', 'testAccessions.txt','Textaccessionsscored.txt')