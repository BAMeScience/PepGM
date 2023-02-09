import pandas as pd
from ete3 import NCBITaxa
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--TaxidAccessionMap',type = str, required = True,help = 'path to taxid-protein accesion in DB mapping file')
parser.add_argument('--Blastp', type = str, required = True, help = 'path to blastp results with output format 6 (tabular)')
parser.add_argument('--ProteinScores', type =str, required = True, help = 'path to csv of protein accessions and their scores')
parser.add_argument('--out', type = str, required = True, help = 'output graph csv path')

args = parser.parse_args()


def GetAllLeafTaxaFromTaxid(Taxid, StrainResolution = True):
        '''
        takes an and add all child taxa into the taxidList attribute
        :input TaxidFile: str, path to tsv file with taxa to be included in graph
        :input StrainResolution, bool, whether to return taxids with strain reslotion or not
        :output: descendant taxids
        '''
        Taxid = int(Taxid)
        ncbi = NCBITaxa()

        TargetTaxon = ncbi.get_taxid_translator([Taxid])[Taxid]
            
        if StrainResolution:
            Taxids = ncbi.get_descendant_taxa(TargetTaxon)
    
        else:
            Taxids = ncbi.get_descendant_taxa(TargetTaxon, collapse_subspecies=True)
            
        return Taxids


def CreateProtCSV(TaxidAccessionMappings,BlastResults,AccessionsScored,out):

    AccessionTaxidFrame = pd.read_csv(TaxidAccessionMappings)
    BlastResults = pd.read_table(BlastResults, names= ['query accession','ref accession','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bit score'])
    ProteinScores = pd.read_csv(AccessionsScored)
    
    AccessionTaxidFrame = AccessionTaxidFrame[AccessionTaxidFrame.taxid != 'no match']
    AccessionTaxidFrame['AllTaxids'] = AccessionTaxidFrame['taxid'].apply(lambda x : GetAllLeafTaxaFromTaxid(x))

    BlastResults = BlastResults.loc[BlastResults['evalue']<=0.01]
    BlastResults = BlastResults.loc[BlastResults['pident']>=90]
    BlastResults = BlastResults[['query accession','ref accession','pident']]
    BlastResults['protein score'] = BlastResults['ref accession'].map(ProteinScores.set_index('All Accessions').to_dict()['Confidence [%]'])
    BlastResults['taxids'] = BlastResults['ref accession'].map(AccessionTaxidFrame.set_index('accession').to_dict()['AllTaxids'])
    BlastResults = BlastResults.explode('taxids')
    BlastResults['score'] = 0.01*BlastResults['pident'].multiply(0.01*BlastResults['protein score'])
    BlastResults.to_csv(out)




#map = '/home/tholstei/repos/PepGM_all/PepGM/results/ProteinGM/sars2/refseqViral_mapped_taxids_accessions.csv'
#blast = '/home/tholstei/repos/PepGM_all/PepGM/results/ProteinGM/sars2/refseqViral_blasted_prots.txt'
#scores = '/home/tholstei/repos/PepGM_all/PepGM/results/ProteinGM/sars2/refseqViral_Protein_accessions_scored.txt'
#out = 'testt.csv'

out = CreateProtCSV(args.TaxidAccessionMap,args.Blastp,args.ProteinScores,args.out)

#(map,blast,scores,out)#