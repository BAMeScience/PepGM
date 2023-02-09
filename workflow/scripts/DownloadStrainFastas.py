
from ete3 import NCBITaxa
from Bio import Entrez
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--TaxidFile', required = True, type = str,help = 'path to taxid txt file with on taxid per line')
parser.add_argument('--out', required = True, type = str, help = 'path where output fasta of strins will be saved')
parser.add_argument('--APIkey', type = str, help = 'your NCBI API key if you have one')
parser.add_argument ('--APImail', help = ' e-mail connected to your NCBI API key' )

args = parser.parse_args()

Entrez.email = args.APImail
Entrez.api_key = args.APIkey


def GetAllLeafTaxaFromTaxids(TaxidFile, StrainResolution = True):
        '''
        takes a list of input taxa and add all child taxa into the taxidList attribute
        :input TaxidFile: str, path to tsv file with taxa to be included in graph
        :input StrainResolution, bool, whether to build the graph with strain resolution or not
        :output: none
        '''
        
        with open(TaxidFile) as Taxids:
            
            TargetTaxa = Taxids.read().splitlines()
            if 'no match' in TargetTaxa:
                TargetTaxa.remove('no match')

        TaxidList = []
        for HighestTaxid in TargetTaxa:

            HighestTaxid = int(HighestTaxid)
            ncbi = NCBITaxa()

    
            TargetTaxon = ncbi.get_taxid_translator([HighestTaxid])[HighestTaxid]
            
            if StrainResolution:
                TaxidList.extend(ncbi.get_descendant_taxa(TargetTaxon))
    
            else:
                TaxidList.extend(ncbi.get_descendant_taxa(TargetTaxon, collapse_subspecies=True))
                TaxidNames= ncbi.translate_to_names(TaxidList)
            
            
        return TaxidList


def FetchTaxonDataToFasta(TaxidList,FastaPath, sourceDB):
        '''
        gets the proteins corresponding to the target taxa from entrez and saves them
        in a json document.
        :input PeptideMapPath: str, path to resultsfile
        :input sourceDB: str, describes which DB ncbi will search for the taxon protein records
        :output: fasta file with all strain specific sequences

        '''
    
        entrezDbName = sourceDB

        # options listed here https://www.ncbi.nlm.nih.gov/books/NBK49540/
        #sourceDBOptions = ['refseq[filter]','swissprot[filter]','protein_all[PROP]'] 
        
        saveFastas= str()
        for Taxid in TaxidList:

            ncbiTaxId = str(Taxid) 

            # Find entries matching the query (only swissprot registered proteins for now)
            entrezQuery = sourceDB +' AND txid%s[ORGN]'%(ncbiTaxId) # AND 
            searchResultHandle = Entrez.esearch(db=entrezDbName, term=entrezQuery, retmax = 1000)
            searchResult = Entrez.read(searchResultHandle)
            searchResultHandle.close()

            #fetch corresponding proteins from NCBI entrez
            if searchResult['Count'] != '0':

                uidList = ','.join(searchResult['IdList'])
                fastaEntry = (Entrez.efetch(db=entrezDbName, id=uidList, rettype='fasta',retmax = 1000).read())
                saveFastas+= fastaEntry

                

        with open(FastaPath, 'w+') as savefile:        
            savefile.write(saveFastas)


Taxidlist = GetAllLeafTaxaFromTaxids(args.TaxidFile)
FetchTaxonDataToFasta(Taxidlist,args.out,sourceDB = 'protein')