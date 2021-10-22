import pandas as pd

#scripts to get all the taxa potentially associated to peptides in my results file. 
import argparse
from LoadMzID import loadfoundProteins

parser = argparse.ArgumentParser(description='gets a list of all taxa present in the peptideShaker output File')

parser.add_argument('--input', help ='input list of found protein accession')
parser.add_argument('--output', help = 'output filepath')
parser.add_argument('--taxonLevel', help = 'taxonomic ranks we want to return (i.e strain/species/family)')

args = parser.parse_args()

# idea from https://bionerdnotes.wordpress.com/2020/03/29/getting-the-taxid-from-sequence-accessions
#but later use index/proper database!!
#!/bin/csh
set file='inputfile.txt'
set library='prot.accession2taxid'
#echo $file
#echo $library
foreach accession (`awk '{print $1}' "$file"`)
    #echo $accession
    grep -m 1 -h "$accession" "$library" | head -1
end