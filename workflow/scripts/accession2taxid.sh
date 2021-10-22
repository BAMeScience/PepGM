#!/bin/sh
set file='$1'
set library='/home/tholstei/repos/PepGM_all/PepGM/resources/ncbi_dumpfiles/prot.accession2taxid'
#echo $file
#echo $library
for accession in awk '{print $1}' "$file"
do
    #echo $accession
    head -1 accession
    head -1 library
    grep -m 1 -h "$accession" "$library" | head -1
done