rule catDatabase:
     input: 
          ResourcesDir + TaxidMapping + "protacc2taxids00",
          ResourcesDir + TaxidMapping + "protacc2taxids01",
          ResourcesDir + TaxidMapping + "protacc2taxids02"

     output: ResourcesDir + TaxidMapping + "protacc2taxids_virus.txt"

     shell: "cat > protacc2taxids_virus.txt"


rule splitToAccessions:
     input: ResourcesDir + TaxidMapping + "protacc2taxids_virus.txt"
     output: ResourcesDir + TaxidMapping + "accessions.txt"
     shell:
        "awk '{{print $1}}' {input} > {output}"


rule splitToTaxids:
     input:
          ResourcesDir + TaxidMapping + "protacc2taxids_virus.txt"
     output:
          ResourcesDir + TaxidMapping + "taxids.txt"
     shell:
          "awk '{{print $2}}' {input} > {output}"


rule hashDatabase:
     input:
        ResourcesDir + TaxidMapping + "accessions.txt"
     output:
          ResourcesDir + TaxidMapping + "accessions_hashed.npy"
     conda: 'envs/graphenv.yml'
     shell:
          "python3 workflow/scripts/hashDatabase.py --i {input}  --o {output}"
