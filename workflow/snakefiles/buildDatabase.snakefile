
rule splitToAccessions:
     input: "resources/taxidMapping/protacc2taxids_virus.txt"
     output: "resources/taxidMapping/accessions.txt"
     shell:
        "awk '{{print $1}}' {input} > {output}"


rule splitToTaxids:
     input:
          "resources/taxidMapping/protacc2taxids_virus.txt"
     output:
          "resources/taxidMapping/taxids.txt"
     shell:
          "awk '{{print $2}}' {input} > {output}"


rule hashDatabase:
     input:
          "resources/taxidMapping/accessions.txt"
     output:
          "resources/taxidMapping/accessions_hashed.npy"
     shell:
          "python3 workflow/scripts/hashDatabase.py --i {input}  --o {output}"
