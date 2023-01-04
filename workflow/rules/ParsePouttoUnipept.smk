###rule for parsing percolator file

rule UnipeptQuery:
    input: MS2RescoreDir+'rescored_searchengine_ms2pip_rt_features.pout'
    params: 
          NumberofTaxa = TaxaNumber,
          targetTaxa = targetTaxa,
          FDR = FDR
    log: ResultsDir + 'UnipeptResponse.log'
    output: 
          ResultsDir + 'UnipeptResponse.json',
          ResultsDir + 'UnipeptPeptides.json'
    conda: 'envs/graphenv.yml'   
    shell: "python3 workflow/scripts/UnipeptGetTaxonomyfromPout.py --UnipeptResponseFile {output[0]} --pep_out {output[1]} --TaxonomyQuery {params.targetTaxa} --NumberOfTaxa {params.NumberofTaxa} --FDR {params.FDR} --PoutFile {input} " 

rule ParseToUnipeptCSV:
    input: 
          ResultsDir + 'UnipeptResponse.json',
          ResultsDir + 'UnipeptPeptides.json'
          
    params: NumberofTaxa = TaxaNumber
    log: ResultsDir + 'ParsetoCSV.log'
    output: ResultsDir + 'GraphDataframe.csv'
    conda: 'envs/graphenv.yml' 
    shell: "python3 workflow/scripts/WeightTaxa.py --UnipeptResponseFile {input[0]} --UnipeptPeptides {input[1]} --out {output} --NumberOfTaxa {params.NumberofTaxa}" 