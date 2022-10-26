###rule for parsing percolator file

rule ParseToUnipeptCSV:
    input: MS2RescoreDir+'rescored_searchengine_ms2pip_rt_features.pout'
    params: 
          NumberofTaxa = TaxaNumber,
          targetTaxa = targetTaxa,
          MinScore = MinPepScore
    output: 
          ResultsDir + 'UnipeptResponse.json',
          ResultsDir + 'GraphDataframe.csv'
    conda: '/home/tholstei/repos/PepGM_all/PepGM/workflow/envs/graphenv.yml'   
    shell: "python3 workflow/scripts/WeightTaxa.py --UnipeptResponseFile {output[0]} --out {output[1]} --TaxonomyQuery {params.targetTaxa} --NumberOfTaxa {params.NumberofTaxa} --MinScore {params.MinScore} --PoutFile {input} " 