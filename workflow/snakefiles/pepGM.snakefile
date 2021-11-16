#all rules relating to the pepGM algorithm itself


rule CreateFactorGraph:
     input: DataDirectory +'{samplename}/{hostname}_{DBname}_Default_PSM_Report.txt'
     output: DataDirectory +'{samplename}/{hostname}_{DBname}_PepGM_graph.graphml'
     conda: WorkflowPath + 'envs/graphenv.yml'
     params:
          samplename = SampleName,
          hostname = HostName,
          DBname = ReferenceDBName,
          targetTaxa = TargetTaxa,
          firstTarget = firstTarget
     shell: 'python3 workflow/scripts/CreatePepGMGraph.py --targetTaxa {params.targetTaxa} --PSM_Report {input} --PeptideMapPath '+DataDirectory+'{params.samplename}/{params.firstTarget}.json --out {output}'


rule RunPepGM:
     input: DataDirectory +'{samplename}/{hostname}_{DBname}_PepGM_graph.graphml'
     output: DataDirectory +'{samplename}/{hostname}_{DBname}_PepGM_Results_{alpha}_{beta}.csv'
     conda: WorkflowPath + 'envs/graphenv.yml'
     params:
          samplename = SampleName,
          hostname = HostName,
          DBname = ReferenceDBName,
          targetTaxa = TargetTaxa,
          firstTarget = firstTarget
     shell: 'python3  workflow/scripts/PepGM.py --targetTaxa {params.targetTaxa} --GraphMLPath {input} --alpha{wildcard.alpha} --beta{wildcard.beta} --out {output}'

