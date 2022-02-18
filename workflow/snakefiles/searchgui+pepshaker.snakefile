

rule AddContaminantsandHost:
     input:
          DatabaseDirectory+'crap.fasta',
          DatabaseDirectory+'{hostname}.fasta',
          DatabaseDirectory+'{DBname}.fasta'

     output: DatabaseDirectory+ '{hostname}+crap+{DBname}.fasta'
     shell:'cat {input} > {output}'

rule RemoveDuplicates:
     input: DatabaseDirectory+ '{hostname}+crap+{DBname}.fasta'
     output: DatabaseDirectory+ '{hostname}+crap+{DBname}_UNI.fasta'
     conda: 'envs/graphenv.yml'
     shell:  'seqkit rmdup -s {input} > {output}'

rule AddDecoys:
     input: DatabaseDirectory+'{hostname}+crap+{DBname}_UNI.fasta'
     output: DatabaseDirectory+'{hostname}+crap+{DBname}_UNI_concatenated_target_decoy.fasta'
     shell: 'java -cp '+SearchGUIDir+'SearchGUI-4.1.1.jar eu.isas.searchgui.cmd.FastaCLI -in {input} -decoy' 

rule SearchSpectra:
     input:
          DataDirectory+'{samplename}/{samplename}'+SpectraFileType, 
          DatabaseDirectory+'{hostname}+crap+{DBname}_UNI_concatenated_target_decoy.fasta',
          DataDirectory+'{samplename}/{samplename}.par'
     params:
          samplename = SampleName,
          hostname = HostName,
          DBname = ReferenceDBName
     output:  ResultsDir+'{samplename}/{hostname}_{DBname}_searchgui_out.zip'
     shell: 'java -cp '+SearchGUIDir+'SearchGUI-4.1.1.jar eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input[0]} -fasta_file {input[1]} -output_folder '+ ResultsDir +'{params.samplename} -id_params {input[2]} -output_default_name {params.hostname}_{params.DBname}_searchgui_out -psm_fdr '+psmFDR+' -peptide_fdr '+peptideFDR+' -protein_fdr '+proteinFDR+' '+searchengines+' 1'

rule RunPeptideShaker:
     input:
          ResultsDir+'{samplename}/{hostname}_{DBname}_searchgui_out.zip',
          DataDirectory+'{samplename}/{samplename}'+SpectraFileType, 
          DatabaseDirectory+'{hostname}+crap+{DBname}_UNI_concatenated_target_decoy.fasta',
     params:
          samplename = SampleName,
          hostname = HostName,
          DBname = ReferenceDBName
     output: ResultsDir+'{samplename}/{hostname}_{DBname}.psdb'
     shell: 'java -cp '+PeptideShakerDir+'PeptideShaker-2.1.1.jar eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.hostname}_{params.DBname} -fasta_file {input[2]} -identification_files {input[0]} -spectrum_files {input[1]} -out {output}'

rule SimplePeptideList:
     input:  ResultsDir+'{samplename}/{hostname}_{DBname}.psdb'
     output: ResultsDir +'{samplename}/{hostname}_{DBname}_Default_PSM_Report.txt'
     params:
          samplename = SampleName,
          hostname = HostName,
          DBname = ReferenceDBName
     shell: 'java -cp '+PeptideShakerDir+'PeptideShaker-2.1.1.jar eu.isas.peptideshaker.cmd.ReportCLI -in {input} -out_reports '+ResultsDir +'{params.samplename} -reports 3'