rule RemoveDuplicates:
     input: DatabaseDirectory+ '{DBname}.fasta'
     output: DatabaseDirectory+ '{DBname}_UNI.fasta'
     conda: 'envs/graphenv.yml'
     shell:  'seqkit rmdup -s {input} > {output}'


rule AddDecoys:
     input: DatabaseDirectory+'{DBname}_UNI.fasta'
     output: DatabaseDirectory+'{DBname}_UNI_concatenated_target_decoy.fasta'
     shell: 'java -cp '+SearchGUIDir+'SearchGUI-4.1.1.jar eu.isas.searchgui.cmd.FastaCLI -in {input} -decoy' 


# check if spectrum should be filtered or not
def SpectrumToUse(condition):
     if condition:
          return ResultsDir+SampleName+'/SpectraFilter/Filtered_'+HostName+SpectraFileType
     else:
          return DataDirectory+SampleName+'/'+SampleName+SpectraFileType

InputSpectrum = SpectrumToUse(FilterSpectra)

rule SearchSpectra:
     input:
          #ResultsDir+SampleName+'/SpectraFilter/Filtered_'+HostName+SpectraFileType, 
          InputSpectrum,
          DatabaseDirectory+'{DBname}_UNI_concatenated_target_decoy.fasta',
          DataDirectory+SampleName+'/'+SampleName+'.par'

     params:
          samplename = SampleName,
          hostname = HostName,
          DBname = ReferenceDBName
     output:  ResultsDir+SampleName+'/{DBname}_searchgui_out.zip'
     shell: 'java -cp '+SearchGUIDir+'SearchGUI-4.1.1.jar eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input[0]} -fasta_file {input[1]} -output_folder '+ ResultsDir +'{params.samplename} -id_params {input[2]} -output_default_name {params.hostname}_{params.DBname}_searchgui_out -psm_fdr '+psmFDR+' -peptide_fdr '+peptideFDR+' -protein_fdr '+proteinFDR+' '+searchengines+' 1'

rule RunPeptideShaker:
     input:
          ResultsDir+SampleName+'/{DBname}_searchgui_out.zip',
          DataDirectory+SampleName+'/'+SampleName+SpectraFileType, 
          DatabaseDirectory+'{DBname}_UNI_concatenated_target_decoy.fasta'
     params:
          samplename = SampleName,
          DBname = ReferenceDBName
     output: ResultsDir+SampleName+'/{DBname}.psdb'
     shell: 'java -cp '+PeptideShakerDir+'PeptideShaker-2.1.1.jar eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.DBname} -fasta_file {input[2]} -identification_files {input[0]} -spectrum_files {input[1]} -out {output}'

rule SimplePeptideList:
     input:  ResultsDir+SampleName+'/{DBname}.psdb'
     output: ResultsDir +SampleName+'/{DBname}_Default_PSM_Report.txt'
     params:
          samplename = SampleName,
          DBname = ReferenceDBName
     shell: 'java -cp '+PeptideShakerDir+'PeptideShaker-2.1.1.jar eu.isas.peptideshaker.cmd.ReportCLI -in {input} -out_reports '+ResultsDir +'{params.samplename} -reports 3'