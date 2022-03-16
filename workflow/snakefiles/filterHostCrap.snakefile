#filtering out the host and crap spectra before performing the actual searchDecoyDatabase

rule AddHostandCrap:
    input:
        DatabaseDirectory+'crap.fasta',
        DatabaseDirectory+'{hostname}.fasta'

    output: DatabaseDirectory+ '{hostname}+crap.fasta'

    shell:'cat {input} > {output}'  

rule RemoveDuplicatesHostandCrap:

    input: DatabaseDirectory+ '{hostname}+crap.fasta'
    output: DatabaseDirectory+ '{hostname}+crap_UNI.fasta'
    conda: 'envs/graphenv.yml'
    shell:  'seqkit rmdup -s {input} > {output}'

rule SearchSpectraHostandCrap:
     input:
          DataDirectory+'{samplename}/{samplename}'+SpectraFileType, 
          DatabaseDirectory+'{hostname}+crap_UNI_concatenated_target_decoy.fasta',
          DataDirectory+'{samplename}/{samplename}.par'
     params:
          samplename = SampleName,
          hostname = HostName,
     output:  ResultsDir+'{samplename}/SpectraFilter/{hostname}_searchgui_out.zip'
     shell: 'java -cp '+SearchGUIDir+'SearchGUI-4.1.1.jar eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input[0]} -fasta_file {input[1]} -output_folder '+ ResultsDir +'{params.samplename} -id_params {input[2]} -output_default_name {params.hostname}_searchgui_out -psm_fdr '+psmFDR+' -peptide_fdr '+peptideFDR+' -protein_fdr '+proteinFDR+' '+searchengines+' 1'

rule RunPeptideShakerHostandCrap:
     input:
          ResultsDir+'{samplename}/SpectraFilter/{hostname}_searchgui_out.zip',
          DataDirectory+'{samplename}/{samplename}'+SpectraFileType, 
          DatabaseDirectory+'{hostname}+crap_UNI_concatenated_target_decoy.fasta',
     params:
          samplename = SampleName,
          hostname = HostName
     output: ResultsDir+'{samplename}/SpectraFilter/{hostname}.psdb'
     shell: 'java -cp '+PeptideShakerDir+'PeptideShaker-2.1.1.jar eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.hostname} -fasta_file {input[2]} -identification_files {input[0]} -spectrum_files {input[1]} -out {output}'

rule SimplePeptideList:
     input:  ResultsDir+'{samplename}/SpectraFilter/{hostname}.psdb'
     output: ResultsDir +'{samplename}/SpectraFilter/{hostname}_Default_PSM_Report.txt'
     params:
          samplename = SampleName,
          hostname = HostName

     shell: 'java -cp '+PeptideShakerDir+'PeptideShaker-2.1.1.jar eu.isas.peptideshaker.cmd.ReportCLI -in {input} -out_reports '+ResultsDir +'{params.samplename} -reports 3'

rule FilterSpectra:
    input: 
        ResultsDir +'{samplename}/SpectraFilter/{hostname}_Default_PSM_Report.txt',
        DataDirectory+'{samplename}/{samplename}'+SpectraFileType
    output: ResultsDir +'{samplename}/SpectraFilter/{hostname}'+SpectraFileType
    conda: 'envs/graphenv.yml'
    shell: "python3 workflow/scripts/filterHostSpectra.py --SpectrumFile {input[1]} --PSMReport {input[0]} --out {output}"
