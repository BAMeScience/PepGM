#rules to produce files needed when no host&crap filtering is performed, but host&crap are to be added to the search
rule AddContaminantsandHostFull:
    input:
        DatabaseDirectory+'crap.fasta',
        DatabaseDirectory+HostName+'.fasta',
        DatabaseDirectory+ReferenceDBName+'.fasta'

    output: DatabaseDirectory+ HostName+'+crap+'+ReferenceDBName+'.fasta'
    shell:'cat {input} > {output}'

rule RemoveDuplicatesFull:
    input: DatabaseDirectory+ HostName+'+crap+'+ReferenceDBName+'.fasta'
    output: DatabaseDirectory+ HostName+'+crap+'+ReferenceDBName+'_UNI.fasta'
    conda: 'envs/graphenv.yml'
    shell:  'seqkit rmdup -s {input} > {output}'

rule AddDecoysFull:
    input: DatabaseDirectory+ HostName+'+crap+'+ReferenceDBName+'_UNI.fasta'
    output: DatabaseDirectory + HostName+'+crap+'+ReferenceDBName+'_UNI_concatenated_target_decoy.fasta'
    shell: 'java -cp ' + SearchGUI + ' eu.isas.searchgui.cmd.FastaCLI -in {input} -decoy'


#rules to produce files necessary for searching after filtering host spectra or to search all spectra but whithout host or crap DB added
rule RemoveDuplicates:
    input: DatabaseDirectory+ ReferenceDBName+'.fasta'
    output: DatabaseDirectory+ReferenceDBName+'_UNI.fasta'
    conda: 'envs/graphenv.yml'
    shell:  'seqkit rmdup -s {input} > {output}'



rule AddDecoys:
    input: DatabaseDirectory + ReferenceDBName+'_UNI.fasta'
    output: DatabaseDirectory + ReferenceDBName+'_UNI_concatenated_target_decoy.fasta'
    shell: 'java -cp ' + SearchGUI + ' eu.isas.searchgui.cmd.FastaCLI -in {input} -decoy'


# check if spectrum should be filtered or not
def SpectrumToUse(condition):
    if condition:
        return ResultsDir + SampleName + '/SpectraFilter/Filtered_' + HostName + '.mgf'
    else:
        return SamplePath


#if the spectra aren't beeing filtered, check whether host and crap should be added to the search DB
def DBToUse(condition):
    if condition:
        return DatabaseDirectory + HostName+'+crap+'+ReferenceDBName+'_UNI_concatenated_target_decoy.fasta'
    else:
        return DatabaseDirectory + ReferenceDBName+'_UNI_concatenated_target_decoy.fasta'


InputSpectrum = SpectrumToUse(FilterSpectra)
InputDB = DBToUse(AddHostandCrapToDB)

rule SearchSpectra:
    input:
        #ResultsDir+SampleName+'/SpectraFilter/Filtered_'+HostName+SpectraFileType,
        InputSpectrum,
        InputDB,
        ParametersFile

    params:
        samplename=SampleName,
        hostname=HostName,
        ResultsDir=ResultsDir,
        DBname=ReferenceDBName
    output: ResultsDir + SampleName + '/searchgui_out.zip'
    shell: 'cp config/config.yaml {params.ResultsDir}{params.samplename}/  &&  java -cp ' + SearchGUI + ' eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input[0]} -fasta_file {input[1]} -output_folder ' + ResultsDir + '{params.samplename} -id_params {input[2]} -output_default_name searchgui_out -psm_fdr ' + psmFDR + ' -peptide_fdr ' + peptideFDR + ' -protein_fdr ' + proteinFDR + ' ' + searchengines + ' 1'

rule RunPeptideShaker:
    input:
        ResultsDir + SampleName + '/searchgui_out.zip',
        InputSpectrum,
        InputDB
    params:
        samplename=SampleName,
        DBname=ReferenceDBName
    output: ResultsDir + SampleName + '/PepShaker.psdb'
    shell: 'java -cp ' + PeptideShaker + ' eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference Out -fasta_file {input[2]} -identification_files {input[0]} -spectrum_files {input[1]} -out {output}'

rule SimplePeptideList:
    input: ResultsDir + SampleName + '/PepShaker.psdb'
    output: ResultsDir + SampleName + '/Out_Default_PSM_Report.txt'
    params:
        samplename=SampleName,
        DBname=ReferenceDBName
    shell: 'java -cp ' + PeptideShaker + ' eu.isas.peptideshaker.cmd.ReportCLI -in {input} -out_reports ' + ResultsDir + '{params.samplename} -reports 3'