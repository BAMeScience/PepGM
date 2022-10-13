#rules to produce files needed when no host&crap filtering is perdormed, but host&crap are to be added to the search
rule AddContaminantsandHostFull:
    input:
        DatabaseDirectory + 'crap.fasta',
        DatabaseDirectory + HostName + '.fasta',
        DatabaseDirectory + '{DBname}.fasta'

    output: DatabaseDirectory + HostName + '+crap+{DBname}.fasta'
    shell: 'cat {input} > {output}'

rule RemoveDuplicatesFull:
    input: DatabaseDirectory + HostName + '+crap+{DBname,[A-Za-z]+}.fasta'
    output: DatabaseDirectory + HostName + '+crap+{DBname}_UNI.fasta'
    conda: 'envs/graphenv.yml'
    shell: 'seqkit rmdup -s {input} > {output}'

rule AddDecoysFull:
    input: DatabaseDirectory + HostName + '+crap+{DBname,[A-Za-z]+}_UNI.fasta'
    output: DatabaseDirectory + HostName + '+crap+{DBname}_UNI_concatenated_target_decoy.fasta'
    shell: 'java -cp ' + SearchGUIDir + 'SearchGUI-4.1.14.jar eu.isas.searchgui.cmd.FastaCLI -in {input} -decoy'


#rules to produce files necessary for searching after filtering host spectra or to search all spectra but whithout host or crap DB added
rule RemoveDuplicates:
    input: DatabaseDirectory + '{DBname,[A-Za-z]+}.fasta'
    output: DatabaseDirectory + '{DBname,[A-Za-z]+}_UNI.fasta'
    conda: 'envs/graphenv.yml'
    shell: 'seqkit rmdup -s {input} > {output}'


rule AddDecoys:
    input: DatabaseDirectory + '{DBname}_UNI.fasta'
    output: DatabaseDirectory + '{DBname,^.[a-zA-Z]$}_UNI_concatenated_target_decoy.fasta'
    shell: 'java -cp ' + SearchGUIDir + 'SearchGUI-4.1.14.jar eu.isas.searchgui.cmd.FastaCLI -in {input} -decoy'


# check if spectrum should be filtered or not
def SpectrumToUse(condition):
    if condition:
        return ResultsDir + SampleName + '/SpectraFilter/Filtered_' + HostName + '.mgf'
    else:
        return SamplePath


#if the spectra aren't beeing filtered, check whether host and crap should be added to the search DB
def DBToUse(condition):
    if condition:
        return DatabaseDirectory + HostName + '+crap+{DBname}_UNI_concatenated_target_decoy.fasta'
    else:
        return DatabaseDirectory + '{DBname}_UNI_concatenated_target_decoy.fasta'


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
    output: ResultsDir + SampleName + '/{DBname}_searchgui_out.zip'
    shell: 'cp config/config.yaml {params.ResultsDir}/{params.samplename}/  &&  java -cp ' + SearchGUIDir + 'SearchGUI-4.1.14.jar eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input[0]} -fasta_file {input[1]} -output_folder ' + ResultsDir + '{params.samplename} -id_params {input[2]} -output_default_name {params.DBname}_searchgui_out -psm_fdr ' + psmFDR + ' -peptide_fdr ' + peptideFDR + ' -protein_fdr ' + proteinFDR + ' ' + searchengines + ' 1'

rule RunPeptideShaker:
    input:
        ResultsDir + SampleName + '/{DBname}_searchgui_out.zip',
        InputSpectrum,
        InputDB
    params:
        samplename=SampleName,
        DBname=ReferenceDBName
    output: ResultsDir + SampleName + '/{DBname}.psdb'
    shell: 'java -cp ' + PeptideShakerDir + 'PeptideShaker-2.2.9.jar eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.DBname} -fasta_file {input[2]} -identification_files {input[0]} -spectrum_files {input[1]} -out {output}'

rule SimplePeptideList:
    input: ResultsDir + SampleName + '/{DBname}.psdb'
    output: ResultsDir + SampleName + '/{DBname}_Default_PSM_Report.txt'
    params:
        samplename=SampleName,
        DBname=ReferenceDBName
    shell: 'java -cp ' + PeptideShakerDir + 'PeptideShaker-2.2.9.jar eu.isas.peptideshaker.cmd.ReportCLI -in {input} -out_reports ' + ResultsDir + '{params.samplename} -reports 3'
