rule SimpleProteinList:
    input: ResultsDir + SampleName + '/{DBname}.psdb'
    output: ResultsDir + SampleName + '/{DBname}_Default_Protein_Report.txt'
    params:
        samplename=SampleName,
        DBname=ReferenceDBName
    shell: 'java -cp ' + PeptideShaker + ' eu.isas.peptideshaker.cmd.ReportCLI -in {input} -out_reports ' + ResultsDir + '{params.samplename} -reports 9'

rule GetProteinAccessions:
    input: ResultsDir + SampleName + '/{DBname}_Default_Protein_Report.txt'
    output:
           ResultsDir + SampleName + '/{DBname}_Protein_accessions.txt',
           ResultsDir + SampleName + '/{DBname}_Protein_accessions_scored.txt'
    conda: 'envs/graphenv.yml'
    shell: 'python3 workflow/scripts/parseProteinReport.py --proteinfile {input[0]} --out {output[0]} --out_score {output[1]}'

rule SearchAccessionsinRefseq:
    input: 
        ResultsDir + SampleName + '/{DBname}_Protein_accessions.txt',
        DatabaseDirectory + '{DBname, [A-Za-z]+}.fasta'

    output: ResultsDir + SampleName + '/{DBname}_Found_proteins.fasta'
    conda: 'envs/graphenv.yml'
    shell: 'seqkit grep -f {input[0]} {input[1]} -o {output}'

rule getStrainFasta:
    input: ResultsDir + SampleName + '/{DBname}_mapped_taxids.txt'
    output: ResultsDir + SampleName + '/{DBname}_AllStrains.fasta'
    params:
        APImail=APImail,
        APIkey=APIkey
    conda: 'envs/graphenv.yml'
    shell: 'python3 workflow/scripts/DownloadStrainFastas.py --TaxidFile {input} --out {output} --APIkey {params.APIkey} --APImail {params.APImail}'

rule createViralStrainBlastDB:
    input: ResultsDir + SampleName + '/{DBname}_Found_proteins.fasta'
    output: ResultsDir + SampleName + '/{DBname, [A-Za-z]+}.psq'
    conda: 'envs/blast.yml'
    shell: 'makeblastdb -in {input} -title ViralblastDB -dbtype prot -out ' +ResultsDir + SampleName + '/{wildcards.DBname}'

rule BlastFoundProteins:
   input: 
        ResultsDir + SampleName + '/{DBname}_AllStrains.fasta',
        ResultsDir + SampleName + '/{DBname}.psq'
   output: ResultsDir + SampleName + '/{DBname}_blasted_prots.txt'
   conda: 'envs/blast.yml'
   shell: 'blastp -query {input[0]} -db '+ResultsDir + SampleName + '/{wildcards.DBname} -out {output} -num_threads 4 -outfmt 6'

rule ProteinGraphCSV:
    input:
        ResultsDir + SampleName + '/{DBname}_blasted_prots.txt',
        ResultsDir + SampleName + '/{DBname}_mapped_taxids_accessions.csv',
        ResultsDir + SampleName + '/{DBname}_Protein_accessions_scored.txt'
    output: ResultsDir + SampleName + '/{DBname}_protein_graph.csv'
    conda: 'envs/graphenv.yml'
    shell: 'python3 workflow/scripts/ProtGMCSVcreation.py --TaxidAccessionMap {input[1]} --Blastp {input[0]} --ProteinScores {input[2]} --out {output}'

     

