# Prerequisite: fasta files (already available for PepGM itself) + config json file for fasta prediction
# Tools that are used: ms2pip
WorkflowPath = os.path.dirname(os.path.realpath(workflow.snakefile))
include: WorkflowPath + '/config.smk'
include: WorkflowPath + '/searchgui+pepshaker.smk'
include: WorkflowPath + '/buildDatabase.smk'
include: WorkflowPath + '/filterHostCrap.smk'
include: WorkflowPath + '/protGM.smk'

from math import ceil


def num_chunks(fasta_file, chunk_size):
    """
    Compute the number of chunks a FASTA file will be split into based on
    the number of sequences and a specified chunk size.
    """
    with open(fasta_file,'r') as f:
        sequence_count = sum(1 for line in f if line.startswith('>'))

    return ceil(sequence_count / chunk_size)


# Define function to create output filename pattern.
def create_output_pattern(num):
    return f"resources/spectralLibrary/RefSeqViral_chunk{num}.fa"


rule all:
    input:
        expand(create_output_pattern("{num}"),num=range(1,num_chunks(DatabaseDirectory + ReferenceDBName + '.fasta',chunk_size) + 1))


rule split_fasta:
    input:
        DatabaseDirectory + ReferenceDBName + '.fasta'
    output:
        chunks=create_output_pattern("{num}")
    params:
        chunk_size=chunk_size  # Use the constant here.
    shell:
        """
        awk 'BEGIN {{n_seq=1; file=sprintf("{output.chunks}", n_seq); chunk_size={params.chunk_size}; current_count=0}} /^>/ {{if(current_count++ % chunk_size == 0 && current_count > 1) file=sprintf("{output.chunks}", ++n_seq);}} {{print >> file;}}' < {input}
        """


rule predict_spectral_library:
