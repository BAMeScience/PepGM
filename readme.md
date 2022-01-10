<h2> PGM for Metaproteomic and virus strain identification <\h2>

A snakemake workflow <br>


_____________________________________________________________________________________________________________
<p> Input: mass spectrum mgf or mzML files, searchGUI parameter file, config file (DB used, target taxa, host)<br>

Output: Organism/ taxonomic classification with a confidence score, visualization of input/output data as phylogenetic tree and barplots<\p>

<br>

<p> The graphical model used was inspired by https://www.openms.de/comp/epifany/ <br>

<p>Workflow description:<br>
1. searchDB cleanup : cRaP DB ist added, host is added (if wanted), duplicate entries are removed using seqkit. generation of target-decoy DB using searchCLI<br>
2. peptide search using searchCLI + PeptideShaker. Generation of a a peptide list<br>
3. all descendant strain of the target taxa are queried in the NCBI protein DB (possibility to filter swissprot only/all/refseq onlyetc) through the NCBI API. scripts: CreatePepGMGraph.py and FactorGraphGeneration.py<br>
4. Donwloaded protein recordes are digested using cp-dt and queried again the protein ID list to generate a bipartite taxon-peptide graph scripts: CreatePepGMGraph.py and FactorGraphGeneration.py<br>
5. The bipartite graph is transformed into a factor graph using convolution trees and conditional probability table factors (CPD). scripts: CreatePepGMGraph.py and FactorGraphGeneration.py<br>
6. for different sets of CPD parameters, the belief propagation algorithm is run until convergence to obtain the posterior probabilites of the taxa. scripts: belief_propagation.py and PepGM.py <br>
7. Through an  empirically deduced metric, the ideal parameter set is inferred. script GridSearchAnalysis.py <br>
8. For this ideal parameter set, we output a results barchart and phylogenetic tree view showcasing the 15 best scoring tax. scripts: BarPlotResults, PhyloTreeView.py<br>

<p> next steps: <br>

- add functions to infer target taxa with possible restrictions to certains hosts by the user -> possibly use taxadb package already available?<br>
- find a way to download proteins for metaproteomic samples / generate a DB locally that stores taxon-peptide relationships? (e.g. was also done for bacterial samples for taxit)<br>
- detection & dampening of oscillations, (less urgent) <br>
- look at the sequence similarity of strain that are get high scores (=is this because the strains a re too similar<br>
- fix deprecated method in BarPlotResults.py <br>
- use e-values from x-tandem directly instead of the peptideshared re-scored scores, check influence of this on the results <br>
- maybe replace cp-dt with something more efficient in the peptide digestion -> instead, try using the cp-dt scores as values of the alpha parameter <br>





<p> Ideas sidelined for now/for later:<br>
- incorporate the connections between co-occuring proteins : use the algorithm proposed in https://advances.sciencemag.org/content/7/17/eabf1211<br>
