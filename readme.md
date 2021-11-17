<h2> PGM for Metaproteomic and virus strain identification <\h2>

A snakemake workflow <br>


_____________________________________________________________________________________________________________
<p> Input: peptide and protein IDs and scores from previous DB searches using viral samples, the DB used for the searches, various model parameters<br>
The input parameters need to be specified in config.yaml 
Output: Organism/ taxonomic classification with a confidence score, visualization of input/output data, other outputs to be defined <\p>

<br>

<p> I am using a model similar too https://www.openms.de/comp/epifany/ <br>

<p> next steps: <br>

- how to evaluate the grid search? <br>
- add functions to infer target taxa with possible restrictions to certains hosts by the user -> possibly use taxadb package already available?<br>
- add possibility to map peptides only to higher taxonomic layer/ if we have an an unknown virus, output viral families -> maybe just infer them from the virus strains identified?<br>
- find a way to download proteins for metaproteomic samples / generate a DB locally that stores taxon-peptide relationships? (e.g. was also done for bacterial samples for taxit)<br>
- try virus strain identification including non-curated Proteins<br>
- strains not in the sample do not seem to get their probability reduced <br> how to (maybe) fix this?<br>
- detection & dampening of oscillations, (less urgent) <br>
- look at the sequence similarity of strain that are get high scores (=is this because the strains a re too similar<br>
- fix deprecated method in BarPlotResults.py <br>
- host peptide filering <br>
- use e-values from x-tandem directly instead of the peptideshared re-scored scores, check influence of this on the results <br>
- maybe replace cp-dt with something more efficient in the peptide digestion) <br>




<p> General mixed TODOs: <br>
- change way of acquiring the peptide-protein graph to allow general user inputs. includes: using scores from DBsearch, downloading proteins for given species/user defined target DB<br>
-Marten danens(?spelling) to see about sars cov 2 strain specific data<br>




<p> Ideas sidelined for now/for later:<br>
- pgmpy python library to implement simple version of inference with a noisy Or CPD at each peptide: this doesn't work with pgmpy because there is no way to "add" the proteins together to create the noisyOr and you can't combine the noisy OR and the Markovmodel in pgmpy <br>
- incorporate the connections between co-occuring proteins : use the algorithm proposed in https://advances.sciencemag.org/content/7/17/eabf1211<br>

<p> solved steps log <br>
- restrict update of messages to the parts of the graph that received new info in the last iterations step, until 17.05<br>
- save ad gml and <br>
- generation of taxonomic layer and passing message between taxonomic and petide layer **till 13.06 **. Is there already a tool that assigns proteins to taxa /peptides to taxa? possibility to get taxonomic tree already as  graph (for later visualization) optional: include LCA, set priors to only account for likely present species (?)<br>
- restructure git repo <br>
- write snakemake pipeline to make my life easier<br>
- fix the problem of to little valid models identified in my xtandem searches...how is that possible?<br>