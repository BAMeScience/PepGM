import re
import networkx as nx
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import Bio.SeqIO
import String_Load 

#This script creates the protein-peptide graph with the protein connected through an interaction network and the peptides initialized wit probabilities from cp-dt peptide/protein is a node attribute, as is petide nae, protein name, protein probabilities etc
#my directories
StringDirectory = "/home/tholstei/VirusGraph/Data/String_files/herpes_sv1.tsv"
FastaDirectory = "/home/tholstei/VirusGraph/Data/Uniprot_sequences/herpes_protein_sequences.fa"
DigestDirectory = 'herpes_digest.fa'

#load the string graph using module string_load
InteractionGraph = String_Load.StringLoad(StringDirectory)
nx.set_node_attributes(InteractionGraph,'protein','type')




#add peptides cleaved by cp-dt into graph, each connected to its protein origin
#get list of protein IDs from the fasta file as they don't show up in cp-dt

with open(FastaDirectory) as FastaFile:
   ProteinIDlist = []
   for FastaSequences in Bio.SeqIO.parse(open(FastaDirectory),"fasta"): #exchange this with simple fasta parser later
      ProteinIDlist.append(FastaSequences.id)
   

#cp-dt digest to data frame, with protein ID's as identifier

#read cp-dt digest petides and score into networkx graph
#DFrows = []                                                              only necessary if we want petide data dataframe, as are al other lines that were commented
with open(DigestDirectory) as ProteinDigest:
   IDindex = 0
   for LineNumber, Line in enumerate(ProteinDigest):
      if Line.find('PEPTIDE') == -1 and Line != ('\n'):                 #find header lines that do not include "PEPTIDE" and aren't just an empty line, assign protein ID
         ParentProt = 'Protein:'+ProteinIDlist[IDindex]
         IDindex +=1 
      if Line.find('PEPTIDE')!=-1:                                      #put peptides into Graph
         #Line = Line.replace('PEPTIDE','')
         Line = re.split(' |:|\n',Line)
         #Line.insert(0,ParentProt)
         #Line = [entry for entry in Line if entry] 
         if float(Line[4])>0.05:                                                                #filter out peptides that have less than 0.05 chance of appearing
            if ParentProt in InteractionGraph:                                                  #extra if makes sure that only proteins other than those in the string database are added for now                                                                     
               if str(Line[2]) not in InteractionGraph:                                                #check if peptide already has other parent. 
                  InteractionGraph.add_node(str(Line[2]),score=float(Line[4]),type= 'peptide')       #if not, add node & score
               else: 
                  score = np.maximum(float(Line[4]),InteractionGraph.nodes[str(Line[2])]['score'])     #if yes, update score to max(score(already existing node),score( new node))
                  nx.set_node_attributes(InteractionGraph,{str(Line[2]):{'score':score}})
               InteractionGraph.add_edge(Line[2],ParentProt, type = 'protein-peptide')                                          
         
         #DFrows.append(Line)
       
#PepNodesDataframe = pd.DataFrame(DFrows, columns=['ProteinID', 'Peptide', 'score'])

#save the graph to test the markov model on it
nx.write_gexf(InteractionGraph,'Interactiongraph.gexf')




   


