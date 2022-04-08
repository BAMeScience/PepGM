''' functions and classes to hold and generate the graphs for the belief propagation algorithm'''
#implementation of belief propagation on a peptide-protein graph
#__________________________________________________________________________________________
import numpy as np
import networkx as nx
from collections import namedtuple 
import pandas as pd
import Bio.SeqIO
import re
import pandas as pd
import json
import subprocess
from ete3 import NCBITaxa
from Bio import Entrez
from LoadSimplePSMResults import loadSimplePepScore

#e-mail for fetching viral protein sequences from entrez 
Entrez.email = 'tanja.holstein@bam.de'
Entrez.tool = 'PepGM'
Entrez.api_key = 'cbd119d8a7d988bc27f5b0a722a7c861f408'



#peptides from cp-dt need to have minimum lenngth <peplength> to be included in the graph
peplength = 3


def normalize(Array):
    normalizedArray = Array/np.sum(Array)
    return normalizedArray

class ProteinPeptideGraph(nx.Graph):
    
    #class for stroing the protein-peptide graph with scores and priors but no factors, e.g for visual representation
    
    def ProteinsFromStringData(self,StringDirectory,pProteinprior):
        StringDirectory = StringDirectory
    
        StringDataFile = open(StringDirectory)
        StringData = pd.read_csv(StringDataFile, delimiter = "\t" )
        Interactions = StringData[['#node1', 'node2', 'combined_score']]  

        # networkx graph 
        Interactions = np.array(Interactions)
        for i in range(len(Interactions)):
           interaction = Interactions[i]
           nodeA = interaction[0]
           nodeB = interaction[1]
           weight = float(interaction[2]) # score as weighted edge where high scores = low weight
           self.add_weighted_edges_from([(nodeA,nodeB,weight)], category = 'protein-protein') # add weighted edge to graph

        nx.set_node_attributes(self,np.asarray([1-pProteinprior,pProteinprior]), name ='InitialBelief')
        nx.set_node_attributes(self,'protein',name = 'category')
        return self
    
    def PeptidesFromCPDT(self,FastaDirectory,DigestDirectory,threshold):
        #read peptide from an output file of CP-DT digested graph and fasta file that was used to generate cp-dt output  
        ProteinIDlist = []
        for FastaSequences in Bio.SeqIO.parse(open(FastaDirectory),"fasta"):
            ProteinIDlist.append(FastaSequences.id)

        with open(DigestDirectory) as ProteinDigest:
           IDindex = 0
           for LineNumber, Line in enumerate(ProteinDigest):
              if Line.find('PEPTIDE') == -1 and Line != ('\n'):                 #find header lines that do not include "PEPTIDE" and aren't just an empty line, assign protein ID
                 ParentProt = ProteinIDlist[IDindex]
                 IDindex +=1 
              if Line.find('PEPTIDE')!=-1:                                      #put peptides into Graph
                 Line = re.split(' |:|\n',Line)
                 if float(Line[4])>threshold and len(str(Line[2]))>peplength:                                                                #filter out peptides that have less than 0.05 chance of appearing
                    if ParentProt in self:                                                  #extra if makes sure that only proteins other than those in the string database are added for now                                                                     
                       if str(Line[2]) not in self:                                                #check if peptide already has other parent. 
                          self.add_node(str(Line[2]),InitialBelief = [1-float(Line[4]),float(Line[4])], category = 'peptide')       #if not, add node & score
                       else: 
                          score = np.maximum(float(Line[4]),self.nodes[str(Line[2])]['InitialBelief'][1])     #if yes, update score to max(score(already existing node),score( new node))
                          nx.set_node_attributes(self,{str(Line[2]):{'InitialBelief_0':1-score,'InitialBelief_1':score}})
                       self.add_edge(Line[2],ParentProt, category = 'protein-peptide')     
        return self

class TaxonGraph(nx.Graph):
    '''
    class with functions to construct a peptide-taxon graph using entrez/ncbi mapping
    to consider for the future: local DB versions!!
    '''
    def __init__(self):
        nx.Graph.__init__(self)
        self.TaxidList = []

    def GetAllLeafTaxa(self,TargetTaxa,StrainResolution = True):

        '''
        Gets all leaf Taxa of the TargetTaxa list. Make sure the taxons in the TargetTaxa list aren't conflicting, i.e no parent/child relationships between them.
        '''
        
        for Taxon in TargetTaxa:
            TargetTaxon = Taxon

            if not isinstance(TargetTaxon,str):
                raise TypeError("Highest taxon node must be a string")
        
            #create the taxonomic graph part and get list of Taxids to fetch
            ncbi = NCBITaxa()

            HighestTaxid = ncbi.get_name_translator([TargetTaxon])[TargetTaxon][0]
            self.add_node(str(HighestTaxid),name = TargetTaxon, rank = ncbi.get_rank([HighestTaxid])[HighestTaxid],category = 'taxon')
            
            if StrainResolution:
                TaxidList = ncbi.get_descendant_taxa(TargetTaxon)
                TaxidNames= ncbi.translate_to_names(TaxidList)
                TaxidNodeTuples = tuple((str(TaxidList[i]),{'name':TaxidNames[i],'rank':ncbi.get_rank([TaxidList[i]])[TaxidList[i]], 'category':'taxon'}) for i in range(len(TaxidList)))
                self.add_nodes_from(TaxidNodeTuples)
            else:
                TaxidList = list(ncbi.get_descendant_taxa(TargetTaxon, collapse_subspecies=True))
                TaxidNames= ncbi.translate_to_names(TaxidList)
                TaxidNodeTuples = tuple((str(TaxidList[i]),{'name':TaxidNames[i],'rank':ncbi.get_rank([TaxidList[i]])[TaxidList[i]],'category':'taxon'}) for i in range(len(TaxidList)))
                self.add_nodes_from(TaxidNodeTuples)


            self.TaxidList = self.TaxidList+TaxidList
    
    def GetAllLeafTaxaFromTaxids(self,TaxidFile,StrainResolution = True):
        
        with open(TaxidFile) as Taxids:
            
            TargetTaxa = Taxids.read().splitlines()
            if 'no match' in TargetTaxa:
                TargetTaxa.remove('no match')
        
        for HighestTaxid in TargetTaxa:

            HighestTaxid = int(HighestTaxid)
            #create the taxonomic graph part and get list of Taxids to fetch
            ncbi = NCBITaxa()

    
            TargetTaxon = ncbi.get_taxid_translator([HighestTaxid])[HighestTaxid]
            self.add_node(str(HighestTaxid),name = TargetTaxon, rank = ncbi.get_rank([HighestTaxid])[HighestTaxid],category = 'taxon')
            
            if StrainResolution:
                TaxidList = ncbi.get_descendant_taxa(TargetTaxon)
                TaxidNames= ncbi.translate_to_names(TaxidList)
                TaxidNodeTuples = tuple((str(TaxidList[i]),{'name':TaxidNames[i],'rank':ncbi.get_rank([TaxidList[i]])[TaxidList[i]], 'category':'taxon'}) for i in range(len(TaxidList)))
                self.add_nodes_from(TaxidNodeTuples)
            else:
                TaxidList = list(ncbi.get_descendant_taxa(TargetTaxon, collapse_subspecies=True))
                TaxidNames= ncbi.translate_to_names(TaxidList)
                TaxidNodeTuples = tuple((str(TaxidList[i]),{'name':TaxidNames[i],'rank':ncbi.get_rank([TaxidList[i]])[TaxidList[i]],'category':'taxon'}) for i in range(len(TaxidList)))
                self.add_nodes_from(TaxidNodeTuples)
                
            self.TaxidList = self.TaxidList+TaxidList



    def FetchTaxonData(self, PeptideMapPath, sourceDB):
        '''
        gets the proteins corresponding to the target taxa from entrez and saves them
        in a json document.
        '''
    
        entrezDbName = 'protein'

        # options listed here https://www.ncbi.nlm.nih.gov/books/NBK49540/
        #sourceDBOptions = ['refseq[filter]','swissprot[filter]','protein_all[PROP]'] 
        
        saveLists = {}
        for Taxid in self.TaxidList:

            ncbiTaxId = str(Taxid) 

            # Find entries matching the query (only swissprot registered proteins for now)
            entrezQuery = sourceDB +' AND txid%s[ORGN]'%(ncbiTaxId) # AND 
            searchResultHandle = Entrez.esearch(db=entrezDbName, term=entrezQuery, retmax = 1000)
            searchResult = Entrez.read(searchResultHandle)
            searchResultHandle.close()

            #fetch corresponding proteins from NCBI entrez
            if searchResult['Count'] != '0':

                uidList = ','.join(searchResult['IdList'])
                proteinList = (Entrez.efetch(db=entrezDbName, id=uidList, rettype='fasta',retmax = 1000).read()).split('\n\n')
                proteinListOfLists = [ i.split('\n') for i in proteinList]
                remove = [list.pop(0) for list in proteinListOfLists]
                proteinList = [''.join(list) for list in proteinListOfLists]
                saveLists[Taxid] = proteinList[:-1] #remove last protein from list as it is empty
               

            
        with open(PeptideMapPath, 'w+') as savefile:        
            json.dump(saveLists,savefile)


    def CreateTaxonPeptideGraph(self,proteinListFile,minScore,minPeplength = 5,maxPeplength = 30,LCA=False):
        
        #read proteinlists from file that recorded the NCBI matches
        with open(proteinListFile) as file:
            ProteinTaxidDict = json.load(file)


        #loop digesting and adding the peptides to the graph, connecting to the corresponding taxid nodes

        for Taxid,AAsequences in ProteinTaxidDict.items():

            for AAsequence in AAsequences:

                # this needs to be replace
                cpdt = subprocess.check_output(['./cp-dt --sequence ' +AAsequence+ ' --peptides'], shell = True).decode('utf-8')
                peptideList = cpdt.split('\n')
                peptideList.pop(0)
                peptideList = [str[9:].split(': ') for str in peptideList]
                del peptideList[-3:]                                    #last three elements from cp-dt oupput are empty
                peptideList = [pep for pep in peptideList if minPeplength <= len(pep[0]) <= maxPeplength]
              
                PeptideNodes = tuple((pep[0],{'InitialBelief_0':1-float(pep[-1]),'InitialBelief_1': float(pep[-1]), 'category':'peptide'})  for pep in peptideList if float(pep[-1])>minScore)
                TaxonPeptideEdges = tuple((Taxid,pep[0]) for pep in PeptideNodes)
                
                #in this version, peptide nodes that already exist and are added again are ignored/ if they attributes differ, they are overwritten. 
                # conserves peptide graph structure, score will come from DB search engines anyways
                self.add_nodes_from(PeptideNodes)
                self.add_edges_from(TaxonPeptideEdges)


    def CreateTaxonPeptidegraphFromPSMresults(self,PSMResultsFile,proteinListFile,minScore,minPeplength = 5,maxPeplength = 30,LCA=False):

        #read proteinlists from file that recorded the NCBI matches
        with open(proteinListFile) as file:
            ProteinTaxidDict = json.load(file)

        Pepnames,Pepscores = loadSimplePepScore(PSMResultsFile)
        PepScoreDict = dict(zip(Pepnames,Pepscores))
        

        for Taxid,AAsequences in ProteinTaxidDict.items():

            self.add_node(Taxid, category='taxon')

            for AAsequence in AAsequences:

                cpdt = subprocess.check_output(['./../../bin/cp-dt --sequence ' +AAsequence+ ' --peptides'], shell = True).decode('utf-8')
                peptideList = cpdt.split('\n')
                peptideList.pop(0)
                peptideList = [str[9:].split(': ') for str in peptideList]
                del peptideList[-3:]                                    #last three elements from cp-dt oupput are empty
                peptideList = [pep for pep in peptideList if minPeplength <= len(pep[0]) <= maxPeplength]

                peptideList = [pep[0] for pep in peptideList]

                SelectedPeps = [pep for pep in Pepnames if pep in peptideList]
                PeptideNodes = tuple((pep,{'InitialBelief_0':1-PepScoreDict[pep]/100,'InitialBelief_1':PepScoreDict[pep]/100, 'category':'peptide'})  for pep in SelectedPeps if PepScoreDict[pep]>minScore)
                TaxonPeptideEdges = tuple((Taxid,pep[0]) for pep in PeptideNodes)
                   

                #in this version, peptide nodes that already exist and are added again are ignored/ if they attributes differ, they are overwritten. 
                # conserves peptide graph structure, score will come from DB search engines anyways
                self.add_nodes_from(PeptideNodes)
                self.add_edges_from(TaxonPeptideEdges)




class Factor:
    #represents noisy OR cpds, has dimension n(parensports)xn(peptide states(=2))
    def __init__(self,CPDarray,VariableArray):
        
        if isinstance(VariableArray, str):
            raise TypeError("VariableArray: Expected type list or array like, got string")

        Factor = namedtuple('Factor', ['array','arrayLabels'])
        self.Factor = Factor(CPDarray,VariableArray)
       
    
   
class Variable:
    #has dimension of petide states ergo 2
    def __init__(self,ProbabilityArray,VariableArray):
        self.Variable = namedtuple('Variable', ['array','arrayLabels'])
        Variable.array = ProbabilityArray
        Variable.arrayLabels = VariableArray

    


#the variable and factor types might be unecessary as i do not need a lot of flexibility in the input

#TODO implement checking& error raising if input aren't np.array/ list or array like  

class FactorGraph(nx.Graph):

    def __init__(self):
        super().__init__()        
        
    def ConstructFromProteinPeptideGraph(self,ProteinPeptideGraph):
        ''''
        Takes a graph of proteins to peptides as input and adds the noisy-OR factors
        '''
        nodelist = list(ProteinPeptideGraph.nodes(data=True))
        self.add_nodes_from(nodelist)
        for node in nodelist:                   
            #create noisy OR cpd per peptide
            if node[1]['category']=='peptide':                        
                degree = ProteinPeptideGraph.degree(node[0])                                  
                neighbors = list(ProteinPeptideGraph.neighbors(node[0]))
                self.add_node(node[0]+' CPD', category = 'factor', ParentNumber = degree)
                self.add_edges_from([(node[0]+' CPD',x) for x in neighbors])
                self.add_edge(node[0]+' CPD',node[0])
               
               
                
        return [self]

    def ConstructFromTaxonGraph(self,TaxonPeptideGraph):
        ''''
        Takes a graph of Taxa to peptides in networkx form as input and adds the factor nodes
        '''
        nodelist = list(TaxonPeptideGraph.nodes(data=True))
        self.add_nodes_from(nodelist)
        for node in nodelist:                   
            #create noisy OR cpd per peptide
            if node[1]['category']=='peptide':                        
                degree = TaxonPeptideGraph.degree(node[0])                                  
                neighbors = list(TaxonPeptideGraph.neighbors(node[0]))
                self.add_node(node[0]+' CPD', category = 'factor',ParentNumber = degree)
                self.add_edges_from([(node[0]+' CPD',x) for x in neighbors])
                self.add_edge(node[0]+' CPD',node[0])

#separate the connected components in the subgraph
def SeparateSubgraphs(graphIN,NodesToKeep):
    '''
    separations of subgraphs (create news graphs for each subgraph)
    
    '''
    newG = CTFactorGraph(nx.Graph())

    # add nodes to keep
    newG.add_nodes_from((n,graphIN.nodes[n]) for n in NodesToKeep)
    #add edges if they are present in original graph
    newG.add_edges_from((n,nbr,d) 
        for n, nbrs in graphIN.adj.items() if n in NodesToKeep
        for nbr, d in nbrs.items() if nbr in NodesToKeep) 


    return newG #ListOfFactorGraphs



class CTFactorGraph(FactorGraph):
    ''''
    This class is a networkx graph representing the full graphical model with all variables, CTrees, and Noisy-OR factors
    '''
    
    def __init__(self,GraphIn,GraphType = 'Taxons'):
        '''
        takes either a graph or a path to a graphML file as input
        
        '''
        super().__init__()

        GraphTypes = ['Proteins','Taxons']
        if GraphType not in GraphTypes:
            raise ValueError("Invalid Graphtype. Expected one of: %s" % GraphTypes)

        if GraphType == 'Taxons':
            self.category = 'taxon'
        elif GraphType == 'Protein':
            self.category = 'protein'
        
        if type(GraphIn) == str:
            GraphIn = nx.read_graphml(GraphIn)
          

        # need these to create a new instance of a CT fractorgraph and not overwrite the previous graph
        self.add_edges_from(GraphIn.edges(data =True))
        self.add_nodes_from(GraphIn.nodes(data=True))
        
    
    def AddCTNodes(self):
        '''
        When creating the CTGraph and not just reading from a previously saved graph format, use this command to add the CT nodes
        '''
        #create the convolution tree nodes and connect them in the graph
        ListOfEdgeAddList = []
        ListOfEdgeRemoveList = []
        ListOfProtLists = []
        ListOfCTs = []
        ListOfFactors =[]
        for node in self.nodes(data = True) :

        #go through all factors with degree>2 and get their protein lists, then generate their conv. trees
            if node[1]['category'] == 'factor' and self.degree[node[0]] >2 :
                
                ProtList =[]
              
                for neighbor in self.neighbors(node[0]):
                    neighbornode = self.nodes[neighbor]
                    if neighbornode['category'] == self.category :
                        
                        ProtList.append(neighbor)
                
                ListOfCTs.append([1])
                ListOfFactors.append(node[0])
                ListOfProtLists.append(ProtList)
                ListOfEdgeAddList.append([('CTree ' + ' '.join(str(ProtList)),x) for x in ProtList])
                ListOfEdgeRemoveList.append([(node[0],x) for x in ProtList ])
                     
                    
        #Fill all info into graph structure, should probably do this inside the loop before, so that i can initialize the messages
        for i in range(len(ListOfCTs)):
            self.add_node('CTree ' + ' '.join(str(ListOfProtLists[i])), ConvolutionTree = str(ListOfCTs[i][0]), category = 'Convolution Tree', NumberOfParents = len(ListOfProtLists[i]))
            self.add_edge('CTree ' + ' '.join(str(ListOfProtLists[i])),ListOfFactors[i], MessageLength = len(ListOfProtLists[i])+1)
            self.add_edges_from(ListOfEdgeAddList[i])
            self.remove_edges_from(ListOfEdgeRemoveList[i])
 

    
    def SaveToGraphML(self,Filename):
        nx.write_graphml(self,Filename)


    def ComputeNetworkAttributes(self):
        '''
        Computes nodes attributes using builtin networkx functions
        Returns degree centrality, closenesscentrality, betweennesscentrality and eigencentrality
        '''
        DegreeCentrality = dict(sorted(nx.degree_centrality(self).items(), key=lambda item: item[1]))
        Closenesscentrality = dict(sorted(nx.closeness_centrality(self).items(), key=lambda item: item[1]))
        BetweennessCentrality = dict(nx.betweenness_centrality(self).items(), key=lambda item: item[1])
        Eigencentrality = dict(nx.eigenvector_centrality(self).items(), key=lambda item: item[1])

        return DegreeCentrality,Closenesscentrality,BetweennessCentrality,Eigencentrality
    
    def FillInFactors(self,alpha, beta):
        '''fills in the noisy or Factors according to detection and error probabilities given'''

        for node in self.nodes(data=True):                   
            #create noisy OR cpd per peptide
            if node[1]['category']=='factor':                        
                # add noisyOR factors
                degree = node[1]['ParentNumber']
                #pre-define the CPD array and fill it with the noisyOR values
                cpdArray = np.full([2,degree+1],1-alpha)         
                ExponentArray = np.arange(0,degree+1)
                cpdArray[0,:] = np.power(cpdArray[0,:],ExponentArray)*(1-beta)
                cpdArray[1,:] = np.add(-cpdArray[0,:],1)
                cpdArray = np.transpose(normalize(cpdArray))

                FactorToAdd = Factor(cpdArray,['placeholder',[node[0]+'0',node[0]+'1']])
               
                #add factor & its edges to network as an extra node  
                nx.set_node_attributes(self,{node[0]:FactorToAdd},'InitialBelief')   
    
    def FillInPriors(self,prior):
        '''fills in the taxon priors according to the given prior'''

        for node in self.nodes(data=True):                   
            #create noisy OR cpd per peptide
            if node[1]['category']=='taxon': 
                nx.set_node_attributes(self,{node[0]:{'InitialBelief_0':1-prior,'InitialBelief_1':prior}}) 






def GenerateCTFactorGraphs(ListOfFactorGraphs,GraphType = 'Taxons'):

    if type(ListOfFactorGraphs) is not list: ListOfFactorGraphs = [ListOfFactorGraphs]
    for Graph in ListOfFactorGraphs:
        CTFactorgraph = CTFactorGraph (Graph,GraphType) #ListOfCTFactorGraphs.append(CTFactorGraph(Graph,GraphType))
    return CTFactorgraph


if __name__== '__main__':
    Taxongraph = TaxonGraph()
    #Taxongraph.GetAllLeafTaxa(['adenoviridae'])
    #Taxongraph.FetchTaxonData('peptidemapapth')
    Taxongraph.CreateTaxonPeptidegraphFromMzID('/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD005104_Herpessimplex_1/human_refseq_Default_PSM_Report.txt','/home/tholstei/repos/PepGM_all/PepGM/resources/SampleData/PXD005104_Herpessimplex_1/herpesviridae.json',0.001)
    #Taxongraph.CreateExample()
    Factorgraph = FactorGraph()
    Factorgraph.ConstructFromTaxonGraph(Taxongraph)

    #Factorgraphs = Factorgraph.SeparateSubgraphs()

    CTFactorgraphs = GenerateCTFactorGraphs(Factorgraph)
    CTFactorgraphs.SaveToGraphML('/home/tholstei/repos/PepGM_all/PepGM/graph.graphml')
