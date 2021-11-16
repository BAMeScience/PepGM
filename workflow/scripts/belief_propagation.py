#implementation of belief propagation on a peptide-protein graph
#__________________________________________________________________________________________
import numpy as np
import networkx as nx
import pandas as pd
import math
from scipy.signal import fftconvolve
import pandas as pd
from Bio import Entrez

import time

#TODO add methods that check the network/Factor graph


#e-mail for fetching viral protein sequences from entrez 
Entrez.email = 'tanja.holstein@bam.de'
Entrez.tool = 'PepGM'

#noisyOR methods
pRandomEmission = 0.4
pDetection = 0.8
pProteinprior = 0.5
Taxonprior = 0.5

#peptides from cp-dt need to have minimum lenngth <peplength> to be included in the graph
peplength = 3

#messagelog = [] #list of triples with from (node, to node, message) that logs all messages sent before the message class is initialized
def normalize(Array):
    normalizedArray = Array/np.sum(Array)
    return normalizedArray



#implementation of the convolution tree according to serang    
# not written by me!!          
class CTNode():

    def __init__(self, jointAbove):
        # normalize for greater precision
        self.jointAbove = normalize(jointAbove)

        self.leftParent = None
        self.rightParent = None

        self.likelihoodBelow = None
    
    #passing msges down: adding variables
    @classmethod
    def createCountNode(cls, lhs, rhs):
        #create a cound node with joint prob for two parents above(vs the init if we have no parents)
        jointAbove = fftconvolve(lhs.jointAbove, rhs.jointAbove)
        result = cls(jointAbove)

        result.leftParent = lhs
        result.rightParent = rhs

        return result


    #passing messages up : subtracting variables
    def messageUp(self, answerSize, otherJointVector):
        startingPoint = len(otherJointVector)-1
        resultcheck = fftconvolve(otherJointVector[::-1], self.likelihoodBelow)
        result = fftconvolve(otherJointVector[::-1], self.likelihoodBelow)[startingPoint:startingPoint+answerSize]

        return normalize(result)

    def messageUpLeft(self):
        return self.messageUp(len(self.leftParent.jointAbove[0]), self.rightParent.jointAbove[0])
    def messageUpRight(self):
        return self.messageUp(len(self.rightParent.jointAbove[0]), self.leftParent.jointAbove[0])

    # once all messages are received
    def posterior(self):
        return normalize(self.jointAbove * self.likelihoodBelow)
    
    def MessagesUp(self):
        check1 = self.jointAbove
        check2 = self.likelihoodBelow
        return self.likelihoodBelow


class ConvolutionTree:
    def __init__(self, nToSharedLikelihoods, proteins):
        self.nToSharedLikelihoods = nToSharedLikelihoods
        self.logLength = int(math.ceil(np.log2(float(len(proteins))))) #length we need
        self.allLayers = []
        self.buildFirstLayer(proteins)
        self.buildRemainingLayers()
        self.propagateBackward()
        self.nProteins = len(proteins)

    def buildFirstLayer(self, proteins):
        # construct first layer (of proteins)
        layer = []
        for prot in proteins:
            protNode = CTNode(prot)
            layer.append( protNode )

        # pad with necessarily absent dummy variables so that the
        # number of variables is a power of 2; this is not the most
        # efficient method for this. because they are absent, they won't influence the
        # total sum, and thus Ds.
        for i in range(0, 2**self.logLength - len(proteins) ):
            # this protein cannot be present, therefor set propbaility array to (0,1)
            layer.append( CTNode([np.array([1,0])]) ) #TODO change this order

        self.allLayers.append(layer)

    def buildRemainingLayers(self):
        # construct layers of count nodes
        for L in range(self.logLength):
            #print('layers needed: ',int(len(self.allLayers[0])/(2**(L+1))))
            mostRecentLayer = self.allLayers[-1]
            layer = []
            for i in range(int(len(self.allLayers[0]) / (2**(L+1)))):
                leftParent = mostRecentLayer[i*2]
                rightParent = mostRecentLayer[i*2+1]
                countNode = CTNode.createCountNode(leftParent, rightParent)
                layer.append( countNode )

            # add connection to remaining nodes (when layer above is not a power of 2)
            self.allLayers.append(layer)

        # final node gets (Ds | N) multiplied into its likelihoodBelow
        finalNode = self.allLayers[-1][0]
        # normalize for greater precision
        finalNode.likelihoodBelow = normalize(self.nToSharedLikelihoods)
        self.LastNode = finalNode
    
    def propagateBackward(self):
      # propagate backward, setting likelihoodBelow.
      # the loop has upper bound at logLength+1
      # because of the layer of proteins
        for L in range(1, self.logLength+1)[::-1]:
            layer = self.allLayers[L]

            for i in range(len(layer)):
                node = layer[i]
                leftParent = node.leftParent
                rightParent = node.rightParent

                leftParent.likelihoodBelow = node.messageUpLeft()
                rightParent.likelihoodBelow = node.messageUpRight()

        self.proteinLayer = self.allLayers[0]

    def posteriorForVariable(self, protInd):
        return self.proteinLayer[protInd].posterior()
    
    def MessageToVariable(self,protInd):
        return self.proteinLayer[protInd].MessagesUp()
        

    def MessageToSharedLikelihood(self):
        return self.LastNode.jointAbove[0][0:(self.nProteins+1)]

    
    

#class to hold and define all sum-product message passing methods
#TODO detect loops & use dampening when messages do not converge
#TODO clean up code & remove attributes I no langer use

class Messages():

    #class that holds the messages of itereation t and iteration t+1 as dictionaries

    def __init__(self,CTGraphIn):

        #if not isinstance(CTGraphIn,FactorGraphGeneration.CTFactorGraph):
            #raise TypeError("Input graph needs to be a CT FactorGraph)")

        self.Msg = {}
        self.MsgNew = {}
        self.MsgLog = {}
        self.Graph = CTGraphIn
        self.MaxVal = None
        self.FullResidual ={}
        self.FullResidualNew = {}
        self.InitialBeliefs = {}
        self.CurrentBeliefs = {}
        self.CurrentBeliefsNew = {}
        self.category = CTGraphIn.category
        #TODO check if I truly need all three of msg new, msglog and msg. chech if i need both fullresidual and fullresidual new, 
    
        
        for node in CTGraphIn.nodes(data = True):
            if node[1]['category'] == 'factor':
                self.InitialBeliefs[node[0]]= node[1]['InitialBelief'].Factor.array
            elif node[1]['category'] == 'peptide' or node[1]['category']==CTGraphIn.category:
                self.InitialBeliefs[node[0]]= np.array([node[1]['InitialBelief_0'],node[1]['InitialBelief_1']])
            else:
                self.InitialBeliefs[node[0]]= np.ones(4) #this entry will never be used as the convolution trees do not hold beliefs
        
        self.CurrentBeliefs = self.InitialBeliefs.copy()
        self.CurrentBeliefsNew = self.InitialBeliefs.copy()
        
       
        for node1,node2,data in CTGraphIn.edges(data = True):
            StartName, EndName = node1, node2

            if 'MessageLength' in data:
                self.Msg[(StartName, EndName)] = np.ones(data['MessageLength'])
            else:
                self.Msg[(StartName, EndName)] = np.array([0.5,0.5])
            
            self.Msg[(EndName, StartName)] = self.Msg[(StartName, EndName)]
            
            if 'MessageLength' in data:
                self.MsgNew[(StartName, EndName)] = np.zeros(data['MessageLength'])
            else:
                self.MsgNew[(StartName, EndName)] = np.array([0,0])
            
            self.MsgNew[(EndName, StartName)] = self.MsgNew[(StartName, EndName)]
        
        self.MsgLog = self.MsgNew.copy()
        
           

    #variables (peptides,proteins,taxa)  
    def  GetIncomingMessageVariable(self,Node,NodeIN):
        #Make sure to only multiply message in again if they have changed. Without checking this, peptide probs were multiplied in again and again
        #if (np.asarray([self.Msg[Node,NodeIN] != np.asarray(self.MsgLog[Node,NodeIN])])).all():
         #   check1 = self.Msg[Node,NodeIN] 
        #    check2 = self.MsgLog[Node,NodeIN]
        returnedMessage = self.Msg[Node,NodeIN]
            #self.MsgLog[Node,NodeIN] = returnedMessage
        return returnedMessage

        #else:
        #    return [1,1]

        

    def  ComputeOutMessageVariable(self,NodeOUT,NodeIN):
         IncomingMessages = []
         NodeBelief = self.CurrentBeliefs[NodeOUT]
         for NodeOUTneighbors in self.Graph.neighbors(NodeOUT):

             if NodeOUTneighbors != NodeIN:# and (np.asarray([self.Msg[NodeOUTneighbors,NodeOUT] != np.asarray(self.MsgLog[NodeOUTneighbors,NodeOUT])])).all():
                IncomingMessages.append(self.GetIncomingMessageVariable(NodeOUTneighbors, NodeOUT))
         
         if not IncomingMessages:
                check = any(NodeBelief == self.InitialBeliefs[NodeOUT])
                if any(NodeBelief == self.InitialBeliefs[NodeOUT]):
                    return NodeBelief
                else:
                    return self.Msg[NodeOUT,NodeIN]

         else:
             IncomingMessages = np.asarray(IncomingMessages).reshape(len(IncomingMessages),2)
             OutMessage = normalize(np.multiply(NodeBelief,[np.prod(IncomingMessages[:,0]),np.prod(IncomingMessages[:,1])]))
             #self.CurrentBeliefsNew[NodeOUT] = OutMessage
             return OutMessage

    
    #factors (Conditional probability tables), handles different dimension of output/input variables
    def  GetIncomingMessageFactor(self,Node,NodeIN):
        check1 = self.Msg[Node,NodeIN] 
        check2 = self.MsgLog[Node,NodeIN]
        #Make sure to only multiply message in again if they have changed. Without checking this, peptide probs were multiplied in again and again
        #if (np.asarray([self.Msg[Node,NodeIN] != np.asarray(self.MsgLog[Node,NodeIN])])).all():
        returnedMessage = self.Msg[Node,NodeIN]
            #self.MsgLog[Node,NodeIN] = returnedMessage
        return returnedMessage
        #else:
        #    return [1,1]

    def  ComputeOutMessageFactor(self,NodeOUT,NodeIN):
         IncomingMessages = []
         NodeBelief = self.CurrentBeliefs[NodeOUT]

         for NodeOUTneighbors in self.Graph.neighbors(NodeOUT):
             if NodeOUTneighbors != NodeIN:
                 if [self.GetIncomingMessageFactor(NodeOUTneighbors, NodeOUT)]:                      #only the messages that have changed get multiplied into the current belief again
                    IncomingMessages.append(self.GetIncomingMessageFactor(NodeOUTneighbors, NodeOUT))
         
         if self.Graph.nodes[NodeIN]['category'] == 'Convolution Tree':
                IncomingMessages.append([1.,1.]) #handles empty & messages with only one value
                IncomingMessages = np.asarray(IncomingMessages).reshape(len(IncomingMessages),2)
                OutMessages = normalize(np.multiply(NodeBelief,[np.prod(IncomingMessages[:,0]),np.prod(IncomingMessages[:,1])]))
                #self.CurrentBeliefsNew[NodeOUT] = OutMessages
            
                return np.add(OutMessages[:,0],OutMessages[:,1])    
         else:
                if  np.asarray(IncomingMessages[0]).shape[0] > 2:
                    IncomingMessages = np.asarray(IncomingMessages).reshape(IncomingMessages[0].shape[0],1)
                    OutMessages = normalize(NodeBelief*IncomingMessages)
                    #self.CurrentBeliefsNew[NodeOUT] = OutMessages
                    return [np.sum(OutMessages[0,:]),np.sum(OutMessages[1,:])] 
                else :
                    IncomingMessages.append([1.,1.])
                    IncomingMessages = np.asarray(IncomingMessages).reshape(len(IncomingMessages),2)
                    OutMessages = normalize(np.multiply(NodeBelief,[np.prod(IncomingMessages[:,0]),np.prod(IncomingMessages[:,1])])) 
                    #self.CurrentBeliefsNew[NodeOUT] = OutMessages
                    return [np.sum(OutMessages[0,:]),np.sum(OutMessages[1,:])] 


    #CTree, computes all out messages in one go
    def ComputeOutMessagesCTtree(self,Node):

        ProtProbList = []
        OldProtProbList = []
        sharedLikelihoods = np.ones(self.Graph.nodes[Node]['NumberOfParents']+1)
        peptides = []
        ProtList = []
        self.CurrentBeliefsNew[Node] = np.ones(4)

        for NodesIN in self.Graph.neighbors(Node):
            if 'CPD' not in str(NodesIN):

                if type(self.Msg[NodesIN,Node])==list:
                    ProtProbList.append([self.Msg[NodesIN,Node]])
                else:
                    ProtProbList.append([self.Msg[NodesIN,Node].tolist()])
            
                if type(self.MsgLog[NodesIN,Node])==list:
                    OldProtProbList.append([self.MsgLog[NodesIN,Node]])
                else:
                    OldProtProbList.append([self.MsgLog[NodesIN,Node].tolist()])     
                 
                ProtList.append(NodesIN)
            else:
                peptides.append(NodesIN)
                sharedLikelihoods = np.multiply(sharedLikelihoods,self.Msg[NodesIN,Node])
                OldSharedLikelihoods = np.multiply(sharedLikelihoods,self.MsgLog[NodesIN,Node])
        
        if all(OldSharedLikelihoods != sharedLikelihoods) and any([(ProtProbList[i][0]) != (OldProtProbList[i][0])for i in range(len(ProtProbList))]):
        #only update when the shared likelihhods or at least on of the protein messages has changed
            CT = ConvolutionTree(sharedLikelihoods,ProtProbList)

            for protein in range(len(ProtList)):
              self.MsgNew[Node,ProtList[protein]] = CT.MessageToVariable(protein)
    
            for pep in peptides:
             self.MsgNew[Node,pep] = CT.MessageToSharedLikelihood()
        
        else:
            for protein in range(len(ProtList)):
              self.MsgNew[Node,ProtList[protein]] = self.Msg[Node,ProtList[protein]]
    
            for pep in peptides:
             self.MsgNew[Node,pep] = self.Msg[Node,pep]




        

    #keeps track of which CTs have been update already in the current computeUpdate() loop
    def CTupdatecheck(self,CT):

        if CT in self.ListOfCTs:
            return False
        else:
            self.ListOfCTs.append(CT)
            return True
   

   #computes the residual between message new/message for a give edge(nodein/nodeOUT)
    def ComputeResidual(self,NodeIN,NodeOUT):
        Msg1 = self.MsgNew[NodeIN,NodeOUT]
        Msg2 = self.Msg[NodeIN,NodeOUT]
        if len(self.MsgNew[NodeIN,NodeOUT]) != len(self.Msg[NodeIN,NodeOUT]):
            Msg2 = [1]*len(self.MsgNew[NodeIN,NodeOUT])
        return np.sum(np.abs(np.subtract(Msg1,Msg2)))  

    #computes new message for a given edge (startname,endname) in the direction startname->endname
    def SingleEdgeDirectionUpdate(self,StartName,EndName):

        if self.Graph.nodes[StartName]['category'] == self.category or self.Graph.nodes[StartName]['category'] == 'peptide':
            self.MsgNew[StartName,EndName] = self.ComputeOutMessageVariable(StartName, EndName) 
            
        if self.Graph.nodes[StartName]['category'] == 'Convolution Tree':
            CTCheck = self.CTupdatecheck(StartName)

            if CTCheck:
                self.ComputeOutMessagesCTtree(StartName)
    
        if self.Graph.nodes[StartName]['category'] == 'factor':
            self.MsgNew[StartName,EndName] = self.ComputeOutMessageFactor(StartName,EndName)
     
                
    

    #compute messages for all edges
    
    def ComputeUpdate(self, localloops = False):

        self.ListOfCTs = [] #keeps track of which CT has already been active

        if not isinstance(localloops,bool):
            raise TypeError("localloops needs to be boolean")


        if localloops and self.MaxVal:

            for EndName in self.Graph.neighbors(self.MaxVal[1]):
                StartName = self.MaxVal[1]
                self.SingleEdgeDirectionUpdate(StartName, EndName)
            


        else:
            for edge in self.Graph.edges():
                #update all edges
                StartName, EndName = edge[0], edge[1]
                self.SingleEdgeDirectionUpdate(StartName, EndName)

                StartName, EndName = edge[1], edge[0]
                self.SingleEdgeDirectionUpdate(StartName, EndName)
            
        
        for edge in self.Graph.edges():
            #compute all residuals of the messages in this loop
            StartName, EndName = edge[1], edge[0]
            self.FullResidual[(StartName, EndName)] = self.ComputeResidual(StartName, EndName)

            StartName, EndName = edge[0], edge[1]
            self.FullResidual[(StartName, EndName)] = self.ComputeResidual(StartName, EndName)
                
                
            

        
    
    
    #send only the message with largest Residual
    def updateResidualMessage(self,Residual):

        self.MaxVal = max(Residual, key = Residual.get)
        self.Msg[self.MaxVal] = self.MsgNew[self.MaxVal]
        check = self.CurrentBeliefsNew[self.MaxVal[0]]
        #self.CurrentBeliefs[self.MaxVal[0]] = self.CurrentBeliefsNew[self.MaxVal[0]]
        return Residual[self.MaxVal]

    
    

    #run the loopy BP, returns number of iterations
    def LoopyLoop(self,maxLoops, tolerance,local = False):
        
        if not isinstance(local,bool):
            raise TypeError("localloops needs to be boolean")

        k = 0
        MaxResidual = 100

        while k < maxLoops and MaxResidual > tolerance:

            #first, do 5 loops where i update all messages
            while k < 5:
                start_t = time.time()
                self.ComputeUpdate()
                self.MsgLog.update(self.Msg)
                self.Msg.update(self.MsgNew)
                #self.CurrentBeliefs.update(self.CurrentBeliefsNew)
                k += 1
                end_t = time.time()
                print( "time per loop" , k, " ", end_t-start_t)

            
            #now start with the residual message passing
            start_t = time.time()
            self.ComputeUpdate(localloops = local)
            MaxResidual = self.updateResidualMessage(self.FullResidual)
            end_t = time.time()
            print( "time per loop" , k, " ", end_t-start_t, "residual max", MaxResidual)
        
            k += 1
        
        #when converged, multiply in all messages to each variable to get the posteriors
        for Variable in self.Graph.nodes():

            if self.Graph.nodes[Variable]['category'] == self.category or self.Graph.nodes[Variable]['category'] == 'peptide':

                IncomingMessages =[]
            
                for VariableNeighbors in self.Graph.neighbors(Variable):
                    IncomingMessages.append(self.GetIncomingMessageVariable(VariableNeighbors, Variable))
                
                IncomingMessages = np.asarray(IncomingMessages).reshape(len(IncomingMessages),2)
                VariableMarginal = normalize(np.multiply(self.InitialBeliefs[Variable],[np.prod(IncomingMessages[:,0]),np.prod(IncomingMessages[:,1])]))
                self.CurrentBeliefs[Variable] = VariableMarginal


        



    
    def DetectOscillations():
        pass




#calibration through message passing of all subgraphs in the List of factor graphs
def CalibrateAllSubgraphs(ListOfCTFactorGraphs, MaxIterations, Tolerance,local = False):

    if not isinstance(ListOfCTFactorGraphs,list):
        raise TypeError("ListOfFactorGraphs needs to be a list of graphs")    
    if not isinstance(local,bool):
        raise TypeError("localloops needs to be boolean")

    ResultsList = []
    ResultsDict = {}
    NodeDict = {}

    

    for Graph in ListOfCTFactorGraphs:

        if Graph.number_of_nodes()>2:

            NodeDict.update(dict(Graph.nodes( data = 'category')))
            InitializedMessageObject = Messages(Graph)
            InitializedMessageObject.LoopyLoop(MaxIterations,Tolerance,local)
            ResultsList.append(InitializedMessageObject.CurrentBeliefs)
            ResultsDict.update(InitializedMessageObject.CurrentBeliefs)




    return ResultsList, ResultsDict, NodeDict

#save the resulstsdictionary from CalibrateAllSubgraphs to a csv file
def SaveResultsToCsv(ResultsDict,NodeDict,NameString):

    if not isinstance(NameString,str):
       raise TypeError("AddNameString needs to a string with Info on your run")
    if not isinstance(ResultsDict,dict):
       raise TypeError("Resultsdict must be diciontary")
    if not isinstance(NodeDict,dict):
       raise TypeError("Resultsdict must be dictionary")
    
    FullResultsDict = {key:[ResultsDict[key][1],NodeDict[key]] for key in ResultsDict.keys()}
    pd.DataFrame.from_dict(data = FullResultsDict, orient = 'index').to_csv(NameString, header = False) #(datetime.now().strftime("%Y-%-m-%d-%H-%M-%S")+'-' +NameString + '-Results.csv', header = False) 



        
def VisualizeResults(NodeDict,ResultsDict,GraphType,Graph,**ResultsList,):

    GraphTypes = ['protein','taxon']
    if GraphType not in GraphTypes:
        raise ValueError("Invalid Graphtype. Expected one of: %s" % GraphTypes)
    
    TaxidNameDict = {}

    #Get translation of taxids to names:
    for nodes in Graph.nodes(data=True):
        if nodes[1]['category']== 'taxon':
            TaxidNameDict[nodes[0]]=nodes[1]['name']
            



    FullResultsDict = {key:[ResultsDict[key][1],NodeDict[key]] for key in ResultsDict.keys()}
    resultsFrame= pd.DataFrame.from_dict(data = FullResultsDict, orient = 'index',columns = ['score','category'])
    proteinResults = resultsFrame.loc[resultsFrame['category']==GraphType].sort_values(by='score')
    proteinResults['score']= pd.to_numeric(proteinResults['score'])
    peptideResults = resultsFrame.loc[resultsFrame['category']=='peptide'].sort_values(by='score')
    proteinResults['score']= pd.to_numeric(proteinResults['score']) 
    proteinResults.rename(index = TaxidNameDict)


    hist1 = px.bar(proteinResults,y = proteinResults.index,x = 'score',orientation = 'h')
    hist2 = px.bar(peptideResults,y = peptideResults.index,x = 'score',orientation = 'h')

    hist1.show()
   #TODO add sliders to make results more legible

   



    
    


        

                   

        





 





if __name__== '__main__':
 
 Taxongraph = TaxonGraph()
 #Taxongraph.GetAllLeafTaxa(['adenoviridae'])
 Taxongraph.FetchTaxonData('peptidemapapth')
 Taxongraph.CreateTaxonPeptidegraphFromMzID('/home/tholstei/repos/VirusGraph/Data/Searchresults/adeno_refseqviruses/adeno_refseq_allviruses_Default_PSM_Report.txt','TaxonGraph_adenoviridae.json',0.001)
 #Taxongraph.CreateExample()
 Factorgraph = FactorGraph()
 Factorgraph.ConstructFromTaxonGraph(Taxongraph,)


 #To save a gml of your graph, uncomment this
 #GmlCTFactorgraphs = CTFactorGraph(Factorgraph)
 #GmlCTFactorgraphs.SaveToGml('TaxonpeptideGraph_herpesviridae_pxd005104')


 Factorgraphs = Factorgraph.SeparateSubgraphs()

 CTFactorgraphs = GenerateCTFactorGraphs(Factorgraphs)

 Resultlist,Resultsdict,Nodetypes = CalibrateAllSubgraphs(CTFactorgraphs,10000,0.006)
 save = SaveResultsToCsv(Resultsdict,Nodetypes,'TaxonnomyTest_adeno_adenovridae_alpha0.4_beta0.9')
 #VisualizeResults(Nodetypes,Resultsdict,'taxon',Taxongraph)






    


