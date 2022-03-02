from ete3 import PhyloTree,Tree,TreeStyle,NodeStyle,faces,AttrFace,CircleFace
from ete3 import NCBITaxa
import pandas as pd
import csv


'''
This script plots the PepGM results onto a phylogenetic tree using the ete3 tree visualization tools
'''

#resultsfile = '/home/tholstei/repos/PepGM_all/PepGM/results/refseqViral/PXD025131_Sars_CoV_2/Prior0.1/chlorocebus_refseqViral_PepGM_Results_a0.01_b0.05_p0.1.csv'
#out = '/home/tholstei/repos/PepGM_all/PepGM/results/refseqViral/PXD025131_Sars_CoV_2/PhyloTree.png'
#host = 'chlorocebus'

#resultsfile = '/home/tholstei/repos/PepGM_all/PepGM/results/refseqViral/PXD002936_avian_bronchitis/Prior0.3/chicken_refseqViral_PepGM_Results_a0.1_b0.1_p0.3.csv'
#out = '/home/tholstei/repos/PepGM_all/PepGM/results/refseqViral/PXD002936_avian_bronchitis/Phylotree.png'
#host = "gallus gallus"

#resultsfile = '/home/tholstei/repos/PepGM_all/PepGM/results/refseqViral/PXD005104_Herpessimplex_1/Prior0.1/human_refseqViral_PepGM_Results_a0.01_b0.4_p0.1.csv'
#out = '/home/tholstei/repos/PepGM_all/PepGM/results/refseqViral/PXD005104_Herpessimplex_1/Phylotree.png'
#host = "homo sapiens"

#resultsfile = '/home/tholstei/repos/PepGM_all/PepGM/results/refseqViral/PXD025130_Sars_CoV_2/Prior0.1/chlorocebus_refseqViral_PepGM_Results_a0.01_b0.05_p0.1.csv'
#out = '/home/tholstei/repos/PepGM_all/PepGM/results/refseqViral/PXD025130_Sars_CoV_2/PhyloTree.png'
#host = 'chlorocebus'

resultsfile = '/home/tholstei/repos/PepGM_all/PepGM/results/refseqViral+cowpoxstrains/PXD003013_Cowpox_BR/Prior0.1/human_refseqViral+cowpoxstrains_PepGM_Results_a0.2_b0.01_p0.1.csv'
out = '/home/tholstei/repos/PepGM_all/PepGM/results/refseqViral+cowpoxstrains/PXD003013_Cowpox_BR/PhyloTree.png'
host = 'homo sapiens'

ncbi = NCBITaxa()

#get csv results file
Results = pd.read_csv(resultsfile, names = ['ID','score','type'])
#remove host taxon from visualization
HostTaxid = ncbi.get_name_translator([host])[host][0]
HostTaxidList = [str(i) for i in ncbi.get_descendant_taxa(HostTaxid)]+[str(HostTaxid)]
        
#keep only taxa in dataframe
TaxIDS = Results.loc[Results['type']=='taxon']
#adjust types
TaxIDS.loc[:,'score'] = pd.to_numeric(TaxIDS['score'],downcast = 'float')
#sort accoring to scores
TaxIDS = TaxIDS.sort_values('score', ascending = False)
#drop host taxa /descendants of host taxa
#TaxIDS = TaxIDS[TaxIDS.ID.isin(HostTaxidList)==False]
TaxIDS.drop(TaxIDS[TaxIDS.ID.isin(HostTaxidList)].index, inplace=True)

#put the 15 highest scoring taxids into minimal connecting phylotree
phylotree = ncbi.get_topology(TaxIDS.ID.tolist()[:15])
#print(phylotree.get_ascii(attributes=["sci_name"]))
#phylotree = PhyloTree(phylotree.write())
#print(phylotree.traverse())
#phylotree.annotate_ncbi_taxa()

#get a dictionary to translate TaxiDs into Taxa names
#TaxaNameDict = ncbi.get_taxid_translator(phylotree.get_species())
#TaxaNames = [TaxaNameDict[int(tax)] for tax in TaxIDS.ID.tolist()[:15]]

#write the phylotree into a tree to be able to manipulate nodesizes



#set score as node weights
for nodeName in phylotree.iter_leaves():
    nodeName.add_features(weight=TaxIDS.loc[TaxIDS['ID']== nodeName.name,'score'].iloc[0] *30)
    #nodeName.spname = TaxaNameDict[int(nodeName.name)]

#annote nodes with psecies names
#for node in phylotree.traverse():
    #print(node.name)
    #print(node.sci_name)
    #node.spname = ncbi.get_name_translator([node.name])[int(node.name)]



def layout(node):
    #if node.is_leaf():
        # Add node name to laef nodes
    N = AttrFace("sci_name", fsize=14, fgcolor="black")
    faces.add_face_to_node(N, node, 0)
    if "weight" in node.features:
        # Creates a sphere face whose size is proportional to node's
        # feature "weight"
        C = CircleFace(radius=node.weight, color="RoyalBlue", style="sphere")
        # Let's make the sphere transparent
        C.opacity = 1
        # And place as a float face over the tree
        faces.add_face_to_node(C, node, 0, position="float")

ts = TreeStyle()
ts.layout_fn = layout
#ts.show_leaf_name = False

phylotree.render(out,tree_style = ts, w=300, units="mm")   




