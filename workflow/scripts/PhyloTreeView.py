import pandas as pd
from ete3 import NCBITaxa
from ete3 import TreeStyle, faces, AttrFace, CircleFace

'''
This script plots the PepGM results onto a phylogenetic tree using the ete3 tree visualization tools
'''


def CreatePhyloTreeView(resultsfile, host, out):
    '''
    Plots the PepGM results onto a phylogenetic tree using the ete3 tree visualization tools
    :params resulsfile: str,file with PepGM .csv results
    :params host: str,host of the virs to be remove from the phylogenetic tree view
    :params out: str, output path
    '''
    ncbi = NCBITaxa()
    Results = pd.read_csv(resultsfile, names=['ID', 'score', 'type'])

    # remove host taxon from visualization
    HostTaxid = ncbi.get_name_translator([host])[host][0]
    HostTaxidList = [str(i) for i in ncbi.get_descendant_taxa(HostTaxid)] + [str(HostTaxid)]

    # keep only taxa in dataframe
    TaxIDS = Results.loc[Results['type'] == 'taxon']
    # adjust types
    TaxIDS.loc[:, 'score'] = pd.to_numeric(TaxIDS['score'], downcast='float')
    # sort accoring to scores
    TaxIDS = TaxIDS.sort_values('score', ascending=False)
    TaxIDS.drop(TaxIDS[TaxIDS.ID.isin(HostTaxidList)].index, inplace=True)

    # put the 15 highest scoring taxids into minimal connecting phylotree
    phylotree = ncbi.get_topology(TaxIDS.ID.tolist()[:15])

    # set score as node weights
    for nodeName in phylotree.iter_leaves():
        nodeName.add_features(weight=TaxIDS.loc[TaxIDS['ID'] == nodeName.name, 'score'].iloc[0] * 30)
        # nodeName.spname = TaxaNameDict[int(nodeName.name)]

    def layout(node):
        # if node.is_leaf():
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

    phylotree.render(out, tree_style=ts, w=300, units="mm")
