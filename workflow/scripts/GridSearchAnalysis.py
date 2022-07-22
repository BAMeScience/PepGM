import os
import re

import matplotlib
import numpy as np
import pandas as pd
from ete3 import NCBITaxa
from matplotlib import pyplot as plt
from scipy.stats import entropy


def ComputeMetric(resultsfolder, host, output, weightsfile):
    """
    Compute a "goodness" metric for the parameters used in the grid search.
    Generates a barplot of the metric and returns the best identified parameters.
    Return a list of the best parameter set in order [alpha,beta,gamma].

    :param resultsfolder: str, path to PepGM resultsfolder
    :param host: str, host to be excluded from the parameter checked taxa
    :param output: str, name of the output .png file
    :param weightsfile: str, path to the .csv file containing potential taxids and their weights

    """

    # predefine necessary lists
    Params = []
    Metrics = []
    SumProportions = []
    Entropies = []
    TaxDistances = []
    WeightCoeffs = []
    ncbi = NCBITaxa()
    matplotlib.use('Agg')

    # file with weights of taxids
    Weights = pd.read_csv(weightsfile)
    Maxweight = Weights.max()['weight']
    AllTaxidsToAdd = []
    # remove 'no match' line from weight DF
    Weights = Weights[Weights.taxid.isin(['no match']) == False]
    Weights.drop(Weights[Weights.taxid.isin(['no match'])].index, inplace=True)
    # add descendant taxa into weight dataframe
    for Taxid in Weights['taxid']:
        TaxidsToAdd = ncbi.get_descendant_taxa(Taxid)
        TaxidsToAdd = [[txd, float((Weights.loc[Weights['taxid'] == Taxid]['weight']).to_string(index=False))] for txd
                       in TaxidsToAdd]
        AllTaxidsToAdd.append(TaxidsToAdd[:])

    AllTaxidsToAdd = [txd_1 for txd in AllTaxidsToAdd for txd_1 in txd]
    AllWeights = pd.concat([Weights, pd.DataFrame(AllTaxidsToAdd, columns=['taxid', 'weight'])], axis=0)

    for folders in os.listdir(resultsfolder):
        if os.path.isdir(resultsfolder + '/' + folders):
            for file in os.listdir(resultsfolder + '/' + folders):

                if file.endswith('.csv'):
                    HostTaxid = ncbi.get_name_translator([host])[host][0]
                    HostTaxidList = [str(i) for i in ncbi.get_descendant_taxa(HostTaxid)] + [str(HostTaxid)]

                    Results = pd.read_csv(resultsfolder + '/' + folders + '/' + file, names=['ID', 'score', 'type'])
                    TaxIDS = Results.loc[Results['type'] == 'taxon']
                    TaxIDS.loc[:, 'score'] = pd.to_numeric(TaxIDS['score'], downcast='float')
                    TaxIDS = TaxIDS.sort_values('score', ascending=False)
                    TaxIDS = TaxIDS[TaxIDS.ID.isin(HostTaxidList) == False]
                    TaxIDS.drop(TaxIDS[TaxIDS.ID.isin(HostTaxidList)].index, inplace=True)
                    # compute the metric

                    # what's the weight of the highest scoring taxid?
                    Weight = AllWeights.loc[AllWeights['taxid'] == int(TaxIDS.ID.head(1).item())]['weight'].head(
                        1).item()
                    WeightCoeff = Weight / Maxweight
                    WeightCoeffs.append(WeightCoeff)

                    # compare the posterior probbilities
                    FirstSum = np.sum(TaxIDS[:3]['score'])
                    NextSum = np.sum(TaxIDS[3:8]['score'])
                    SumProportion = FirstSum / NextSum
                    SumProportions.append(SumProportion)
                    Entropy = entropy(TaxIDS[:10]['score'])
                    Entropies.append(Entropy)

                    # compute the pairwise taxonomic distance of the first 4 Taxa
                    FourTaxa = np.array(TaxIDS[:4]['ID'])
                    tree = ncbi.get_topology(FourTaxa)
                    DistanceSum = (tree.get_distance(FourTaxa[0], FourTaxa[1])) ** 2 + tree.get_distance(FourTaxa[1],
                                                                                                         FourTaxa[
                                                                                                             2]) + tree.get_distance(
                        FourTaxa[2], FourTaxa[3])
                    TaxDistances.append(DistanceSum)
                    Matching = (1 / (Entropy)) * SumProportion * (1 / DistanceSum) * WeightCoeff

                    Metrics.append(Matching)

                    # get the corresponding parameters
                    params = re.findall('\d+\.\d+', file)
                    Params.append([params[0], params[1], params[2]])

    zipLists = zip(Metrics, Params)
    SortedPairs = sorted(zipLists, reverse=True)

    tuples = zip(*SortedPairs)
    Metrics, Params = [list(tuple) for tuple in tuples]

    figure1 = plt.figure(figsize=(30, 15), tight_layout=True)

    plt.barh(list(range(20)), Metrics[0:20], color='mediumvioletred')
    plt.xticks(fontsize=35)
    plt.yticks(list(range(20)), Params[0:20], fontsize=35)

    print(output)
    plt.savefig(output)
    plt.close()

    return Params[0]
