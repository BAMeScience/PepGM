""" This module formats grid search results and generates a phylo tree."""
import argparse

from GridSearchAnalysis import *
from PhyloTreeView import *


def SaveReducedCSV(results_gridsearch, host, output):
    """
    Save best results.
    :param results_gridsearch: str, path of cvs to be reduced
    :param host: str, if included in the search results, the host will be excluded from the results list
    :param output: str, output path

    """
    host_taxid = ncbi.get_name_translator([host])[host][0]
    host_taxid_list = [str(i) for i in ncbi.get_descendant_taxa(host_taxid)] + [str(host_taxid)]

    results_df = pd.read_csv(results_gridsearch, names=['ID', 'score', 'type'])
    # filter tax ids
    tax_ids = results_df.loc[results_df['type'] == 'taxon']
    # convert to numeric
    tax_ids.loc[:, 'score'] = pd.to_numeric(tax_ids['score'], downcast='float')
    tax_ids = tax_ids.sort_values('score', ascending=False)
    # drop all non host taxids
    tax_ids = tax_ids[tax_ids.ID.isin(host_taxid_list) == False]
    tax_ids.drop(tax_ids[tax_ids.ID.isin(host_taxid_list)].index, inplace=True)
    tax_ids.to_csv(output, index=False)


def MoveBestResultsPlot(filepath, out):
    """
    Copy bar plot to results folder.
    :param filepath: list, the three grid search parameters identified as best for the sample at hand
    :param out: output path where the file will be copied

    """
    os.system('cp ' + filepath + ' ' + out)


if __name__ == "__main__":
    # init argparser
    parser = argparse.ArgumentParser(description='Downstream analysis of grid search results.')
    required = parser.add_argument_group('required arguments')
    required.add_argument('--results', required=True, help='path to folder with raw results from belief propagation (.csv)')
    required.add_argument('--out', required=True, help='output path to save parameter grid search results')
    required.add_argument('--host', required=True, help='name of the host: to be excluded from parameter checked taxa')
    required.add_argument('--refdb', required=True, help='name of the reference database')
    required.add_argument('--weights', required=True, help='path to file with weighted taxids')
    args = parser.parse_args()

    # init NCBI API
    ncbi = NCBITaxa()

    # score grid search results with empirical metric
    parameter = ComputeMetric(args.results, args.host, args.out, args.weights)
    results = args.results + '/Prior' + str(parameter[2]) + '/' + args.refdb + '_PepGM_Results_a' + str(
        parameter[0]) + '_b' + str(parameter[1]) + '_p' + str(parameter[2])

    # save the reduced results as csv
    SaveReducedCSV(results + '.csv', args.host, args.results + '/PepGm_Results.csv')
    MoveBestResultsPlot(results + '.png', args.results + '/PepGM_ResultsPlot.png')

    # save a phylogenetic tree view of grid search results
    CreatePhyloTreeView(results + '.csv', args.host, args.results + '/PhyloTreeView.png')
