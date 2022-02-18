import argparse
import datetime
import linecache
import os

import mmh3
import numba as nb
import numpy as np
import pandas as pd
import psutil

# argparser preliminaries
parser = argparse.ArgumentParser(description = 'Find ')
parser.add_argument('-rq', '--path_to_raw_query', help='path to raw PSM report')
parser.add_argument('-q', '--path_to_query', help='path to processed query')
parser.add_argument('-d', '--path_to_database', help='path to hashed database')
parser.add_argument('-t', '--path_to_mapped_taxids', help='path to results')
args = parser.parse_args()


def save_results_to_txt(path, results):
    """
    Save results.

    :param path: str, where to save results
    :param results: lst, results from database query
    """
    output = open(path, 'w')
    [output.write(str(result)) for result in results]
    output.close()


def preprocess_query(input_path, output_path):
    """
    Extract protein accessions from PeptideShaker output (PSM report) and save them.
    :param input_path: str, input path to raw query
    :param output_path: str, output path
    :return: -
    """
    # error bad lines should be remove when development is done
    psm_report = pd.read_csv(input_path, sep = '\t', error_bad_lines=False)
    # extract accession numbers
    accessions_raw = psm_report['Protein(s)'].tolist()
    proteins = []
    # split into sublists
    for accession in accessions_raw:
        accession.split(',')
        proteins.append(accession)
    # merge sublists
    proteins = [j for i in [protein.split(',') for protein in proteins] for j in i]
    save_results_to_txt(output_path, proteins)


def hash_query(path):
    """
    Hash input query.

    :param path: str, path to query
    :return: lst, hashed query
    """
    accessions = []
    with open(path, 'r') as f:
        for line in f:
            line.strip()
            accessions.append(line)
    query = np.array([mmh3.hash64(accession, signed=False, seed=18)[0] for accession in accessions])
    return query


@nb.jit(nopython=True, parallel=True)
def query_database(database, query):
    """
    Loop to intersect.

    :param database: lst, database
    :param query: lst, query
    :return: lst, line numbers where query and database are matching
    """
    lst = [0] * len(query)
    for i in nb.prange(len(database)):
        for j in nb.prange(len(query)):
            # if there is a match append the corresponding line from database
            if database[i] == query[j]:
                lst[j] = i + 1
    return lst


def lines_to_taxids(match, path='../../resources/taxids.txt'):
    """
    Map matches to taxids.

    :param match: lst, index of matches
    :param path: str, path to taxid database
    :return lst, taxids
    """
    taxids = []
    for idx in match:
        taxids.append((linecache.getline(path, idx)))
    taxids_unique = []
    # remove duplicates
    for taxid in tax_ids:
        if taxid not in taxids_unique:
            taxids_unique.append(x)
    return taxids_unique


# '/home/fkistner/pepgm/resources/chicken_refseq_Default_PSM_Report.txt'
# '/home/fkistner/pepgm/results/chicken_refseq_query_accessions.txt'
preprocess_query(args.path_to_raw_query, args.path_to_query)

# '/home/fkistner/pepgm/results/chicken_refseq_query_accessions.txt'
query_accessions = hash_query(args.path_to_query)

# '/home/fkistner/pepgm/resources/accessions_hashed.npy'
database_accessions = np.load(args.path_to_database)
lookup = query_database(database_accessions, query_accessions)
taxids = lines_to_taxids(lookup)

# '/home/fkistner/pepgm/results/mapped_taxids.txt'
save_results_to_txt(args.path_to_mapped_taxids, taxids)
print('Found targets.')
