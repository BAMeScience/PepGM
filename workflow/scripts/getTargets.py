import datetime
import linecache
from pathlib import Path

import mmh3
import os, psutil
import numba as nb
import numpy as np


# preliminaries
path_to_resources = Path('../../resources/')
path_to_sample = Path('../../resources/test/')


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
    tax_ids = []
    for idx in match:
        tax_ids.append((linecache.getline(path, idx)))
    return tax_ids


def save_results_to_txt(path, results):
    """
    Save results.

    :param path: str, where to save results
    :param results: lst, results from database query
    """
    f_out = open(path, 'w')
    for result in results:
        f_out.write(str(result))
    f_out.close()


if __name__ == "__main__":
    #TODO: ask T where samples are to be saved
    query_accessions = hash_query(path_to_sample / 'sample_n100_head.txt')

    database_accessions = np.load(path_to_resources / 'accessions_hashed.npy')
    lookup = query_database(database_accessions, query_accessions)

    taxids = lines_to_taxids(lookup)

    save_results_to_txt(path_to_resources / 'results_taxids.txt', taxids)
    print('Found targets.')
