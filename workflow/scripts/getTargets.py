import argparse
import linecache
import os

import mmh3
import numba as nb
import numpy as np
import pandas as pd


# argparser preliminaries
parser = argparse.ArgumentParser(description = 'Find ')
parser.add_argument('-rq', '--path_to_raw_query', help='path to raw query')
parser.add_argument('-q', '--path_to_query', help='path to processed query')
parser.add_argument('-d', '--path_to_database', help='path to hashed database')
parser.add_argument('-r', '--path_to_mapped_taxids', help='path to results')
parser.add_argument('-t', '--path_to_taxids', nargs='*')
args = parser.parse_args()


def save_to_txt(path, results):
    """
    Save results.

    :param path: str, where to save results
    :param results: lst, results from database query
    """
    f = open(path, 'w')
    [f.write(str(result)) for result in results]
    f.close()


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
    # write in output file
    accessions_final = open(output_path, "w")
    [accessions_final.write(element + "\n") for element in proteins]
    accessions_final.close()
    # save value counts
    df_accessions = pd.DataFrame()
    df_accessions['accession'] = proteins
    print(df_accessions)
    df = df_accessions.value_counts().rename_axis('accession').to_frame('counts')
    df.to_csv('/home/fkistner/pepgm/results/PXD002936_avian_bronchitis/chicken_refseq_query_valuecounts.csv')
    return df_accessions

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


def lines_to_taxids(match, path, df):
    """
    Map matches to taxids.

    :param match: lst, index of matches
    :param path: str, path to taxid database
    :return lst, mapped taxids
    """
    taxids = []
    for idx in match:
        if idx != 0:
            taxids.append((linecache.getline(path, idx)))
        else:
            taxids.append('Nan')
    df['taxids'] = taxids
    df_counts = df.value_counts()
    df.to_csv('/home/fkistner/pepgm/results/PXD002936_avian_bronchitis/chicken_refseq_accession_taxids.csv')
    df_counts.to_csv('/home/fkistner/pepgm/results/PXD002936_avian_bronchitis/chicken_refseq_count_taxids.csv')


    """
    taxids_unique = []
    # remove duplicates
    for taxid in taxids:
        if taxid not in taxids_unique:
            taxids_unique.append(taxid)
    return taxids_unique
    """

# prepare
df_accessions = preprocess_query(args.path_to_raw_query, args.path_to_query)
query_accessions = hash_query(args.path_to_query)
database_accessions = np.load(args.path_to_database)
# lookup
lookup = query_database(database_accessions, query_accessions)
taxids = lines_to_taxids(lookup, args.path_to_taxids[0], df_accessions)
# save
# save_to_txt(args.path_to_mapped_taxids, taxids)
