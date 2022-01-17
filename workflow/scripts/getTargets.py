import mmh3
import datetime
import numpy as np
import csv
import numba as nb
from numba.typed import List
import argparse
import linecache


parser = argparse.ArgumentParser(description='get a list of all taxa which are present in peptideShaker output')
parser.add_argument('-q', '--query',
                    help ='path to input list of protein accession')
parser.add_argument('-d', '--database',
                    help='path to input protein accession database')
parser.add_argument('-r', '--runtime', nargs='?',
                    help='decide if runtime analytics printed')
parser.add_argument('-t', '--taxids',
                    help='path to taxid database')
parser.add_argument('-sd', '--save_database', nargs='?',
                    help='path to directory where to save database in .npy format')
parser.add_argument('-sr', '--save_results', nargs='?',
                    help='path to directory where to save results')
args = parser.parse_args()


def hashDatabase(path, save=True, print_runtime=args.runtime, seed=18):
    """
    Perform mmh3 hashing on input database.
    MurmurHash3 is a non-cryptographic hashing algorithm.
    For more information on mmh3 algorithm visit https://pypi.org/project/mmh3/.

    :param path: str, input path to protein accession database (accessions need to be separated by \n)
    :param save: bool, choose to save database in numpy format
    :param print_runtime: bool, choose to print runtime analytics
    :param seed: int, hashing seed: select the same integer for identical hashing
    :return: -
    """
    print('Database hashing in process...')
    print('Preparing...')
    # initialize database
    database = np.empty([sum(1 for _ in open(path))])
    # save start time for runtime analytics
    start = datetime.datetime.now()
    print('Hashing...')
    # hash protein accessions - this is where the magic happens
    with open(path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            database[line_num-1] = mmh3.hash64(line, signed=False, seed=seed)[0]
    print('Hashing done.')
    if print_runtime:
        print('Time for database hashing: ', datetime.datetime.now() - start)
    print('---')
    if save:
        np.save(args.save_database, database)
    return database


@nb.jit(nopython=True, parallel=True)
def queryDatabase(database, query):
    """
    Loop to intersect.

    :param database: lst, hashed protein accession numbers
    :param query: lst, query accession numbers
    :return: lst, matching lines of database
    """
    lst = [0] * len(query)
    for i in nb.prange(len(database)):
        for j in nb.prange(len(query)):
            # if there is a match append the corresponding line from database
            if database[i] == query[j]:
                lst[j] = i + 1
    return lst


def hashQuery(path):
    """
    Hash query accession numbers.

    :param path: str, path to query accession numbers
    :return: lst, hashed query
    """
    accessions = []
    with open(path, 'r') as f:
        for line in f:
            line.strip()
            accessions.append(line)
    query = np.array([mmh3.hash64(accession, signed=False, seed=18)[0] for accession in accessions])
    return query


def saveResults(results):
    """
    Save taxon IDs in a text file.
    :param results: lst, taxids
    """
    output = open(args.save_results, 'w')
    for element in results:
        output.write(str(element))
    output.close()


def linesToTargets(match, path):
    """
    Convert index of match to taxon IDs.
    :param match: lst, index of matches
    :param path: str, path to taxid database
    :return lst, taxids
    """
    taxids = []
    for idx in match:
        taxids.append((linecache.getline(path, idx)))
    return taxids


def getTargets(create_database=False):
    """
    Main.
    :param create_database: bool, chose whether to build new or use existing database
    """

    if create_database:
        database = hashDatabase(args.database)
    else:
        database = np.load(args.database)
    print('Length of database: ', len(database))


    print('Query hashing...')
    query_accessions = hashQuery(args.query)
    print('Query hashing done.')
    print('---')


    print('Looking up...')
    start = datetime.datetime.now()
    match = queryDatabase(database, query_accessions)
    print('Lookup done.')
    print('Time for lookup: ', datetime.datetime.now() - start)
    print('---')


    start = datetime.datetime.now()
    taxids = linesToTargets(match, args.taxids)
    print('Lookup shape: ', len(match))
    print('Time for lookup: ', datetime.datetime.now() - start)
    print('---')


    print('Writing into output file...')
    saveResults(taxids)
    print('Done.')

getTargets()
