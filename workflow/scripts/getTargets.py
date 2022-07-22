import argparse
import linecache

import mmh3
import numba as nb
import numpy as np
import pandas as pd


def save(output_path, result):
    """
    Save in a \n separated text file.

    :param output_path: str, output path
    :param result: lst, items to be saved
    """
    with open(output_path, 'w') as f:
        f.write('\n'.join(result))


def preprocess_query(input_path, output_path):
    """
    Extract protein accessions from PeptideShaker output (PSM report).

    :param input_path:
    :param output_path:
    :return: df, dataframe with split accessions and their weights
    """
    raw = pd.read_csv(input_path, sep='\t', error_bad_lines=False, usecols=['Protein(s)'])
    raw.columns = ['accession']
    raw['weight'] = 1 / (raw.accession.str.count(',') + 1)
    raw['accession'] = raw.accession.apply(lambda x: x.split(','))
    raw['psmid'] = raw.index + 1
    df = raw.explode('accession', ignore_index=True)
    # extract accessions and save them
    query = df.accession.tolist()
    save(output_path, query)
    return df


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
    match = [0] * len(query)
    for i in nb.prange(len(database)):
        for j in nb.prange(len(query)):
            # if there is a match append the corresponding line from database
            if database[i] == query[j]:
                match[j] = i + 1
    return match


def lines2taxids(match, path):
    """
    Map matches to taxids.

    :param match: lst, index of matches
    :param path: str, path to taxid database
    :return lst, mapped taxids
    """
    taxids = []
    for idx in match:
        if idx != 0:
            taxids.append((linecache.getline(path, idx).strip()))
        else:
            taxids.append('no match')
    return taxids


def score(df, taxids, output_path):
    """
    Score taxids according to their confidence and select the ones which are top scoring.

    :param df: df, contains protein 2 taxid mapping and the respective weights
    :param taxids: lst, mapped taxids
    :param output_path: str, output path
    :param subset: int, top <subset> scoring taxids
    """
    df['taxid'] = taxids
    df_score = df.groupby('taxid')['weight'].sum().reset_index()
    df_score = df_score.sort_values(by=['weight'], ascending=False)
    threshold = df_score.weight.median()
    top_scoring_df = df_score.loc[df_score["weight"] >= threshold]
    top_scoring_df.to_csv(output_path[:-4] + '_weights.csv', index=False)
    top_scoring_taxids = top_scoring_df.taxid.tolist()
    save(output_path, top_scoring_taxids)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find ')
    parser.add_argument('-rq', '--path_to_raw_query', help='path to raw query')
    parser.add_argument('-q', '--path_to_query', help='path to processed query')
    parser.add_argument('-d', '--path_to_database', help='path to hashed database')
    parser.add_argument('-r', '--path_to_mapped_taxids', help='path to results')
    parser.add_argument('-t', '--path_to_taxids', nargs='*')
    args = parser.parse_args()

    # prepare
    df_accession = preprocess_query(args.path_to_raw_query, args.path_to_query)
    query_accessions = hash_query(args.path_to_query)
    database_accessions = np.load(args.path_to_database)
    # lookup
    lookup = query_database(database_accessions, query_accessions)
    taxids = lines2taxids(lookup, args.path_to_taxids[0])
    # score
    score(df_accession, taxids, args.path_to_mapped_taxids)
