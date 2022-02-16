import datetime
import gzip
import subprocess
from pathlib import Path

import mmh3
import numpy as np
from requests import get


# preliminaries
path_to_resources = Path('../../resources/')
command_to_split = "awk '{print $1}' '../../resources/protacc2taxids_virus.txt' > '../../resources/accessions.txt'; " \
                   "awk '{print $2}' '../../resources/protacc2taxids_virus.txt' > '../../resources/taxids.txt'"


def split_to_accessions_and_taxids(command):
    """
    Bash wrapper for splitting two-column file into two files with a single column.

    """
    subprocess.run(command, shell=True)


def hash_database(path, seed=18):
    """
    Perform mmh3 hashing on input database.
    MurmurHash3 is a non-cryptographic hashing algorithm.
    For more information on mmh3 algorithm visit https://pypi.org/project/mmh3/.

    :param path: str, input path to protein accession database (accessions need to be separated by \n)
    :param seed: int, hashing seed: use identical integer for identical hashing results
    :return: database, hashed database
    """
    # initialize database
    database = np.empty([sum(1 for _ in open(path))])
    # save start time for runtime analytics
    # hash protein accessions - this is where the magic happens
    with open(path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            database[line_num-1] = mmh3.hash64(line, signed=False, seed=seed)[0]
    np.save(path_to_resources / 'accessions_hashed.npy', database)
    return database


if __name__ == "__main__":
    split_to_accessions_and_taxids(command_to_split)
    hash_database(path_to_resources / 'accessions.txt')
    print('Created database.')
