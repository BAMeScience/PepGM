import argparse
import json

import pandas as pd


def init_argparser():
    """Init argument parser."""
    parser = argparse.ArgumentParser()

    parser.add_argument('--UnipeptResponseFile', type=str, required=True, help='path to Unipept response .json file')
    parser.add_argument('--NumberOfTaxa', type=int, required=True, help='number of taxa to include in the output')
    parser.add_argument('--out', type=str, required=True, help='path to csv out file')
    parser.add_argument('--UnipeptPeptides', type=str, required=True, help='path to Unipept response .json file')
    parser.add_argument('--ProteomeSizePerTaxID', type=str, required=True, help='path to proteome size per taxID file')

    args = parser.parse_args()

    return args


def GetProteinCount(proteins_per_taxon):
    """
    Convert tab-separated taxon-protein counts to a dictionary.
    Parameters
    ----------
    proteins_per_taxon: str,
        Path to file that contains the counts.

    """
    protein_counts_per_taxid = {}

    # this tsv was previously created
    with open(proteins_per_taxon, 'r') as f:
        for line in f:
            taxid, count = line.strip().split('\t')
            protein_counts_per_taxid[taxid] = int(count)

    return protein_counts_per_taxid


def WeightTaxa(UnipeptResponse, PeptScoreDict, MaxTax, ProteinsPerTaxon, chunks=True, N=1):
    """
    Weight inferred taxa based on their (1) degeneracy and (2) their proteome size.
    Parameters
    ----------
    UnipeptResponse: str
        Path to Unipept response json file
    PeptScoreDict: dict
        Dictionary that contains peptide to score & number of PSMs map
    MaxTax: int
        Maximum number of taxons to include in the graphical model
    ProteinsPerTaxon: str
        Path to the file that contains the size of the proteome per taxID (tab-separated)
    chunks: bool
        Allow memory-efficient reading of large json files
    N: int
        tbd

    Returns
    -------
    dataframe
        Top scoring taxa

    """
    with open(PeptScoreDict, 'r') as file:
        PeptScoreDictload = json.load(file)

    if chunks:
        with open(UnipeptResponse, 'r') as file:
            UnipeptDict = {"peptides": []}
            for line in file:
                try:
                    print()
                    UnipeptDict["peptides"].extend(json.loads(line)["peptides"])
                except:
                    # TODO: Pieter fixes internal server error
                    # in the meantime, we work with the incomplete mapping
                    # UnipeptDict["peptides"] = [json.loads(line)["peptides"]]
                    print()
                    continue

    else:
        with open(UnipeptResponse, 'r') as file:
            UnipeptDict = json.load(file)

    # Convert a JSON object into a Pandas DataFrame
    # record_path Parameter is used to specify the path to the nested list or dictionary that you want to normalize
    UnipeptFrame = pd.json_normalize(UnipeptDict, record_path=['peptides'])
    # Merge psm_score and number of psms
    UnipeptFrame = pd.concat([UnipeptFrame,
                              pd.json_normalize(UnipeptFrame['sequence'].map(PeptScoreDictload))], axis=1)
    # Score the degeneracy of a taxa, i.e.,
    # how conserved a peptide sequence is between taxa.
    # Divide the number of PSMs of a peptide by the number of taxa the peptide is associated with
    UnipeptFrame['weight'] = UnipeptFrame['psms'].div([len(element) for element in UnipeptFrame['taxa']])
    UnipeptFrame = UnipeptFrame.explode('taxa', ignore_index=True)

    # Sum up the weights of a taxon and sort by weight
    TaxIDWeights = UnipeptFrame.groupby('taxa')['weight'].sum().reset_index().sort_values(by=['weight'],
                                                                                          ascending=False)
    # Retrieve the proteome size per taxid as a dictionary
    # This file was previously prepared by filtering a generic accession 2 taxid mapping file
    # to swissprot (i.e., reviewed) proteins only

    ProteomeSize = GetProteinCount(ProteinsPerTaxon)
    TaxIDWeights["proteome_size"] = TaxIDWeights['taxa'].astype(str).map(ProteomeSize)

    # Since large proteomes tend to have more detectable peptides,
    # we adjust the weight by dividing by the size of the proteome i.e.,
    # the number of proteins that are associated with a taxon
    TaxIDWeights["scaled_weight"] = TaxIDWeights["weight"] / (TaxIDWeights["proteome_size"]) ** N

    # Filter to taxa with a weight greater than the median weight
    # However, if len > 50, take the top 50 taxa
    TopTaxa = TaxIDWeights.loc[TaxIDWeights["scaled_weight"] >= TaxIDWeights.scaled_weight.median()]

    if len(TopTaxa.taxa) < 50:
        return UnipeptFrame[UnipeptFrame['taxa'].isin(TopTaxa.taxa)]
    else:
        TopTaxaSorted = TopTaxa.sort_values(by="scaled_weight", ascending=False)
        return UnipeptFrame[UnipeptFrame['taxa'].isin(TopTaxaSorted.taxa[0:MaxTax])]


if __name__ == '__main__':
    args = init_argparser()
    DF = WeightTaxa(args.UnipeptResponseFile, args.UnipeptPeptides, args.NumberOfTaxa,
                    args.ProteomeSizePerTaxID)
    DF.to_csv(args.out)
