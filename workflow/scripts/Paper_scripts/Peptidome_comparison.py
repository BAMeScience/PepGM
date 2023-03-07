import re

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sbn
from Bio import SeqIO
from ete3 import NCBITaxa

ncbi = NCBITaxa()


# scripts that compute the peptidome similarity between taxa using only the peptides that were included in the PepGm graph
def GetNhighestTaxa(resultscsv, N):
    IDs = pd.read_csv(resultscsv, dtype=str)
    return IDs.ID.to_list()[0:N]


def GetPeptidesperTaxon(Graphin, Taxa):
    graph = nx.read_graphml(Graphin)
    PeptidomeDict = {}
    for node in graph.nodes(data=True):
        if node[1]['category'] == 'taxon' and node[0] in Taxa:
            neighbors = graph.neighbors(node[0])
            PeptidomeDict.update({node[0]: [n[:-4] for n in neighbors]})

    return PeptidomeDict


def digest(peptide, n_missed_sites=2):
    """
    In-situ trypsin digestion of peptides. Creates all possible cleavage products with up to n missed cleavage sites.
    :param n_missed_sites: int, number of allowed missed cleavage sites
    :param peptide: str, input peptide sequence
    :return: lst, digested peptides
    """

    # in-situ trypsin digestion (without errors)
    # cut after every R and K except a P follows
    trypsin_pattern = re.compile(r'(?<=[RK])(?=[^P])')
    digested_peps = list(re.split(trypsin_pattern, peptide))

    # in-situ trypsin digestion (with up to n errors)
    # prepare index for slicing
    # note that the last element of a list slice is not contained in the slice
    cut = 2
    counter = 0
    # initialize limit
    N = len(digested_peps) - 1

    for n in range(0, n_missed_sites):
        # shorten the list to avoid indexing error
        for i in range(0, N - counter):
            missed_peptide = "".join(digested_peps[i:i + cut])
            # sequentially append to list
            digested_peps.append(missed_peptide)
        counter += 1
        cut += 1

    return digested_peps


def GetPeptidesPerTaxon(FastaIN):
    sequences = set()
    with open(FastaIN) as file:
        for sequence in SeqIO.parse(file, "fasta"):
            sequences.update(digest(str(sequence.seq)))

    return sequences


def ComputeTaxonPeptidomeSimilarity(PeptidomeDict, ListOfFastaFiles, ListOfTaxa):
    SimMatrixMax = []
    SimMatrixMin = []
    Taxa1 = []
    Taxa2 = ListOfTaxa
    for taxon1 in PeptidomeDict.keys():
        Taxa1.append(taxon1)
        SimMatrixMaxRow = []
        SimMatrixMinRow = []
        for sequences2 in ListOfFastaFiles:
            peptides1 = set(PeptidomeDict[taxon1])
            peptides2 = GetPeptidesPerTaxon(sequences2)
            shared = len(peptides1.intersection(peptides2))
            try:
                SimMax = shared / (len(peptides1))
            except:
                SimMax = 0
            try:
                SimMin = SimMin = shared / (len(peptides1))
            except:
                SimMin = 0

            SimMatrixMaxRow.append(SimMax)
            SimMatrixMinRow.append(SimMin)
        SimMatrixMax.append(SimMatrixMaxRow)
        SimMatrixMin.append(SimMatrixMinRow)

    return SimMatrixMax, SimMatrixMin, Taxa1, Taxa2


def ComputeDetectedPeptidomeSimilarity(PeptidomeDict):
    SimMatrixMax = []
    SimMatrixMin = []
    Taxa1 = []
    Taxa2 = []
    for taxon1 in PeptidomeDict.keys():
        Taxa1.append(taxon1)
        SimMatrixMaxRow = []
        SimMatrixMinRow = []
        for taxon2 in PeptidomeDict.keys():
            Taxa2.append(taxon2)
            peptides1 = set(PeptidomeDict[taxon1])
            peptides2 = set(PeptidomeDict[taxon2])
            shared = len(peptides1.intersection(peptides2))
            try:
                SimMax = shared / (max(len(peptides1), len(peptides2)))
            except:
                SimMax = 0
            try:
                SimMin = SimMin = shared / (min(len(peptides1), len(peptides2)))
            except:
                SimMin = 0

            SimMatrixMaxRow.append(SimMax)
            SimMatrixMinRow.append(SimMin)
        SimMatrixMax.append(SimMatrixMaxRow)
        SimMatrixMin.append(SimMatrixMinRow)

    return SimMatrixMax, SimMatrixMin, Taxa1, Taxa2


if __name__ == "main":
    # ListofTaxa = ['human adenovirus 2','human adenovirus 6', 'human adenovirus 5', 'simian adenovirus 34','human adenovirus 57']
    # ListofFastas = ['/home/tholstei/repos/PepGM_all/PepGM/resources/Strainfiles/HUman_adenovirus_2.fasta','/home/tholstei/repos/PepGM_all/PepGM/resources/Strainfiles/human_adenovirus_5.fasta','/home/tholstei/repos/PepGM_all/PepGM/resources/Strainfiles/human_adenovirus_6.fasta',
    # '/home/tholstei/repos/PepGM_all/PepGM/resources/Strainfiles/simian_adenovirus_34.fasta','/home/tholstei/repos/PepGM_all/PepGM/resources/Strainfiles/human_adenovirus_57.fasta']

    # ListofTaxa=['Beaudette CK', 'Beaudette US', 'Beaudette', 'strain M41']
    # ListofFastas = ['/home/tholstei/repos/PepGM_all/PepGM/resources/Strainfiles/beaudette_CK.fasta','/home/tholstei/repos/PepGM_all/PepGM/resources/Strainfiles/beaudette_US.fasta',
    #               '/home/tholstei/repos/PepGM_all/PepGM/resources/Strainfiles/avian_beaudette.fasta','/home/tholstei/repos/PepGM_all/PepGM/resources/Strainfiles/avian_M41.fasta']

    # ListofTaxa = ['Hendra virus horse', 'Nipah Henipavirus']
    # ListofFastas = ['/home/fkistner/pepgm/resources/Strains/hendra_horse.fasta',
    #                 '/home/fkistner/pepgm/resources/Strains/nipah_henipavirus.fasta']

    # ListofTaxa = ["Human adenovirus 2", "Human adenovirus 6", "Human adenovirus 5", "Simian Adenovirus 34", "Human Adenovirus 57"]
    # ListofFastas = ['/home/fkistner/pepgm/resources/Strains/HUman_adenovirus_2.fasta',
    #                 "/home/fkistner/pepgm/resources/Strains/human_adenovirus_6.fasta",
    #                 "/home/fkistner/pepgm/resources/Strains/human_adenovirus_5.fasta",
    #                 "/home/fkistner/pepgm/resources/Strains/simian_adenovirus_34.fasta",
    #                 "/home/fkistner/pepgm/resources/Strains/human_adenovirus_57.fasta"]

    # ListofTaxa = ["Human alphaherpesvirus 1 strain 17",
    #               "Human alphaherpesvirus 1 strain RH2",
    #               "Human alphaherpesvirus1 strain KOS"]
    # ListofFastas = ["/home/fkistner/pepgm/resources/Strains/herpes_17.fasta",
    #                 "/home/fkistner/pepgm/resources/Strains/herpes_RH2.fasta",
    #                 "/home/fkistner/pepgm/resources/Strains/herpes_KOS.fasta"]
    #
    # ListofTaxa = ["Cowpox virus (Brighton Red)", "Acanthocystis turfacea chlorella virus 1", "Brazilian marseillevirus"]
    # ListofFastas = ["/home/fkistner/pepgm/resources/Strains/cowpox_014913.fasta",
    #                 "/home/fkistner/pepgm/resources/Strains/acanthocystis.fasta",
    #                 "/home/fkistner/pepgm/resources/Strains/brazilian.fasta"]

    # ListofTaxa = ["SarS Coronavirus 2", "SarS Coronavirus Tor2", "Bat coronavirus BM48-31/BGR/2008"]
    # ListofFastas =  ["/home/fkistner/pepgm/resources/Strains/sars-cov-2.fasta",
    #                  "/home/fkistner/pepgm/resources/Strains/sars-cov-tor2.fasta",
    #                  "/home/fkistner/pepgm/resources/Strains/bat_coronavirus.fasta"]

    #
    # ListofTaxa = ["Cowpox virus (Brighton Red)", "Orthopoxvirus Abatino"]
    # ListofFastas = ["/home/fkistner/pepgm/resources/Strains/cowpox_003013.fasta",
    #                 "/home/fkistner/pepgm/resources/Strains/orthopox_abatino.fasta"]
    #

    # ListofTaxa = ["Avian infectious bronchitis virus (strain Beaudette US)",
    #               "Avian infectious bronchitis virus (strain Beaudette)",
    #               "Avian infectious bronchitis virus (strain Beaudette CK)",
    #               "Avian infectious bronchitis virus (strain Beaudette M41)"]
    #
    # ListofFastas = ["/home/fkistner/pepgm/resources/Strains/beaudette_US.fasta",
    #                 "/home/fkistner/pepgm/resources/Strains/avian_beaudette.fasta",
    #                 "/home/fkistner/pepgm/resources/Strains/beaudette_CK.fasta",
    #                 "/home/fkistner/pepgm/resources/Strains/avian_M41.fasta"]

    # graphpath = '/home/fkistner/pepgm/results/sars2/refseqViral_PepGM_graph.graphml'
    # outmax = '/home/fkistner/pepgm/results/plots/sars2_peptidome_sim_MAX.png'
    # taxonresults = '/home/fkistner/pepgm/results/sars2/PepGm_Results.csv'

    Taxa = GetNhighestTaxa(taxonresults, 15)
    PeptiDict = GetPeptidesperTaxon(graphpath, Taxa)
    SimMax, SimMin, Taxa1, Taxa2 = ComputeTaxonPeptidomeSimilarity(PeptiDict, ListofFastas, ListofTaxa)
    Taxid1 = ncbi.get_taxid_translator(Taxa1)
    TaxonList1 = [Taxid1[int(tax)] for tax in Taxa]

    # Heatmap
    sbn.set(rc={'figure.figsize': (13.7, 10.27)})
    matrix = np.triu(np.ones_like(SimMax))
    np.fill_diagonal(matrix, 0)
    ax1 = sbn.heatmap(SimMax, yticklabels=TaxonList1, xticklabels=ListofTaxa, cmap="YlGnBu",
                      annot=[[round(n, 2) for n in inner_list] for inner_list in SimMax])
    # Rotate x-axis labels
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(outmax, dpi=600)
    plt.close()

# ax1 = sbn.heatmap(SimMin, xticklabels = TaxonList, yticklabels = TaxonList,cmap="YlGnBu")
# plt.tight_layout()
# plt.savefig(outmin)
# plt.close()
