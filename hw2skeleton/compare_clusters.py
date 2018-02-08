from .utils import Atom, Residue, ActiveSite
from .similarity import SimilarityComparator
import numpy as np
from tqdm import tqdm

def sum_of_distances(clusters, sim_matrix):
    """

    :param clusters: list, list of ActiveSite objects where each item in the list is a clustering
    :param sim_matrix: np.array NxN, similarity matrix for N ActiveSites
    :return: float, sum of distances within a cluster
    """

    Q = 0.0
    for c in clusters:

        for i in range(len(c)):
            for j in range(len(c)):
                Q += sim_matrix[i, j]

    return Q / len(clusters)


def rand_index(clusterings_1, clusterings_2, active_sites):
    """

    :param clusterings_1: list, list of ActiveSite objects where each item in the list is a clustering
    :param clusterings_2: list, list of ActiveSite objects where each item in the list is a clustering
    :param active_sites: list, list of ActiveSite objects that were clustered in both clusterings_1 and *_2
    :return: float, rand index used as a metric for comparing similarities of clustering assignments
    """
    tp = 0
    tn = 0
    fp = 0
    fn = 0

    # iterate over all pairs of sample
    for a1 in active_sites:
        for a2 in active_sites:
            if a1.name == a2.name:
                continue

            # iterate over all pairs of clusterings to compare assignments of a1 and a2
            for c1 in clusterings_1:
                for c2 in clusterings_2:

                    if a1.name in c1 and a2.name in c1:
                        # a1 and a2 are assigned to same cluster in both sets
                        if a1.name in c2 and a2.name in c2:
                            tp += 1

                        # a1 and a2 are assigned to different clusters in both sets
                        else:
                            fp += 1

                    if a1.name in c1 and a2.name not in c1:
                        # a1 and a2 are assigned to different clusters in both sets
                        if a1.name in c2 and not a2.name in c2:
                            tn += 1

                        # a1 and a2 are assigned to the same cluster in one set, but different
                        # clusters in another set
                        if a1.name in c2 and a2.name in c2:
                            fn += 1

    return (tp + tn) / (tp + fp + tn + fn)

