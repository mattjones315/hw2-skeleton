from .utils import Atom, Residue, ActiveSite
from .cluster import create_similarity_matrix, cluster_by_partitioning, cluster_hierarchically, \
        convert_indices_to_active_sites
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

    print(active_sites)

    # iterate over all pairs of sample
    for a1 in active_sites:
        for a2 in active_sites:
            if a1.name == a2.name:
                continue

            # iterate over all pairs of clusterings to compare assignments of a1 and a2
            for c1 in clusterings_1:

                if a1.name in c1 and a2.name in c1:

                    same = False

                    for c2 in clusterings_2:

                        # a1 and a2 are assigned to same cluster in both sets
                        if a1.name in c2 and a2.name in c2:
                            same = True

                    if same:
                        tp += 1
                    else:
                        fp += 1


                if a1.name in c1 and a2.name not in c1:
                    diff = False

                    for c2 in clusterings_2:

                        if a1.name in c2 and not a2.name in c2:
                            diff = True

                    if diff:
                        tn += 1
                    else:
                        fn += 1

    return (tp + tn) / (tp + fp + tn + fn)

def benchmark_clusters(active_sites, alg = "P", sim_matrix=None, metric="RMSD", N=1000):

    if sim_matrix is None:

        sim_matrix = create_similarity_matrix(active_sites, metric)

    ssds = []

    print("Benchmarking " + metric + " similarity")
    for i in tqdm(range(N)):

        if alg == "P":

            assignments = cluster_by_partitioning(active_sites, sim_matrix)
            ssds.append(sum_of_distances(assignments, sim_matrix))

        if alg == "H":

            assignments = cluster_hierarchically(active_sites, sim_matrix)
            ssds.append(sum_of_distances(assignments, sim_matrix))

    return np.array(ssds)

def benchmark_rand(active_sites, alg="P", sim_matrix=None, metric="RMSD", N=100, step=10):

    if sim_matrix is None:
        sim_matrix = create_similarity_matrix(active_sites, metric)

    rs = {}

    for i in tqdm(range(step, N+step, step)):

        if alg == "P":
            cs = []

            for j in range(i):
                assignments = cluster_by_partitioning(active_sites, sim_matrix)
                clusters = convert_indices_to_active_sites(assignments, active_sites)
                cs.append(clusters)

            rs[i] = compute_avg_rand(cs, active_sites)

    return rs

def compute_avg_rand(all_clusters, active_sites):

    rand = []

    for i in range(len(all_clusters)):

        for j in range(len(all_clusters)):

            if i != j:

                rand.append(rand_index(all_clusters[i], all_clusters[j], active_sites))

    return np.mean(np.array(rand))


def rand_versus(active_sites, sim_matrix=None, metric="RMSD", N=30, step=5):

    if sim_matrix is None:
        sim_matrix = create_similarity_matrix(active_sites, metric)

    rs = {}

    for i in tqdm(range(step, N+step, step)):

        a1 = cluster_by_partitioning(active_sites, sim_matrix, K=i)
        a2 = cluster_hierarchically(active_sites, sim_matrix, K=i)

        c1 = convert_indices_to_active_sites(a1, active_sites)
        c2 = convert_indices_to_active_sites(a2, active_sites)

        rs[i] = rand_index(c1, c2, active_sites)

    return rs