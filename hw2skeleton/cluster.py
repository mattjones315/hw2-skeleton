from .utils import Atom, Residue, ActiveSite
from .similarity import SimilarityComparator
import numpy as np
from tqdm import tqdm


def compute_similarity(site_a, site_b, comparator="RMSD"):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """

    # Fill in your code here!

    sim_comp = SimilarityComparator(method=comparator)
    similarity = sim_comp.compare(site_a, site_b)

    return similarity


def convert_indices_to_active_sites(clusters, active_sites):

    n_clusters = [[] for i in range(len(clusters))]

    for c in range(len(clusters)):
        clust = clusters[c]
        for elem in clust:
            n_clusters[c].append(active_sites[elem].name)

    return n_clusters

def cluster_by_partitioning(active_sites, sim_matrix=None, K=5, MAX_ITER=1000):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """

    def check_prev(c1, c2):

        if c1 is None or c2 is None:
            return False

        for k in range(len(c1)):
            if set(c1[k]) != set(c2[k]):
                return False

        return True

    def find_centroids(clusters, sim_matrix, K):

        k_ii = np.zeros((K), dtype=int)

        for k in range(K):

            clust = clusters[k]
            if len(clust) == 0:
                continue

            curr = clust[0]
            curr_avg = np.mean(sim_matrix[curr, clust])
            for elem in clust:
                avg = np.mean(sim_matrix[elem, clust])
                if avg < curr_avg:
                    curr = elem
                    curr_avg = avg

            k_ii[k] = int(curr)

        return k_ii

    if active_sites is None:
        return None

    N = len(active_sites)

    if N <= K:
        return [[a.name] for a in active_sites]


    if sim_matrix is None:
        sim_matrix = create_similarity_matrix(active_sites, "hydrophobicity")

    # initialize clustering
    k_ii = np.random.choice(N, size=K, replace=False)

    prev_assignment = None
    _iter = 0
    while _iter < MAX_ITER:
        curr_assignment = [[] for i in range(K)]

        for a in range(N):
            clust = np.argmin(sim_matrix[a, k_ii])
            curr_assignment[clust].append(a)

        # If the assignment is stable, return the current assignment
        if check_prev(prev_assignment, curr_assignment):
            return curr_assignment

        k_ii = find_centroids(curr_assignment, sim_matrix, K)
        prev_assignment = curr_assignment

        _iter += 1

    return curr_assignment


def cluster_hierarchically(active_sites, sim_matrix=None, linkage="single", K=10):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    def find_clusters_to_merge(clusts, sim_matrix):

        C = len(clusts)
        linkage_mat = np.zeros((C,C))
        if linkage == "single":

            # Find single linkage between all clusters
            for i in range(C):
                for j in range(C):
                    if i == j:
                        val = float("inf")
                    else:
                        val = np.min(sim_matrix[np.ix_(clusts[i], clusts[j])])
                    linkage_mat[i, j] = val

        to_merge = np.unravel_index(np.argmin(linkage_mat), (C, C))
        return to_merge

    if active_sites is None:
        return None

    # Fill in your code here!
    N = len(active_sites)
    K = min(K, N)

    if sim_matrix is None:
        sim_matrix = create_similarity_matrix(active_sites, "hydrophobicity")

    clusters = []
    for a in range(N):
        clusters.append([a])

    C = len(clusters)

    while C > K:

        merge_ii = find_clusters_to_merge(clusters, sim_matrix)
        clusters[merge_ii[0]] = list(clusters[merge_ii[0]] + clusters[merge_ii[1]])
        clusters = [clusters[i] for i in range(C) if i != merge_ii[1]]
        C = len(clusters)

    return clusters


def create_similarity_matrix(active_sites, comparator="RMSD"):
    """

    :param active_sites: List of ActiveSite objects of length N
    :param comparator: str, method to use when comparing active sites
    :return: np.array, NxN similarity matrix
    """
    A = len(active_sites)

    sim_mat = np.zeros((A, A))

    print("Creating Similarity Matrix")
    for i in tqdm(range(A)):
        for j in range(A):
            if i == j:
                sim_mat[i, j] = 0.0

            else:
                sim_mat[i, j] = compute_similarity(active_sites[i], active_sites[j], comparator)

    return sim_mat

