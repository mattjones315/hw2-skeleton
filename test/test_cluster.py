from hw2skeleton import cluster
from hw2skeleton import io
from hw2skeleton import compare_clusters
import os


def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")

    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    assert cluster.compute_similarity(activesite_a, activesite_b, comparator="hydrophobicity") == 11.1556


def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    assert cluster.cluster_by_partitioning(active_sites) == [["276"], ["4629"], ["10701"]]

    # check empty active sites doesn't crash
    assert cluster.cluster_by_partitioning(None) is None


def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))
    assignment = cluster.cluster_hierarchically(active_sites)

    assert cluster.convert_indices_to_active_sites(assignment, active_sites) == [["276"], ["4629"], ["10701"]]

    pdb_ids = [276, 4629]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb" % id)
        active_sites.append(io.read_active_site(filepath))

    # check empty active sites doesn't crash
    assert cluster.cluster_hierarchically(None, K=10) is None


def test_rand():

    pdb_ids = [276, 4629, 10701, 10814, 13052, 14181]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb" % id)
        active_sites.append(io.read_active_site(filepath))

    c1 = [["276", "4629"], ["10701", "10814"], ["13052", "14181"]]
    c2 = [["10701", "10814"], ["276", "4629"], ["13052", "14181"]]

    assert compare_clusters.rand_index(c1, c2, active_sites) == 1
