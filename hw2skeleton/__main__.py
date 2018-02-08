import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically, create_similarity_matrix, \
        convert_indices_to_active_sites
from .compare_clusters import sum_of_distances, rand_index
import pickle as pic

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

active_sites = read_active_sites(sys.argv[2])

sim_matrix = create_similarity_matrix(active_sites, "hydrophobicity")
#pic.dump(sim_matrix, open("sim_matrix.pkl", "wb"))

#im_matrix = pic.load(open("sim_matrix.pkl", "rb"))

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    assignments = cluster_by_partitioning(active_sites, sim_matrix)
    clustering = convert_indices_to_active_sites(assignments, active_sites)
    write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    assignments = cluster_hierarchically(active_sites, sim_matrix)
    clusterings = convert_indices_to_active_sites(assignments, active_sites)
    write_clustering(sys.argv[3], clusterings)

if sys.argv[1][0:2] == "-C":
    print("Comparing Kmeans to Hierarchical")
    assignments_1 = cluster_by_partitioning(active_sites, sim_matrix)
    assignments_2 = cluster_by_partitioning(active_sites, sim_matrix)
    clustering_1 = convert_indices_to_active_sites(assignments_1, active_sites)
    clustering_2 = convert_indices_to_active_sites(assignments_2, active_sites)

    print(rand_index(clustering_1, clustering_2, active_sites))

