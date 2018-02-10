import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically, create_similarity_matrix, \
        convert_indices_to_active_sites
from .compare_clusters import sum_of_distances, rand_index, benchmark_clusters, benchmark_rand, rand_versus
import pickle as pic

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

active_sites = read_active_sites(sys.argv[2])

sim_matrix = create_similarity_matrix(active_sites, "hydrophobicity")
pic.dump(sim_matrix, open("sim_matrix.pkl", "wb"))

#sim_matrix = pic.load(open("sim_matrix.pkl", "rb"))

# Choose clustering algorithm
if sys.argv[1][0:3] == "-PC":
    print("Benchmarking Partitioning Algorithm")
    rmsd_ssds = benchmark_clusters(active_sites, alg="P", metric="RMSD", N=1000)
    hyd_ssds = benchmark_clusters(active_sites, alg="P", metric="hydrophobicity", N=1000)

    pic.dump(rmsd_ssds, open("kmeans_rmsd_ssds.pkl", "wb"))
    pic.dump(hyd_ssds, open("kmeans_hyd_ssds.pkl", "wb"))

if sys.argv[1][0:3] == "-HC":
    print("Benchmarking Hierarchical Algorithm")
    rmsd_ssds = benchmark_clusters(active_sites, alg="H", metric="RMSD", N=50)
    hyd_ssds = benchmark_clusters(active_sites, alg="H", metric="hydrophobicity", N=50)

    pic.dump(rmsd_ssds, open("kmeans_rmsd_ssds.pkl", "wb"))
    pic.dump(hyd_ssds, open("kmeans_hyd_ssds.pkl", "wb"))

if sys.argv[1][0:3] == "-PR":
    print("Computing Average Rand Index, K-means")
    rmsd_rand = benchmark_rand(active_sites, alg="P", metric="RMSD", N=30, step=10)
    hyd_rand = benchmark_rand(active_sites, alg="P", metric="hydrophobicity", N=30, step=10)

    pic.dump(rmsd_rand, open("kmeans_rmsd_rand.pkl", "wb"))
    pic.dump(hyd_rand, open("kmeans_hyd_rand.pkl", "wb"))

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
    hyd_rand = rand_versus(active_sites, N=20, metric="hydrophobicity")
    rmsd_rand = rand_versus(active_sites, N=20, metric="RMSD")


    pic.dump(rmsd_rand, open("compare_rmsd_rand.pkl", "wb"))
    pic.dump(hyd_rand, open("compare_hyd_rand.pkl", "wb"))


