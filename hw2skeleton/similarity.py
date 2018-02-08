from .utils import Atom, Residue, ActiveSite
import numpy as np

HYDROPHOBICITY = {"ILE": .73, "PHE": .61, "VAL": .54, "LEU": .53,
                  "TRP": .37, "MET": .26, "ALA": .25, "GLY": .16,
                  "CYS": .04, "TYR": .02, "PRO": -.07, "THR": -0.18,
                  "SER": -.26, "HIS": -.40, "GLU": -0.62, "ASN": -.64,
                  "GLN": -0.69, "ASP": -0.72, "LYS": -1.10, "ARG": -1.80}


class SimilarityComparator:
    """
    Creates a similarity comparator, defined by METHOD. The object will be able to compare two active sites
    and return a float representing the degree of similarity between the two active site sequences.

    """

    def __init__(self, method = "RMSD"):

        self.method = method

    def compare(self, a1, a2):
        """
        Compares two active sites according to SELF.METHOD.

        :param a1: Active site 1
        :param a2: Active site 2
        :return: float, degree of similarity between a1 and a2
        """

        if self.method == "RMSD":

            return self.RMSD(a1, a2)

        elif self.method == "hydrophobicity":
            return self.score_hydrophobicity(a1, a2)

        raise ValueError("Similarity method not recognized.")

    # MUST FIX
    def TM(self, a1, a2):
        """
        Performs an iterative structural mapping algorithm to find the best scoring alignment, which will be used
        to calculate the similarity score between active sites. Adapted from MaxCluster protein structure comparison
        tool (http://www.sbg.bio.ic.ac.uk/maxcluster/index.html)

        :param a1: ActiveSite object 1
        :param a2: ActiveSite object 2
        :return: returns a optimized sequence independent alignment score
        """
        r1 = a1.getResidues()
        r2 = a2.getResidues()

        N = min(len(r1), len(r2))

        distmat = np.zeros((len(r1), len(r2)))

        for i in range(len(r1)):
            for j in range(len(r2)):

                dij = r1[i].rmsd(r2[j])
                distmat[i, j]= dij

        print(distmat.shape)

    def RMSD(self, a1, a2, _type="min"):
        """
        Performs a naive RMSD search, and returns a single RMSD value (min, max, average, sum, or median). Will search over
        all possible alignments, calculate total RMSD, and then return the _TYPE of RMSD specified.

        :param a1: ActiveSite object 1
        :param a2: ActiveSite object 2
        :param _type: str, type of RMSD to return. Supported operations: "min", "max", "average", "median"
        :return: returns _type(RMSD) between a1 and a2
        """
        r1 = a1.getResidues()
        r2 = a2.getResidues()

        if len(r1)  < len(r2):
            _short = r1
            _long = r2
        else:
            _short = r2
            _long = r1

        N = len(_short)
        M = len(_long)
        all_rmsds = []

        # We have M-N possible alignments
        for i in range(M - N + 1):
            curr_rmsd = 0.0

            # Calculate RMSD for the current alignment
            for j in range(N):
                curr_rmsd += _short[j].rmsd(_long[i+j])

            # Add current RMSD to set of alignment RMSDs
            all_rmsds.append(curr_rmsd)

        if _type == "min":
            return min(all_rmsds)

        elif _type == "max":
            return max(all_rmsds)

        elif _type == "average":
            return np.mean(all_rmsds)

        elif _type == "sum":
            return np.sum(all_rmsds)

        elif _type == "median":
            return np.median(all_rmsds)

        else:
            raise ValueError("Type of RMSD not recognized")

    def score_hydrophobicity(self, a1, a2):
        """
        Comparison based on charge and hydrophobicity of residues of each active site.


        :param a1: ActiveSite object 1
        :param a2: ActiveSite object 2
        :return: float, similarity between a1 and a2
        """

        # Convert each active site sequence to a hydrophbocity score

        h1 = 0.0
        h2 = 0.0

        r1 = a1.getResidues()
        r2 = a2.getResidues()

        for r in r1:
            h1 += HYDROPHOBICITY[r.type]

        for r in r2:
            h2 += HYDROPHOBICITY[r.type]

        return (h1 - h2)**2


