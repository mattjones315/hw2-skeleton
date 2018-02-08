# Some utility classes to represent a PDB structure
import numpy as np

BACKBONE = [0,1,2,]

class Atom:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type):
        self.type = type
        self.coords = (0.0, 0.0, 0.0)

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return "{0}, {1}".format(self.type, self.coords)
        #return self.type

class Residue:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type, number):
        self.type = type
        self.number = number
        self.atoms = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return "{0}\t{1}\n".format(self.type, str(self.atoms))
        #return "{0} {1}".format(self.type, self.number)

    def getCoords(self):
        return [a.coords for a in self.atoms]

    def getBackboneCoords(self):

        return [self.atoms[i].coords for i in BACKBONE]

    def rmsd(self, r2):
        c1 = self.getBackboneCoords()
        c2 = r2.getBackboneCoords()
        dist_vec = np.array([euclideanDist(np.array(c1[i]), np.array(c2[i])) for i in range(len(c1))])
        rmsd = np.sum(dist_vec) / len(c2)
        return np.sqrt(rmsd)


class ActiveSite:
    """
    A simple class for an active site
    """

    def __init__(self, name):
        self.name = name
        self.residues = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.name

    def getResidues(self):
        return self.residues

def euclideanDist(c1, c2):
    """

    :param c1: list, set of atomic coordinates
    :param c2: list, second set of coordinates to compare to
    :return: float, euclidean distance between c1 and c2
    """

    return np.sqrt(np.sum( (c1 - c2)**2 ))


