# Coordination of solute in cylinder analysis class
import numpy as np
import pandas as pd
from numpy.linalg import norm
import mdtraj as md

from analysis import AllAtOnceAnalysis

def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays. (Faster than itertools)

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

class FirstShellCoordinationCylinder(AllAtOnceAnalysis):

    def __init__(self, solutesel, refsel="protein",
                 ligand_sel=None, ligand_colnames=None,
                 rdf_cut=0.3,
                 radius=1.0,):
        """Compute histogram of solute molecules coordinates (X, Y, or Z) within
        a cylinder along a given axis.

        :Arguments:
          *solutesel*
            Solute selection string
          *refsel*
            cylinder C.O.M. string
          *ligand_sel*
            list of selections that represent coordination ligands to solute
          *rdf_cut*
            distance between ligands to be considered hydrogen-bonded
          *radius*
            radius of cylinder

        """

        self.solutesel = solutesel
        self.refsel = refsel
        self.ligand_sel = ligand_sel
        if ligand_colnames != None:
            self.ligand_colnames = ligand_colnames
        else:
            self.ligand_colnames = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:len(ligand_sel)]
        self.rdf_cut = rdf_cut
        self.r2 = radius**2

    def results(self):

        # This is a helper method that just does an array lookup
        def _fill_coord(row):
            return ionprot_dists_sum[int(row["Frame"]-1), int(row["IonIDX"])]

        # Extract the positions of all indices and take the mean
        ref_xyz = self.trj.xyz[:,self.ref_idx,:].mean(axis=1)

        # Subtract the reference point (swapaxes is needed for broadcast)
        ion_xyz = self.trj.xyz[:,self.solute_idx,:].swapaxes(0,1) - ref_xyz

        # Do a boolean check to determine which ions are within the cylinder of radius 1.0
        ion_xyz2 = ion_xyz**2
        inside_cylinder = (ion_xyz2[:,:,0]+ion_xyz2[:,:,1]) < self.r2

        # Extract the unique ID for each ion
        ion_index, time_index = np.where(inside_cylinder)

        ion_index_in_whole_system = self.solute_idx[ion_index]
        resids = [self.trj.topology.atom(i).residue.resSeq for i in ion_index_in_whole_system]

        # Construct a dataframe with this positional data
        ion_positions_df = pd.DataFrame(ion_xyz[inside_cylinder]*10, columns=["X","Y","Z"])
        ion_positions_df["ResidID"] = resids
        ion_positions_df["IonIDX"] = ion_index
        ion_positions_df["Frame"] = time_index+1
        time_sorted_ions = ion_positions_df.sort_values("Frame") #.reset_index(drop=True)

        for lig_name, lig_subset in zip(self.ligand_colnames, self.ligand_idx):
            ionprot_pairs=cartesian((self.solute_idx, lig_subset))
            ionprot_dists = md.compute_distances(self.trj, atom_pairs=ionprot_pairs)
            ionprot_dists_sum = (ionprot_dists.reshape(ionprot_dists.shape[0],
                                                       len(self.solute_idx), -1) < self.rdf_cut).sum(axis=2)
            time_sorted_ions[lig_name]=time_sorted_ions.apply(_fill_coord, axis=1)

        time_sorted_ions.drop('IonIDX', axis=1, inplace=True)

        return time_sorted_ions

    def _update_selections(self):
        self.solute_idx=self.trj.topology.select(self.solutesel)
        self.ref_idx=self.trj.topology.select(self.refsel)
        self.ligand_idx=[self.trj.topology.select(sel) for sel in self.ligand_sel]

