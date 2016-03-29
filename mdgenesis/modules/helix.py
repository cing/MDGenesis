import mdtraj as md
import pandas as pd
import numpy as np
from analysis import AllAtOnceAnalysis

def angle_between_vectors(A, B):
    """
    Returns the angle theta between two vectors A and B
    A . B = |a| * |b| * cos(theta)
    """
    A_dot_B = np.einsum('ij,ij->i', A, B)
    a = np.linalg.norm(A, axis=1)
    b = np.linalg.norm(B, axis=1)
    theta = np.degrees(np.arccos(A_dot_B/a/b))
    #theta[theta>90] *= -1
    #theta[theta>90] += 180
    #return np.abs(theta)
    return theta

def vector_projection_onto_plane(v, n):
    """
    Returns the vector v_prime obtained by orthogonally projecting a
    vector v onto a plane defined by vectors A and B
    n = A x B
    a = n x (v X n)
    v_prime = [(a . v) / (a . a)] * a
    """
    # Vector orthogonal to plane
    #n = np.cross(A, B)
    # Vector in the direction of the projection
    a = np.cross(n, np.cross(v, n))
    # Projected vector
    temp_dot = np.einsum('ij,ij->i', a, v)
    temp_norm = np.einsum('ij,ij->i', a, a)
    v_prime = (temp_dot/temp_norm)[:, np.newaxis] * a

    return v_prime

# Take a sliced trajectory that only contains the ca atoms of a helix
# and returns a principal axis vector
def principal_axis(traj):
    mat = md.compute_inertia_tensor(traj)
    evals,evecs = np.linalg.eig(mat) #mat[0] is the 0th frame, must be modified to rapidly compute for whole traj

    # Need to make everything N*3 shape
    l = evals.shape[0]

    # Need the index of the largest eigenvalue only, in N*3 shape
    idx = np.argsort(evals)[:,0][:, np.newaxis]
    idx2 = np.arange(l)[:, np.newaxis]
    idx3 = (np.array([[0, 1, 2]]*l))

    return -evecs[idx2, idx3, idx]

class HelixTilt(AllAtOnceAnalysis):

    def __init__(self, helix_selection, principal_axis_selection=None):
        """ Helix tilt is computed for all helical selections in the
            atom selection list.
        """

        self._paxis_selection = principal_axis_selection
        self._helix_selection = helix_selection

    def results(self):
        # Determine the principal axis of the channel
        if self._paxis_selection != None:
            p_slice = self.trj.atom_slice(self._principal_axis_atoms)
            paxis = principal_axis(p_slice)
        else:
            paxis = np.array([[0,0,1]]*len(self.trj))

        temp_results = []
        for helix_atoms in self._helix_atoms:
            a_slice = self.trj.atom_slice(helix_atoms)
            haxis = principal_axis(a_slice)
            temp_results.append(angle_between_vectors(paxis, haxis))

        return pd.DataFrame(temp_results)

    def _update_selections(self):
        if self._paxis_selection != None:
            self._principal_axis_atoms = self.trj.topology.select(self._paxis_selection)
        self._helix_atoms = [self.trj.topology.select(a) for a in self._helix_selection]

class HelixRotation(AllAtOnceAnalysis):

    def __init__(self, helix_selection, reference_selection, principal_axis_selection=None):
        """ Helix rotation is computed for a helical selection with respect
            to the center of mass of the principal axis and a reference atom.
            If helix selection includes multiple atoms, a center of mass is calculated
        """

        self._reference_selection = reference_selection
        self._paxis_selection = principal_axis_selection
        self._helix_selection = helix_selection

    def results(self):
        p_slice = self.trj.atom_slice(self._principal_axis_atoms)
        paxis = principal_axis(p_slice)
        paxis_com = md.compute_center_of_mass(p_slice)

        temp_results = []
        for ref_atom, helix_atoms in zip(self._reference_atoms, self._helix_atoms):
            if len(ref_atom) > 1:
                ref_pos = self.trj.xyz[:,ref_atom,:].mean(axis=1)
            else:
                ref_pos = self.trj.xyz[:,ref_atom[0],:]

            helix_com = self.trj.xyz[:, helix_atoms, :].mean(axis=1)

            proj_of_com = vector_projection_onto_plane(paxis_com - helix_com, paxis)
            proj_of_ref = vector_projection_onto_plane(ref_pos - helix_com, paxis)

            a = angle_between_vectors(proj_of_com, proj_of_ref)
            temp_results.append(a)

        return pd.DataFrame(temp_results)

    def _update_selections(self):
        if self._paxis_selection != None:
            self._principal_axis_atoms = self.trj.topology.select(self._paxis_selection)
        else:
            self._principal_axis_atoms = self.trj.topology.select("protein and name CA")

        self._helix_atoms = [self.trj.topology.select(a) for a in self._helix_selection]
        self._reference_atoms = [self.trj.topology.select(a) for a in self._reference_selection]
