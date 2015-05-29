# RMSD analysis class
import numpy as np
#import numpy.linalg

#from MDAnalysis import *
#from MDAnalysis.core.AtomGroup import Residue, AtomGroup
#import MDAnalysis.core.rms_fitting
from analysis import PerFrameAnalysis

'''
# RMSD with alignment, written by David Caplan
class FrameData(object):
    atoms = None
    coordinates = None
    masses = None
    com = None
    rmsd = None

    def __init__(self, atom_group, allocate_only=False):
        self.atoms = atom_group
        if allocate_only:
            self.coordinates = atom_group.coordinates().copy()
        else:
            self.masses = atom_group.masses()
            self.com = atom_group.centerOfMass().astype(numpy.float32)
            self.coordinates = atom_group.coordinates() - self.com

class RMSD(object):
    _selection = None
    _rmsds = []

    def _rmsd(self, a,b):
        """Returns RMSD between two coordinate sets a and b."""
        return numpy.sqrt(numpy.sum(numpy.power(a-b,2))/a.shape[0])

    def __init__(self, selection):
        self._selection = selection

    def prepare(self, ref, trj):
        ref_atoms = ref.selectAtoms(self._selection)
        trj_atoms = trj.selectAtoms(self._selection)
        self.fit_ref = FrameData(ref.selectAtoms('backbone'))
        self.fit_trj = FrameData(trj.selectAtoms('backbone'))
        self.rmsd_ref = FrameData(ref_atoms)
        self.rmsd_trj = FrameData(trj_atoms, allocate_only=True)
        # print "Done RMSD prepare."

    def process(self, ts):
        # print "RMSD Fitting Frame %5d" % (ts.frame)
        x_com = self.fit_trj.atoms.centerOfMass().astype(numpy.float32)
        self.fit_trj.coordinates[:] = self.fit_trj.atoms.coordinates() - x_com
        R = numpy.matrix(MDAnalysis.core.rms_fitting.rms_rotation_matrix(self.fit_trj.coordinates, self.fit_ref.coordinates, self.fit_ref.masses),dtype=numpy.float32)
        ts._pos   -= x_com
        ts._pos[:] = ts._pos * R
        ts._pos   += self.fit_ref.com
        self._rmsds.append(self._rmsd(self.rmsd_ref.atoms.coordinates(),self.rmsd_trj.atoms.coordinates()))
        # print self._rmsds[-1]

    def results(self):
        return self._rmsds

'''

class RMSD(PerFrameAnalysis):

    def __init__(self, selection_str):
        """Keep in mind that this does not do alignment and the selections
            are fixed (as they should be), it's the laziest RMSD
        script of all time!
        """
        self._selection_str = selection_str

    def process(self, frame):
        self.frames_processed += 1
        a = self._selection.atoms.coordinates()
        b = self._refselection.atoms.coordinates()
        self.framedata.append(np.sqrt(np.sum(np.power(a-b,2))/a.shape[0]))

    def _update_selections(self):
        self._selection = self.u.selectAtoms(self._selection_str)
        self._refselection = self.ref.selectAtoms(self._selection_str)

# TODO: Make this computation resumable by storing coord_sum and coord_sqsum
#       in some kind of MDGenesis metadata columns. FrameData is currently unused.
class RMSF(PerFrameAnalysis):

    def __init__(self, selection_str, coordinate_indices=None):
        """Selection string describes the atoms to be selected for RMSF analysis.
        coordinate range is a list of sequential indices.
        """
        self._selection_str = selection_str
        self._used_frames = 0

        # Fancy indexing in numpy lets you specify a list of whatever
        # columns you want. Just pass [0,1] to get X,Y and pass [2] to get Z.
        if coordinate_indices is None:
            self._coordinate_indices = range(3)
        else:
            self._coordinate_indices = coordinate_indices

    def results(self):
        """ Post-processing is needed to return the results since they are not
            appended. """
        f = self.frames_processed
        ci = self._coordinate_indices

        if len(ci) > 0:
            return np.sqrt((self._coord_sqsum[:,ci]/f -
                           (self._coord_sum[:,ci]/f)**2).sum(1))
        else:
            return np.sqrt((self._coord_sqsum/f -
                           (self._coord_sum/f)**2).sum(1))

    def process(self,ts):
        self.frames_processed += 1
        self._coord_sum += self._selection.coordinates()
        self._coord_sqsum += self._selection.coordinates()**2

    def _update_selections(self):
        self._selection = self.u.selectAtoms(self._selection_str)
        self._coord_sum = np.zeros([len(self._selection), 3])
        self._coord_sqsum = np.zeros([len(self._selection), 3])

class RadiusOfGyration(PerFrameAnalysis):

    def __init__(self, selection):
        self._selection_str =  selection

    def process(self, frame):
        self.framedata.append(self._selection.radiusOfGyration())

    def _update_selections(self):
        self._selection = self.u.selectAtoms(self._selection_str)
