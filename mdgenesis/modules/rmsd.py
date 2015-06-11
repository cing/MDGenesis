# RMSD analysis class
import numpy as np
import pandas as pd
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

    def _loadcheckpoint(self, framedata, intdata):
        self.intdata = intdata
        if framedata.empty:
            self.framedata = pd.DataFrame(columns=["rmsd"])
        else:
            self.framedata = framedata

    def process(self, frame, frameid):
        a = self._selection.atoms.coordinates()
        b = self._refselection.atoms.coordinates()
        rmsd = np.sqrt(np.sum(np.power(a-b,2))/a.shape[0])
        rmsd_df = pd.DataFrame(rmsd, columns=["rmsd"], index=[frameid])
        self.framedata = self.framedata.append(rmsd_df)
        return True

    def _update_selections(self):
        self._selection = self.u.selectAtoms(self._selection_str)
        self._refselection = self.ref.selectAtoms(self._selection_str)

class RMSF(PerFrameAnalysis):

    # Fancy indexing in numpy lets you specify a list of whatever
    # columns you want. Just pass [0,1] to get X,Y and pass [2] to get Z.
    def __init__(self, selection_str, coordinate_indices=range(3)):
        """Selection string describes the atoms to be selected for RMSF analysis.
        coordinate range is a list of sequential indices.
        """
        self._selection_str = selection_str

        self._ci = coordinate_indices
        self._clabels = [["x","y","z"][ind] for ind in coordinate_indices]
        self._rmsf_sum_labels = ["rmsf_sum_"+lbl for lbl in self._clabels]
        self._rmsf_sqsum_labels = ["rmsf_sqsum_"+lbl for lbl in self._clabels]

    def _loadcheckpoint(self, framedata, intdata):

        self.framedata = framedata
        if intdata.empty:
            # We add one to store a resid index labelled zero for total_frames storage
            self.intdata = pd.DataFrame(np.zeros([len(self._selection)+1,
                                                  2*len(self._ci)]),
                                        columns=self._rmsf_sum_labels+self._rmsf_sqsum_labels,
                                        index=np.hstack([np.array([0]),self._sids]))

            self.intdata['total_frames'] = pd.Series(np.zeros([len(self._selection)+1]),
                                                     index=np.hstack([np.array([0]),self._sids]))

        else:
            self.intdata = intdata

    def results(self):
        """ Post-processing is needed to return the results since they are not
            appended. """

        f = float(self.intdata["total_frames"][0])

        # dropna is needed to drop the damn frame count row with resid 0
        rmsf_sqsum=self.intdata[self._rmsf_sqsum_labels].dropna()/f
        rmsf_sum=(self.intdata[self._rmsf_sum_labels].dropna()/f)**2
        rmsf_sqsum.columns=self._clabels
        rmsf_sum.columns=self._clabels

        return np.sqrt((rmsf_sqsum - rmsf_sum).sum(1))

    def process(self, frame, frameid):
        self.intdata.ix[0, "total_frames"] += 1
        self.intdata[self._rmsf_sum_labels] += pd.DataFrame(self._selection.coordinates()[:,self._ci],
                                                 index=self._sids,
                                                 columns=self._rmsf_sum_labels)
        self.intdata[self._rmsf_sqsum_labels] += pd.DataFrame(self._selection.coordinates()[:,self._ci]**2,
                                                   index=self._sids,
                                                   columns=self._rmsf_sqsum_labels)

        return True

    def _update_selections(self):
        self._selection = self.u.selectAtoms(self._selection_str)
        # We don't use resids because our protein has multiple subunits, it would
        # make a lot more sense to do that though...
        self._sids = np.array([sel.number for sel in self._selection])

class RadiusOfGyration(PerFrameAnalysis):

    def __init__(self, selection_str):
        """Keep in mind that this does not do alignment and the selections
            are fixed (as they should be), it's the laziest RGYR
            script of all time!
        """
        self._selection_str = selection_str

    def _loadcheckpoint(self, framedata, intdata):
        self.intdata = intdata
        if framedata.empty:
            self.framedata = pd.DataFrame(columns=["rgyr"])
        else:
            self.framedata = framedata

    def process(self, frame, frameid):
        rgyr_df = pd.DataFrame(self._selection.radiusOfGyration(),
                               columns=["rgyr"], index=[frameid])
        self.framedata = self.framedata.append(rgyr_df)
        return True

    def _update_selections(self):
        self._selection = self.u.selectAtoms(self._selection_str)
