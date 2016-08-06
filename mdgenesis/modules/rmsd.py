# RMSD analysis class
import numpy as np
import pandas as pd
#import numpy.linalg

#from MDAnalysis import *
#from MDAnalysis.core.AtomGroup import Residue, AtomGroup
import MDAnalysis.analysis.rms

from analysis import PerFrameAnalysis
from analysis import AllAtOnceAnalysis

class MultiRMSD(AllAtOnceAnalysis):

    def __init__(self, ref, select="protein and name CA",
                 group_selections=None, group_names=None):
        """ RMSD is computed after an alignment on "select" to the reference
            and then RMSD is computed for all other selections.
        """
        self.select = select
        self.group_selections = group_selections
        self.ref = ref

        self.colnames=["frame","time","align"]
        if group_selections != None:
            if group_names != None:
                if len(group_names) == len(group_selections):
                    self.colnames.extend(group_names)
                else:
                    ngroups = len(group_selections)
                    lcols = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:ngroups]
                    self.colnames.extend(list("ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:ngroups]))
            else:
                ngroups = len(group_selections)
                lcols = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:ngroups]
                self.colnames.extend(list("ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:ngroups]))

    def results(self):
        R = MDAnalysis.analysis.rms.RMSD(self.u, self.ref,
                                         select=self.select,
                                         groupselections=self.group_selections)
        R.run()
        return pd.DataFrame(R.rmsd, columns=self.colnames)

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
