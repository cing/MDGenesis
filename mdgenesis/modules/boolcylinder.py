# Cylinder analysis class
import numpy as np
import pandas as pd
from numpy.linalg import norm
from collections import defaultdict

from analysis import PerFrameAnalysis

class MultiCylinderBoolHistogram(PerFrameAnalysis):

    def __init__(self, solutesel, refsel="protein",
                 radius=5, extension=0, soluteaxis=2,
                 histmin=-50, histmax=50, histbins=200, jointz=True):
        """Return frames where histogram bin occupancy occured in
           all given refsels.

        :Arguments:
          *solutesel*
            Solute selection string
          *refsel*
            list of cylinder C.O.M. strings
          *radius*
            radius of cylinder
          *extension*
            extension of cylinder centered at the C.O.M. value
          *soluteaxis*
            coordination of solute molecules to histogram
          *histmin*
            minimum value of your box, relative to refsel
          *histmax*
            maximum value of your box, relative to refsel
          *histbins*
            number of histogram bins along the solute axis of the box
          *jointz*
            use the z-value of refsel 1 as the C.O.M. for all refsels

        """

        self.solutesel = solutesel
        self.refsel = refsel
        self.refnames = ["cyl"+str(r) for r in range(len(refsel))]
        self.r2 = radius**2
        self.extension = extension

        self.saxis = soluteaxis
        self.osaxis = list(set(range(3))-set([self.saxis]))

        self.histbins = histbins
        self.histmin = histmin
        self.histmax = histmax
        self.jointz = jointz

        self._edges = np.linspace(self.histmin, self.histmax,
                                  num=self.histbins+1)[:self.histbins]

    # We don't need to initialize framedata because it's only made
    # when results() is called.
    def _loadcheckpoint(self, framedata, intdata):
        self.framedata = framedata
        self.intdata = intdata

    def process(self, frame, frameid):
        """ Process a single trajectory frame """

        self._update_selections()
        all_resids_per_bin = []

        for rc in self._ref_coms:
            resids_per_bin = defaultdict(list)

            # Find out what is inside a cylinder
            radius_bool = np.sum(((self._scoord - rc)**2)[:,self.osaxis],axis=1) <= self.r2
            height_bool1 = (self._scoord - rc)[:,self.saxis] < self._height*0.5
            height_bool2 = (self._scoord - rc)[:,self.saxis] > -self._height*0.5
            cylinder_bool = radius_bool & height_bool1 & height_bool2

            # Determine what histogram bin it should goto
            bin_digit = np.digitize(self._scoord[cylinder_bool][:,self.saxis] -
                                    rc[self.saxis], self._edges)

            # Track only 1 residue per bin
            for resid, bin in zip(self._sids[cylinder_bool], bin_digit):
                resids_per_bin[bin] = resid
                #resids_per_bin[bin].append(resid)

            all_resids_per_bin.append(resids_per_bin)

        # Construct a row if resids exist in the same bin of all selections.
        unique_bins = set([res for ar in all_resids_per_bin for res in ar.keys()])
        temp_intdata = []
        for bin in sorted(unique_bins):
            if all([bin in resids for resids in all_resids_per_bin]):
                resids = [resids[bin] for resids in all_resids_per_bin]
                temp_intdata.append([frameid, self._edges[bin]] + resids)

        self.framedata = self.framedata.append(pd.DataFrame(temp_intdata,
                                                            columns=["frameid","bin"]+self.refnames))

        return True

    def _update_selections(self):
        self._ref_coms = np.array([self.u.selectAtoms(rs).centerOfMass() for rs in self.refsel])
        if self.jointz:
            self._ref_coms[:, self.saxis]=self._ref_coms[0, self.saxis]
        self._height = self.extension

        self._sids = np.array([sol.resid for sol in self.u.selectAtoms(self.solutesel)])
        self._scoord = self.u.selectAtoms(self.solutesel).coordinates()