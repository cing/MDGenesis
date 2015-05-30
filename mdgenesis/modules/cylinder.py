# Cylinder analysis class
import numpy as np
import pandas as pd
from numpy.linalg import norm

from analysis import PerFrameAnalysis

class CylinderHistogram(PerFrameAnalysis):

    def __init__(self, topsel, bottomsel, solutesel, refsel="protein",
                 radius=5, extension=0, soluteaxis=2,
                 histmin=-50, histmax=50, histbins=200):
        """Compute histogram of solute molecules coordinates (X, Y, or Z) within
        a cylinder defined by a top and bottom selection and a cylinder radius.

        :Arguments:
          *topsel*
            Top of cylinder selection string
          *bottomsel*
            Bottom of cylinder selection string
          *solutesel*
            Solute selection string
          *comsel*
            cylinder C.O.M. string
          *radius*
            radius of cylinder
          *extension*
            extension of cylinder
          *soluteaxis*
            coordination of solute molecules to histogram
          *histmin*
            minimum value of your box, relative to refsel
          *histmax*
            maximum value of your box, relative to refsel
          *histbins*
            number of histogram bins along the solute axis of the box

        """

        self.topsel = topsel
        self.bottomsel = bottomsel
        self.solutesel = solutesel
        self.refsel = refsel
        self.r = radius
        self.extension = extension
        self.saxis = soluteaxis

        self.histbins = histbins
        self.histmin = histmin
        self.histmax = histmax

        #self._all_edges = np.linspace(histmin, histmax, num=histbins)

    # We don't need to initialize framedata because it's only made
    # when results() is called.
    def _loadcheckpoint(self, frames_processed, framedata, intdata):
        if len(intdata) == 0:
            self.frames_processed = 0
            self.intdata = pd.DataFrame(np.zeros(self.histbins, dtype=np.uint64),
                                        columns=["bincount"])
        else:
            self.frames_processed = frames_processed
            self.framedata = framedata
            self.intdata = intdata

    def process(self, frame):
        """ Process a single trajectory frame """
        self._update_selections()
        self.frames_processed += 1

        # distance to the helix axis defined by top_com and bottom_com
        d = norm(np.cross(self._scoord - self._top_com,
                          self._scoord - self._bottom_com),
                 axis=1)/self._height

        # boolean mask for distances below the radius, hist'd on solute-axis
        h = np.histogram(self._scoord[d <= self.r][:,self.saxis] - self._ref_com[self.saxis],
                         range=[self.histmin, self.histmax],
                         bins=self.histbins, normed=False)[0]

        # optionally output the resids of solute in cylinder (useful for VMD)
        #print self._sids[dist <= self.radius]

        self.intdata += pd.DataFrame(h, columns=["bincount"])

    def results(self):
        if self.frames_processed > 0:
            return self.intdata/float(self.frames_processed)
        else:
            return self.intdata

    def _update_selections(self):
        self._top_com = self.u.selectAtoms(self.topsel).centerOfMass()
        self._bottom_com = self.u.selectAtoms(self.bottomsel).centerOfMass()
        self._ref_com = self.u.selectAtoms(self.refsel).centerOfMass()

        self._midpoint = (self._top_com+self._bottom_com)/2.0
        self._height = np.linalg.norm(self._top_com-self._bottom_com)
        self._sradius = (self._height/2) + self.extension

        #self._sids = np.array([sol.resid for sol in self.u.selectAtoms(self.solutesel)])
        self._scoord = self.u.selectAtoms(self.solutesel).coordinates()