# Cylinder analysis class
import MDAnalysis.KDTree.NeighborSearch as ns
import numpy as np
#import numpy.linalg

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
        self.radius = radius
        self.extension = extension
        self.soluteaxis = soluteaxis

        self.histbins = histbins
        self.histmin = histmin
        self.histmax = histmax

        self._frames_processed = 0.0
        self._all_histograms = np.zeros(histbins, dtype=np.uint64)
        self._all_edges = np.linspace(histmin, histmax, num=histbins)

    def process(self, frame):
        """ Process a single trajectory frame """
        self._update_selections()
        self._frames_processed += 1.0

        solute_in_sphere = set(self._solutesearch.search(self._midpoint,
                                                         self._sradius,
                                                         level="R"))

        solute_axis_within_cylinder = []
        # This is a processing step where we embed a cylinder within the sphere
        # that connects the top and bottom points (no matter what angle it is).
        # Keep in mind that solute center of mass is used no matter what the
        # selection is (since it happens on the residue level)
        for solute in solute_in_sphere:
            s_com = solute.centerOfMass()
            s_com_wrt_ref = s_com - self._ref_com
            s_dist_to_cylinder_axis = np.linalg.norm(np.cross(s_com - self._top_com, s_com - self._bottom_com))
            if s_dist_to_cylinder_axis/self._height <= self.radius:
                solute_axis_within_cylinder.append(s_com_wrt_ref[self.soluteaxis])

        self._all_histograms += np.histogram(np.array(solute_axis_within_cylinder),
                                             range=np.array([self.histmin, self.histmax]),
                                             bins=np.array(self.histbins), normed=False)[0]

    def results(self):
        if self._frames_processed > 0:
            return np.array(self._all_histograms)/self._frames_processed
        else:
            return np.array(self._all_histograms)

    def _update_selections(self):
        self._top_com = self.u.selectAtoms(self.topsel).centerOfMass()
        self._bottom_com = self.u.selectAtoms(self.bottomsel).centerOfMass()
        self._ref_com = self.u.selectAtoms(self.refsel).centerOfMass()

        self._midpoint = (self._top_com+self._bottom_com)/2.0
        self._height = np.linalg.norm(self._top_com-self._bottom_com)
        self._sradius = (self._height/2) + self.extension

        self._solute = self.u.selectAtoms(self.solutesel)
        self._solutesearch = ns.AtomNeighborSearch(self._solute)