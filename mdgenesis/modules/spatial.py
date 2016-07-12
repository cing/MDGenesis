from analysis import PerFrameAnalysis
from itertools import product
import numpy as np
import pandas as pd
import logging
logger = logging.getLogger('distance_between_points')

class DistanceBetweenPoints(PerFrameAnalysis):
    """ Calculate distance between the C.O.M. of two selections
    """

    def __init__(self, selection1, selection2):
        """Calculate distance between two selections.

        :Arguments:
          *selection1*
            Selection string for first selection
          *selection2*
            Selection string for second selection

        """

        self.selection1 = selection1
        self.selection2 = selection2

        if not (self.selection1 and self.selection2):
            raise Exception('DistanceBetweenPoints: invalid selections')

    def _loadcheckpoint(self, framedata, intdata):
        self.intdata = intdata
        if framedata.empty:
            self.framedata = pd.DataFrame(columns=["dist"])
        else:
            self.framedata = framedata

    def process(self, frame, frameid):
        r = self._s1.center_of_mass() - self._s2.center_of_mass()
        d = np.sqrt(np.sum(r*r))
        dist_df = pd.DataFrame(d, columns=["dist"], index=[frameid])
        self.framedata = self.framedata.append(dist_df)
        return True

    def _update_selections(self):
        self._s1 = self.u.selectAtoms(self.selection1)
        self._s2 = self.u.selectAtoms(self.selection2)

class CenterOfMassPosition(PerFrameAnalysis):
    """ Calculates the coordinates of an arbitrary number of selections
        relative to some selection.
    """

    def __init__(self, selections, refsel=None, selaxis=range(3)):
        """Calculate distance between two selections.

        :Arguments:
          *selections*
            A list of atom selection strings
          *refsel*
            A selection for which the C.O.M. will provide or a list of selections
          *selaxis*
            Returns a specific axis (0,1,2) for (X,Y,Z) rather than all coordinates.

        """

        if len(selections) == 0:
            raise Exception('Position: no selection strings provided')

        self.selections = selections
        self.selnames = ["pos"+str(r) for r in range(len(selections))]
        self.refsel = refsel
        self.saxis = selaxis
        self._clabels = [["_x","_y","_z"][ind] for ind in selaxis]
        self._poslabels = ["".join(p) for p in list(product(self.selnames,
                                                            self._clabels))]

    def _loadcheckpoint(self, framedata, intdata):
        self.intdata = intdata
        if framedata.empty:
            self.framedata = pd.DataFrame(columns=self._poslabels)
        else:
            self.framedata = framedata

    def process(self, frame, frameid):
        """ Process a single trajectory frame """
        if self.refsel == None:
            rel_pos = [a.center_of_mass() for a in self._selection_atoms]
        elif len(self.refsel) == len(self._selection_atoms):
            rel_pos = [a.center_of_mass()-ref.center_of_mass()
                       for a,ref in zip(self._selection_atoms, self._refsel_atoms)]
        else:
            ref_com = self._refsel_atoms.center_of_mass()
            rel_pos = [a.center_of_mass()-ref_com for a in self._selection_atoms]

        p = np.hstack([p[self.saxis] for p in rel_pos]).reshape(1,len(self._poslabels))
        pos_df = pd.DataFrame(p, columns=self._poslabels, index=[frameid])

        self.framedata = self.framedata.append(pos_df)
        return True

    def _update_selections(self):
        self._selection_atoms = [self.u.selectAtoms(sel) for sel in self.selections]

        if len(self.refsel) == len(self.selections):
            self._refsel_atoms = [self.u.selectAtoms(ref) for ref in self.refsel]
        elif self.refsel != None:
            self._refsel_atoms = self.u.selectAtoms(self.refsel)
