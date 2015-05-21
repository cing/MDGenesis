from analysis import PerFrameAnalysis
import numpy
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

    def process(self, frame):
        """ Process a single trajectory frame """
        r = self._s1.centerOfMass() - self._s2.centerOfMass()
        d = numpy.sqrt(numpy.sum(r*r))
        self.framedata.append(d)

    def _update_selections(self):
        self._s1 = self.u.selectAtoms(self.selection1)
        self._s2 = self.u.selectAtoms(self.selection2)

class CenterOfMassPosition(PerFrameAnalysis):
    """ Calculates the coordinates of an arbitrary number of selections
        relative to some selection.
    """

    def __init__(self, selections, refsel=None, selaxis=None):
        """Calculate distance between two selections.

        :Arguments:
          *selections*
            A list of atom selection strings
          *refsel*
            A selection for which the C.O.M. will provide
          *selaxis*
            Returns a specific axis (0,1,2) for (X,Y,Z) rather than all coordinates.

        """

        if len(selections) == 0:
            raise Exception('Position: no selection strings provided')

        self.selections = selections
        self.refsel = refsel
        self.selaxis = selaxis

    def process(self, frame):
        """ Process a single trajectory frame """
        if self.refsel != None:
            ref_com = self._refsel_atoms.centerOfMass()
            rel_pos = [a.centerOfMass()-ref_com for a in self._selection_atoms]
        else:
            rel_pos = [a.centerOfMass() for a in self._selection_atoms]

        if self.selaxis != None:
            self.framedata.append([pos[self.selaxis] for pos in rel_pos])
        else:
            self.framedata.append([coord for pos in rel_pos for coord in pos])

    def _update_selections(self):
        self._selection_atoms = [self.u.selectAtoms(sel) for sel in self.selections]
        if self.refsel != None:
            self._refsel_atoms = self.u.selectAtoms(self.refsel)

