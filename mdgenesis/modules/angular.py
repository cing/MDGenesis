from analysis import PerFrameAnalysis
import numpy as np
import pandas as pd
import logging
logger = logging.getLogger('dihedral_angle')

class DihedralAngle(PerFrameAnalysis):

    def __init__(self, selection1, selection2, selection3, selection4):
        """Calculate dihedral angle between four atoms in an atom selection

        :Arguments:
          *selection1*
            Selection string for first atom
          *selection2*
            Selection string for second atom
          *selection3*
            Selection string for third atom
          *selection4*
            Selection string for fourth atom

        """

        self.selection1 = selection1
        self.selection2 = selection2
        self.selection3 = selection3
        self.selection4 = selection4

        if not (self.selection1 and self.selection2 and self.selection3 and self.selection4):
            raise Exception('DihedralAngle: invalid selections')

    def process(self, frame, frameid):
        """ Process a single trajectory frame """
        dihe_df = pd.DataFrame(self._sel.dihedral.value(), columns=["dihe"], index=[frameid])
        self.framedata = self.framedata.append(dihe_df)
        return True

    def _update_selections(self):
        self._sel = self.u.selectAtoms(self.selection1)
        self._sel += self.u.selectAtoms(self.selection2)
        self._sel += self.u.selectAtoms(self.selection3)
        self._sel += self.u.selectAtoms(self.selection4)

class PhiPsiAngle(PerFrameAnalysis):

    def __init__(self, resnums, atoms=[["C","N","CA","C"],["N","CA","C","N"]]):
        """Calculates the phi and psi backbone angles based on an
           input selection of a residue.

        :Arguments:
          *resnums*
            A list of residue numbers within the protein

        """

        self.resnums = resnums
        self.atoms = atoms

    def process(self, frame, frameid):
        """ Process a single trajectory frame """
        angs = [bbatoms.dihedral.value() for bbatoms in self._phi_sel+self._psi_sel]
        dihe_df = pd.DataFrame(np.array(angs).T, columns=["phi","psi"], index=[frameid])
        self.framedata = self.framedata.append(dihe_df)
        return True

    def _update_selections(self):
        self._phi_sel = []
        self._psi_sel = []
        for resnum in self.resnums:
            temp_phi_sel = self.u.selectAtoms("protein and resid " + str(resnum-1) + " and name " + self.atoms[0][0])
            temp_phi_sel += self.u.selectAtoms("protein and resid " + str(resnum) + " and name " + self.atoms[0][1])
            temp_phi_sel += self.u.selectAtoms("protein and resid " + str(resnum) + " and name " + self.atoms[0][2])
            temp_phi_sel += self.u.selectAtoms("protein and resid " + str(resnum) + " and name " + self.atoms[0][3])
            temp_psi_sel = self.u.selectAtoms("protein and resid " + str(resnum) + " and name " + self.atoms[1][0])
            temp_psi_sel += self.u.selectAtoms("protein and resid " + str(resnum) + " and name " + self.atoms[1][1])
            temp_psi_sel += self.u.selectAtoms("protein and resid " + str(resnum) + " and name " + self.atoms[1][2])
            temp_psi_sel += self.u.selectAtoms("protein and resid " + str(resnum+1) + " and name " + self.atoms[1][3])

            if len(temp_phi_sel) != 4:
                raise Exception("PhiPsiAngle: atoms %s not found around target resid" % self.atoms[0])
            if len(temp_psi_sel) != 4:
                raise Exception("PhiPsiAngle: atoms %s not found around target resid" % self.atoms[1])

            self._phi_sel.append(temp_phi_sel)
            self._psi_sel.append(temp_psi_sel)

class Chi1Chi2Angle(PerFrameAnalysis):

    def __init__(self, resnums, atoms=[["N","CA","CB","CG"],["CA","CB","CG","CD"]]):
        """Calculates the chi1/chi2 angles of a sidechain based on an
           input selection of a residue.

        :Arguments:
          *resnums*
            A list of residue numbers within the protein

        """

        self.resnums = resnums
        self.atoms = atoms

    def process(self, frame, frameid):
        """ Process a single trajectory frame """
        angs = [bbatoms.dihedral.value() for bbatoms in self._chi1_sel+self._chi2_sel]
        dihe_df = pd.DataFrame(np.array(angs).T, columns=["chi1","chi2"], index=[frameid])
        self.framedata = self.framedata.append(dihe_df)

    def _update_selections(self):
        self._chi1_sel = []
        self._chi2_sel = []
        for resnum in self.resnums:
            temp_chi1_sel = self.u.selectAtoms("protein and resid " + str(resnum-1) + " and name " + self.atoms[0][0])
            temp_chi1_sel += self.u.selectAtoms("protein and resid " + str(resnum) + " and name " + self.atoms[0][1])
            temp_chi1_sel += self.u.selectAtoms("protein and resid " + str(resnum) + " and name " + self.atoms[0][2])
            temp_chi1_sel += self.u.selectAtoms("protein and resid " + str(resnum) + " and name " + self.atoms[0][3])
            temp_chi2_sel = self.u.selectAtoms("protein and resid " + str(resnum) + " and name " + self.atoms[1][0])
            temp_chi2_sel += self.u.selectAtoms("protein and resid " + str(resnum) + " and name " + self.atoms[1][1])
            temp_chi2_sel += self.u.selectAtoms("protein and resid " + str(resnum) + " and name " + self.atoms[1][2])
            temp_chi2_sel += self.u.selectAtoms("protein and resid " + str(resnum+1) + " and name " + self.atoms[1][3])

            if len(temp_chi1_sel) != 4:
                raise Exception("Chi1Chi2Angle: atoms %s not found around target resid" % self.atoms[0])
            if len(temp_chi2_sel) != 4:
                raise Exception("Chi1Chi2Angle: atoms %s not found around target resid" % self.atoms[1])

            self._chi1_sel.append(temp_phi_sel)
            self._chi2_sel.append(temp_psi_sel)
