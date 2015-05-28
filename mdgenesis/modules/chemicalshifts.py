import mdtraj as md
from analysis import AllAtOnceAnalysis

class ShiftX2(AllAtOnceAnalysis):

    def __init__(self, ph=7.0, temp=300, skip=100, atom_selection=None):
        """ DSSP is computed for all atoms, but specify individual
            to be returned with selection.
        """

        self._ph = ph
        self._temp = temp
        self._skip = skip
        self._atom_selection = atom_selection

    def results(self):
        temp_results = md.nmr.chemical_shifts_shiftx2(self.traj[::self._skip],
                                                      pH=self._ph,
                                                      temperature=self._temp)

        if self._atom_selection is not None:
            i = temp_results.index.get_level_values('name').isin(self._atom_selection)
            return temp_results[i]
        else:
            return temp_results
