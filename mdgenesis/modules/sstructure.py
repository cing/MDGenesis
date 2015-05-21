import mdtraj as md
from analysis import AllAtOnceAnalysis

class DSSP(AllAtOnceAnalysis):

    def __init__(self):
        """ DSSP is computed for all atoms, but specify individual
            to be returned with selection.
        """

    def results(self):
        self.framedata = md.compute_dssp(self.u)[:,self._selection]
        return self.framedata

    def _update_selections(self):
        self._selection = [res.index for res in self.u.topology.residues if res.is_protein]