import mdtraj as md
import pandas as pd
from analysis import AllAtOnceAnalysis

class DSSP(AllAtOnceAnalysis):

    def __init__(self, simplify=True):
        """ DSSP is computed for all atoms, but specify individual
            to be returned with selection.
        """
        self.simplify = simplify

    def results(self):
        self.framedata = md.compute_dssp(self.trj, self.simplify)[:,self._selection]
        resid_list = [res.resSeq for res in self.trj.topology.residues if res.is_protein]
        return pd.DataFrame(self.framedata, columns=resid_list)

    def _update_selections(self):
        self._selection = [res.index for res in self.trj.topology.residues if res.is_protein]
