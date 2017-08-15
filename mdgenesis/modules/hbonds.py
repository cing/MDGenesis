import pandas as pd
import MDAnalysis.analysis.hbonds
from analysis import AllAtOnceAnalysis

class HBonds(AllAtOnceAnalysis):

    def __init__(self, selection_str1="protein", selection_str2="protein"):
        """ Hydrogen bonds are computed for all frames, and we return
            the frequency per type for that trajectory. Using default dist/angle
            for CHARMM atom names.
        """
        self._selection_str1 = selection_str1
        self._selection_str2 = selection_str2

    def results(self):
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(self.u, self._selection_str1, self._selection_str2,
                                                            distance=3.0, angle=120.0, step=2, start=None, stop=None)
        h.run()
        h.generate_table()
        df = pd.DataFrame(h.count_by_type())
        df["donor_resnm"] = df["donor_resnm"].astype(str)
        df["donor_heavy_atom"] = df["donor_heavy_atom"].astype(str)
        df["donor_atom"] = df["donor_atom"].astype(str)
        df["acceptor_resnm"] = df["acceptor_resnm"].astype(str)
        df["acceptor_atom"] = df["acceptor_atom"].astype(str)
        return df
