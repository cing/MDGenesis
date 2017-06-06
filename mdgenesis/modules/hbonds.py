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
                                                            distance=3.0, angle=120.0)
        h.run()
        h.generate_table()
        return pd.DataFrame(h.count_by_type())
