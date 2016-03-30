import numpy as np
import mdtraj as md
import pandas as pd
from analysis import AllAtOnceAnalysis
from MDAnalysis.analysis.hole import HOLEtraj

class PoreRadius(AllAtOnceAnalysis):

    def __init__(self, step=1, executable="/Users/cing/Code/hole2/exe/hole"):
        """ HOLE is computed for all atoms, with the profile returned
            as per usual.
        """
        self.step = step
        self.executable = executable

    def results(self):
        H = HOLEtraj(self.u, step=self.step, executable=self.executable,
                     cvect=[0,0,1], cpoint=True)
        H.run()
        # For all frames, make a Series with average rxncoord
        # and mean across all time points.
        temp_series=[pd.Series(H.profiles[x]["radius"],
                               index=np.round(H.profiles[x]["rxncoord"],1)) for x in H.profiles]

        return pd.DataFrame(temp_series).T
