import mdtraj as md
import pandas as pd
from analysis import AllAtOnceAnalysis
from MDAnalysis.analysis.hole import HOLEtraj

class PoreRadius(AllAtOnceAnalysis):

    def __init__(self, executable="/Users/cing/Code/hole2/exe/hole2"):
        """ HOLE is computed for all atoms, with the profile returned
            as per usual.
        """
        self.executable = executable

    def results(self):
        H = MDAnalysis.analysis.hole.HOLEtraj(self.u, executable=self.executable,
                                              cvect=[0,0,1], cpoint=True)
        H.run()
        return pd.DataFrame(H.profiles)
