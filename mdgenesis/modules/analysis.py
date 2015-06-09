"""
Base class for Analysis, the components of BatchAnalysis.

"""

import numpy as np
import pandas as pd

class PerFrameAnalysis(object):
    """ Base class for analysis that runs independently on each step
    as analyzed by the process method, returning framedata each time.
    Most analysis done using MDAnalysis should use this approach.
    """

    def __init__(self):
        """Base class for all timeseries based analysis. It stores
        a datafile called timeseries, that's about it. Subclasses
        will define the processing method.

        The timeseries accessible as the attribute :attr:`TSAnalysis.framedata`.
        """

        pass

    def prepare(self, trj, u=None, ref=None, intdata=pd.DataFrame(),
                framedata=pd.DataFrame()):
        """ Prepares the analysis routine and loads intermediate data
            if it exists. trj must be an iterable made outside of MDGenesis
        """

        self._loadcheckpoint(framedata, intdata)

        self.trj = trj
        self.u = u
        self.ref = ref

        self._update_selections()

    def process(self, frame):
        """ Process a single trajectory frame """
        self.framedata.append(True)

    def _update_selections(self):
        pass

    def _loadcheckpoint(self, framedata, intdata):
        self.framedata = framedata
        self.intdata = intdata

    def results(self):
        """ Returns an array of your analysis """
        return self.framedata

    def intresults(self):
        """ Returns an array of intermediate data """
        return self.intdata

class AllAtOnceAnalysis(object):
    """ Base class for analysis that runs all at once without
    explicit iteration over the frames of the set. Most analysis done using
    MDTraj should use this approach. These modules are not fault-tolerant.
    """

    def __init__(self):
        """Base class for all-at-once analysis. It stores
        a datafile which may be a timeseries, that's about it. Subclasses
        will define the processing method.

        The timeseries accessible as the attribute :attr:`TSAnalysis.framedata`.
        """

        pass

    def prepare(self, trj, u=None, ref=None):
        """ Prepares the analysis routine and loads intermediate data
            if it exists. trj must be an iterable made outside of MDGenesis """

        self.trj = trj
        self.u = u
        self.ref = ref
        self._update_selections()
        self.framedata = []  # final result

    def process(self, frame):
        raise NotImplementedError()

    def _update_selections(self):
        pass

    def results(self):
        """ Returns an array of your analysis """
        return self.framedata

    def intresults(self):
        raise NotImplementedError()
