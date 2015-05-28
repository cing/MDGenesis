"""
Base class for Analysis, the components of BatchAnalysis.

"""

import numpy as np

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

    def run(self, trj, u=None, frames_processed=0, intdata=None):
        """ Analyze trajectory and produce timeseries. Use when you don't
            care about saving intermediate data during long calculations
            or only analyzing a subset of frames.
        """
        self.framedata = []
        self.prepare(trj=trj, u=u, ref=ref,
                     frames_processed=frames_processed, intdata=intdata)
        for frame in self.trj:
            logger.debug("Analyzing frame %d" % ts.frame)
            self.process(frame)
        return self.framedata

    def prepare(self, trj, u=None, start=0, stop=-1,
                ref=None, frames_processed=0, intdata=None):
        """ Prepares the analysis routine and loads intermediate data
            if it exists. trj must be an iterable made outside of MDGenesis """

        self._loadcheckpoint(frames_processed, intdata)

        if stop != -1:
            self.trj = trj[start:stop]
        else:
            self.trj = trj[start:]
        self.u = u
        self.ref = ref

        self.framedata = []  # final result

    def process(self, frame):
        """ Process a single trajectory frame """
        self.framedata.append(True)

    def _update_selections(self):
        pass

    def _loadcheckpoint(self, frames_processed, intdata):
        self.frames_processed = frames_processed
        self.intdata = intdata

    def results(self):
        """ Returns an array of your analysis """
        return np.array(self.framedata)

    def intresults(self):
        """ Returns an array of intermediate data """
        return np.array(self.intdata)

    def framecount(self):
        """ Returns the number of frames processed """
        return self.frames_processed

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

    def run(self, trj, u=None, ref=None):
        """ Analyze trajectory and return results array. """
        self.framedata = []
        self.prepare(trj=trj, u=u, ref=ref)
        return self.results()

    def prepare(self, trj, u=None, ref=None, start=0, stop=-1):
        """ Prepares the analysis routine and loads intermediate data
            if it exists. trj must be an iterable made outside of MDGenesis """
        if stop != -1:
            self.trj = trj[start:stop]
        else:
            self.trj = trj[start:]
        self.u = u
        self.ref = ref
        #self.u.trajectory.rewind()  # Is this EVER needed?
        self._update_selections()
        self.framedata = []  # final result

    def process(self, frame):
        raise NotImplementedError()

    def _update_selections(self):
        pass

    def results(self):
        """ Returns an array of your analysis """
        return np.array(self.framedata)

    def intresults(self):
        raise NotImplementedError()
