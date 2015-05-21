"""
Base class for Analysis, the components of BatchAnalysis.

"""

import numpy

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

    def run(self, trj):
        """ Analyze trajectory and produce timeseries. """
        self.framedata = []
        self.prepare(trj=trj, ref=ref)
        for ts in self.u.trajectory:
            logger.debug("Analyzing frame %d" % ts.frame)
            self.process(ts.frame)
        return self.framedata

    def prepare(self, trj=None, ref=None):
        """ Prepare the trajectory (trj is a Universe object). No reference object is needed. """
        self.u = trj
        self.ref = ref
        #self.u.trajectory.rewind()
        self._update_selections()  # Is this EVER needed?
        self.framedata = []  # final result

    def results(self):
        """ Returns an array containing the total count of nearby atoms """
        return numpy.array(self.framedata)

    def process(self, frame):
        """ Process a single trajectory frame """
        self.framedata.append(True)

    def _update_selections(self):
        pass

class AllAtOnceAnalysis(object):
    """ Base class for analysis that runs all at once without
    explicit iteration over the frames of the set. Most analysis done using
    MDTraj should use this approach.
    """

    def __init__(self):
        """Base class for all-at-once analysis. It stores
        a datafile which may be a timeseries, that's about it. Subclasses
        will define the processing method.

        The timeseries accessible as the attribute :attr:`TSAnalysis.framedata`.
        """

        pass

    def run(self, trj, ref=None):
        """ Analyze trajectory and produce timeseries. """
        self.framedata = []
        self.prepare(trj=trj, ref=ref)
        return self.results()

    def prepare(self, trj=None, ref=None, start=0, stop=-1):
        """ Prepare the trajectory (trj is a Universe object). No reference object is needed. """
        if stop != -1:
            self.u = trj[start:stop]
        else:
            self.u = traj[start:]
        self.ref = ref
        #self.u.trajectory.rewind()  # Is this EVER needed?
        self._update_selections()
        self.framedata = []  # final result

    def process(self, frame):
        raise NotImplementedError()

    def results(self):
        """ Returns an array containing the total count of nearby atoms """
        return numpy.array(self.framedata)

    def _update_selections(self):
        pass