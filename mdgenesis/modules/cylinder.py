# Cylinder analysis class
import numpy as np
import pandas as pd
from numpy.linalg import norm

from analysis import PerFrameAnalysis

class CylinderTiltHistogram(PerFrameAnalysis):

    def __init__(self, topsel, bottomsel, solutesel, refsel="protein",
                 radius=5, extension=0, soluteaxis=2,
                 histmin=-50, histmax=50, histbins=200):
        """Compute histogram of solute molecules coordinates (X, Y, or Z) within
        a cylinder defined by a top and bottom selection and a cylinder radius.
        This cylinder can tilt and even change height during simulation if
        the topsel and bottomsel were to fluctuate, so make sure your histogram
        range is larger than you expect. Note that in this case, the edges of
        the histogram may be undersampled so consider omitting them.

        :Arguments:
          *topsel*
            Top of cylinder selection string
          *bottomsel*
            Bottom of cylinder selection string
          *solutesel*
            Solute selection string
          *refsel*
            cylinder C.O.M. string
          *radius*
            radius of cylinder
          *extension*
            extension of cylinder
          *soluteaxis*
            coordination of solute molecules to histogram
          *histmin*
            minimum value of your box, relative to refsel
          *histmax*
            maximum value of your box, relative to refsel
          *histbins*
            number of histogram bins along the solute axis of the box

        """

        self.topsel = topsel
        self.bottomsel = bottomsel
        self.solutesel = solutesel
        self.refsel = refsel
        self.r = radius
        self.extension = extension
        self.saxis = soluteaxis

        self.histbins = histbins
        self.histmin = histmin
        self.histmax = histmax

        #self._all_edges = np.linspace(histmin, histmax, num=histbins)

    # We don't need to initialize framedata because it's only made
    # when results() is called.
    def _loadcheckpoint(self, framedata, intdata):
        self.framedata = framedata
        if intdata.empty:
            self.intdata = pd.DataFrame(np.zeros([self.histbins,2],
                                                 dtype=np.int64),
                                        columns=["bincount",
                                                 #"boolcount",
                                                 "total_frames"])
        else:
            self.intdata = intdata

    def process(self, frame, frameid):
        """ Process a single trajectory frame """
        self._update_selections()

        # distance to the helix axis defined by top_com and bottom_com
        d = norm(np.cross(self._scoord - self._top_com,
                          self._scoord - self._bottom_com),
                 axis=1)/self._height

        # boolean mask for distances below the radius, hist'd on solute-axis
        h = np.histogram(self._scoord[d <= self.r][:,self.saxis] - self._ref_com[self.saxis],
                         range=[self.histmin, self.histmax],
                         bins=self.histbins, normed=False)[0]

        # optionally output the resids of solute in cylinder (useful for VMD)
        #print self._sids[dist <= self.radius]

        # Yep, I make three entirely new dataframes for each frame!
        self.intdata.ix[0, "total_frames"] += 1
        #self.intdata["boolcount"] += pd.Series(h) > 0
        self.intdata["bincount"] += pd.Series(h)

        # If you do not return True, then the frame will be re-analyzed next
        # time, but heads up, that could mess up your intermediate data!
        return True

    def results(self):
        frames_processed = float(self.intdata["total_frames"][0])
        edges = pd.DataFrame(np.linspace(self.histmin, self.histmax,
                                         num=self.histbins+1)[:self.histbins],
                             columns=["edges"])

        if frames_processed > 0:
            return pd.concat([self.intdata["bincount"]/frames_processed,
                              edges], axis=1)
        else:
            return pd.concat([self.intdata["bincount"], edges], axis=1)

    def _update_selections(self):
        self._top_com = self.u.select_atoms(self.topsel).center_of_mass()
        self._bottom_com = self.u.select_atoms(self.bottomsel).center_of_mass()
        self._ref_com = self.u.select_atoms(self.refsel).center_of_mass()

        self._midpoint = (self._top_com+self._bottom_com)/2.0
        self._height = np.linalg.norm(self._top_com-self._bottom_com)
        self._sradius = (self._height/2) + self.extension

        #self._sids = np.array([sol.resid for sol in self.u.selectAtoms(self.solutesel)])
        self._scoord = self.u.select_atoms(self.solutesel).coordinates()

class CylinderHistogram(PerFrameAnalysis):

    def __init__(self, solutesel, refsel="protein",
                 radius=5, extension=0, soluteaxis=2,
                 histmin=-50, histmax=50, histbins=200):
        """Compute histogram of solute molecules coordinates (X, Y, or Z) within
        a cylinder along a given axis.

        :Arguments:
          *solutesel*
            Solute selection string
          *refsel*
            cylinder C.O.M. string
          *radius*
            radius of cylinder
          *extension*
            extension of cylinder centered at the C.O.M. value
          *soluteaxis*
            coordination of solute molecules to histogram
          *histmin*
            minimum value of your box, relative to refsel
          *histmax*
            maximum value of your box, relative to refsel
          *histbins*
            number of histogram bins along the solute axis of the box

        """

        self.solutesel = solutesel
        self.refsel = refsel
        self.r2 = radius**2
        self.extension = extension
        self.saxis = soluteaxis
        self.osaxis = list(set(range(3))-set([self.saxis]))

        self.histbins = histbins
        self.histmin = histmin
        self.histmax = histmax

        #self._all_edges = np.linspace(histmin, histmax, num=histbins)

    # We don't need to initialize framedata because it's only made
    # when results() is called.
    def _loadcheckpoint(self, framedata, intdata):
        self.framedata = framedata
        if intdata.empty:
            self.intdata = pd.DataFrame(np.zeros([self.histbins,2],
                                                 dtype=np.int64),
                                        columns=["bincount",
                                                 #"boolcount",
                                                 "total_frames"])
        else:
            self.intdata = intdata

    def process(self, frame, frameid):
        """ Process a single trajectory frame """
        self._update_selections()

        radius_bool = np.sum(((self._scoord - self._ref_com)**2)[:,self.osaxis],axis=1) <= self.r2
        height_bool1 = (self._scoord - self._ref_com)[:,self.saxis] < self._height*0.5
        height_bool2 = (self._scoord - self._ref_com)[:,self.saxis] > -self._height*0.5
        cylinder_bool = radius_bool & height_bool1 & height_bool2

        # boolean mask for distances below the radius, hist'd on solute-axis
        h = np.histogram(self._scoord[cylinder_bool][:,self.saxis] - self._ref_com[self.saxis],
                         range=[self.histmin, self.histmax],
                         bins=self.histbins, normed=False)[0]

        # optionally output the resids of solute in cylinder (useful for VMD)
        #print self._sids[cylinder_bool]

        # Yep, I make three entirely new dataframes for each frame!
        self.intdata.ix[0, "total_frames"] += 1
        #self.intdata["boolcount"] += pd.Series(h) > 0
        self.intdata["bincount"] += pd.Series(h)

        # If you do not return True, then the frame will be re-analyzed next
        # time, but heads up, that could mess up your intermediate data!
        return True

    def results(self):
        frames_processed = float(self.intdata["total_frames"][0])
        edges = pd.DataFrame(np.linspace(self.histmin, self.histmax,
                                         num=self.histbins+1)[:self.histbins],
                             columns=["edges"])

        if frames_processed > 0:
            return pd.concat([self.intdata["bincount"]/frames_processed,
                              edges], axis=1)
        else:
            return pd.concat([self.intdata["bincount"], edges], axis=1)

    def _update_selections(self):
        self._ref_com = self.u.select_atoms(self.refsel).center_of_mass()
        self._height = self.extension

        #self._sids = np.array([sol.resid for sol in self.u.select_atoms(self.solutesel)])
        self._scoord = self.u.select_atoms(self.solutesel).coordinates()

class ConditionalCylinderHistogram(CylinderHistogram):

    def __init__(self, solutesel, refsel="protein",
                 radius=5, extension=0, soluteaxis=2,
                 histmin=-50, histmax=50, histbins=200,
                 conditional_sel=None, conditional_bool=None):
        """Compute histogram of solute molecules coordinates (X, Y, or Z) within
        a cylinder along a given axis based on a boolean conditional. Useful
        for examining hydration based on sidechain conformation, etc.

        :Arguments:
          *solutesel*
            Solute selection string
          *refsel*
            cylinder C.O.M. string
          *radius*
            radius of cylinder
          *extension*
            extension of cylinder centered at the C.O.M. value
          *soluteaxis*
            coordination of solute molecules to histogram
          *histmin*
            minimum value of your box, relative to refsel
          *histmax*
            maximum value of your box, relative to refsel
          *histbins*
            number of histogram bins along the solute axis of the box
          *conditional_sel*
            string for eval() to select/compute quantity of interest
          *conditional_bool*
            string for eval() that returns True/False (use vars in conditional_sel)

        """

        self.solutesel = solutesel
        self.refsel = refsel
        self.r2 = radius**2
        self.extension = extension
        self.saxis = soluteaxis
        self.osaxis = list(set(range(3))-set([self.saxis]))

        self.histbins = histbins
        self.histmin = histmin
        self.histmax = histmax

        self.conditional_sel = conditional_sel
        self.conditional_bool = conditional_bool

        #self._all_edges = np.linspace(histmin, histmax, num=histbins)

    # We don't need to initialize framedata because it's only made
    # when results() is called.
    def _loadcheckpoint(self, framedata, intdata):
        self.framedata = framedata
        if intdata.empty:
            self.intdata = pd.DataFrame(np.zeros([self.histbins,3],
                                                 dtype=np.int64),
                                        columns=["bincount_true",
                                                 "bincount_false",
                                                 "total_frames"])
        else:
            self.intdata = intdata

    def process(self, frame, frameid):
        """ Process a single trajectory frame """
        self._update_selections()

        radius_bool = np.sum(((self._scoord - self._ref_com)**2)[:,self.osaxis],axis=1) <= self.r2
        height_bool1 = (self._scoord - self._ref_com)[:,self.saxis] < self._height*0.5
        height_bool2 = (self._scoord - self._ref_com)[:,self.saxis] > -self._height*0.5
        cylinder_bool = radius_bool & height_bool1 & height_bool2

        # boolean mask for distances below the radius, hist'd on solute-axis
        h = np.histogram(self._scoord[cylinder_bool][:,self.saxis] - self._ref_com[self.saxis],
                         range=[self.histmin, self.histmax],
                         bins=self.histbins, normed=False)[0]

        # Yep, I make three entirely new dataframes for each frame!
        self.intdata.ix[0, "total_frames"] += 1

        if (self.conditional_sel != None) & (self.conditional_bool != None):
            eval(self.conditional_sel)
            if eval(self.conditional_bool):
                self.intdata["bincount_true"] += pd.Series(h)
            else:
                self.intdata["bincount_false"] += pd.Series(h)
        else:
            self.intdata["bincount_true"] += pd.Series(h)

        # If you do not return True, then the frame will be re-analyzed next
        # time, but heads up, that could mess up your intermediate data!
        return True

    def results(self):
        frames_processed = float(self.intdata["total_frames"][0])
        edges = pd.DataFrame(np.linspace(self.histmin, self.histmax,
                                         num=self.histbins+1)[:self.histbins],
                             columns=["edges"])

        if frames_processed > 0:
            return pd.concat([self.intdata["bincount_true"]/frames_processed,
                              self.intdata["bincount_false"]/frames_processed,
                              edges], axis=1)
        else:
            return pd.concat([self.intdata["bincount_true"],
                              self.intdata["bincount_false"], edges], axis=1)

class CylinderCount(PerFrameAnalysis):

    def __init__(self, solutesel, minval=-10, maxval=10,
                 refsel="protein",
                 radius=5, extension=0, soluteaxis=2):
        """Compute a count of solute molecules coordinates within
        a cylinder along a given axis.

        :Arguments:
          *solutesel*
            Solute selection string
          *minval*
            minimum value to consider count, rel. to refsel
          *maxval*
            maximum value to consider count, rel. to refsel
          *refsel*
            cylinder C.O.M. string
          *radius*
            radius of cylinder
          *extension*
            extension of cylinder centered at the C.O.M. value
          *soluteaxis*
            coordination of solute molecules to histogram

        """

        self.solutesel = solutesel
        self.minval = minval
        self.maxval = maxval
        self.refsel = refsel
        self.r2 = radius**2
        self.extension = extension
        self.saxis = soluteaxis
        self.osaxis = list(set(range(3))-set([self.saxis]))

    def _loadcheckpoint(self, framedata, intdata):
        self.intdata = intdata
        if framedata.empty:
            self.framedata = pd.DataFrame(columns=["count"])
        else:
            self.framedata = framedata

    def process(self, frame, frameid):
        """ Process a single trajectory frame """
        self._update_selections()

        radius_bool = np.sum(((self._scoord - self._ref_com)**2)[:,self.osaxis],axis=1) <= self.r2
        height_bool1 = (self._scoord - self._ref_com)[:,self.saxis] < self.maxval
        height_bool2 = (self._scoord - self._ref_com)[:,self.saxis] > self.minval
        cylinder_bool = radius_bool & height_bool1 & height_bool2

        c = np.sum(cylinder_bool)
        c_df = pd.DataFrame(c, columns=["count"], index=[frameid])
        self.framedata = self.framedata.append(c_df)
        return True

    def _update_selections(self):
        self._ref_com = self.u.select_atoms(self.refsel).center_of_mass()

        #self._sids = np.array([sol.resid for sol in self.u.select_atoms(self.solutesel)])
        self._scoord = self.u.select_atoms(self.solutesel).coordinates()
