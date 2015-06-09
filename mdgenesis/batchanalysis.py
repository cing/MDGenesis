from collections import namedtuple, OrderedDict
from MDAnalysis import collection
import numpy as np
import pandas as pd
import mdsynthesis as mds
import datetime as dt

class BatchAnalysis():
    '''
    A set of analyis routines of type "DCD Timeseries (MDAnalysis-only)",
    "Frame-by-Frame (Sequence)", or "All-at-once". BatchAnalysis
    must be instantiated for each trajectory file, which should have
    a unique filename. It will create a Sim object in MDSynthesis for
    that trajectory and append all analysis data to HDF5 files contained therein.
    Performing analysis on multiple trajectories is not yet supported.
    '''

    def __init__(self, simname, universe=None):

        self.Analysis = namedtuple("Analysis", "func sync checkpoint")

        self._simname = simname
        self.sim = self.open_or_create(universe)
        self.reset()

    def open_or_create(self, universe=None):
        """ 
        Create the MDSynthesis Sim instance, with or without an MDAnalysis
        universe list of ['path/to/topology', 'path/to/trajectory']
        """
        #TODO: Store categories during creation that correspond to the
        #      types of analysis performed (?)
        if universe is not None:
            return mds.Sim(self._simname, universe=universe)
        else:
            return mds.Sim(self._simname)

    def reset(self):
        # the following dicts store the actual analyses which are processed
        self._sequential = OrderedDict()
        self._dcd_timeseries = OrderedDict()
        self._allatonce = OrderedDict()

    def add_dcd_timeseries(self, path, func):
        self._dcd_timeseries[path] = self.Analysis(func=func)

    def add_to_sequence(self, path, func, sync=False, checkpoint=0):
        self._sequential[path] = self.Analysis(func=func, sync=sync,
                                               checkpoint=checkpoint)

    def add_allatonce(self, path, func, sync=False):
        self._allatonce[path] = self.Analysis(func=func, sync=sync)

    # TODO: There isn't a mechanism to avoid setting a stop value!
    #       In other words, then doesn't really work for a persistent process
    def run(self, trj, u=None, ref=None, start=0, stop=1000000, skip=1):
        '''This executes all loaded analysis types'''
        self._trj = trj
        self._u = u
        self._ref = ref

        if len(self._dcd_timeseries) > 0:
            print "Starting timeseries analysis..."
            self.run_dcd_timeseries(start, stop, skip=skip)
            print "Done timeseries analysis."

        if len(self._sequential) > 0:
            print "Running sequential analyses..."
            self.run_sequential_analysis(start, stop, skip=skip)
            print "Done sequential analysis."

        if len(self._allatonce) > 0:
            print "Running all at once analyses..."
            self.run_allatonce_analysis(start, stop, skip=skip)
            print "Done all-at-once analysis."

    def run_dcd_timeseries(self, start, stop, skip=1):

        collection.clear()
        for path, analysis in self._dcd_timeseries.items():
            print " Adding timeseries: %s" % path
            collection.addTimeseries(analysis.func)

        print " Computing..."
        collection.compute(self._u, start=start, stop=stop, skip=skip)
        print " Done computing."

        print "Loading data..."
        for i, path in enumerate(self._dcd_timeseries.keys()):
            print " loading table %s with %d values..." % (path, len(collection[i][0]))
            self.sim.data[path] = collection[i][0]
            self.sim.data[path+"/analysis_stats"] = np.array([start, stop, skip, -1, -1])

    def unanalyzed_frames(self, path, start, stop, skip):
        """ Does the sole job of returning a dataframe of unanalyzed frames """
        # This is the only place where we decide what frames are done
        possible_frames = set(np.arange(start, stop, skip))
        analyzed_frames = self.sim.data[path+'/analysis_log']
        #print path, start, stop, skip, possible_frames, analyzed_frames
        if analyzed_frames is not None:
            return possible_frames - set(analyzed_frames.index)
        else:
            return possible_frames

    # TODO: why load framedata?
    def read_stored_frames(self, path):
        intdata = self.sim.data[path+'/intermediate_data']
        framedata = self.sim.data[path]
        if intdata is None:
            intdata = pd.DataFrame()
        if framedata is None:
            framedata = pd.DataFrame()
        return intdata, framedata

    def run_sequential_analysis(self, start, stop=1000000, skip=1):

        frames_to_analyze = {}
        for path, analysis in self._sequential.items():
            print " Preparing %s" % path
            intdata, framedata = self.read_stored_frames(path)
            analysis.func.prepare(self._trj, u=self._u, ref=self._ref,
                                  intdata=intdata,
                                  framedata=framedata)

            if not analysis.sync:
                print " Sync not set, clearing %s/analysis_log" % path
                self.sim.data.remove(path+'/analysis_log')

        for path, analysis in self._sequential.items():
            print " Processing %s" % path
            # TODO: consider doing this more often to exploit persistence layer
            #       in a better way!
            frames_to_process = self.unanalyzed_frames(path, start, stop, skip)
            if len(frames_to_process) > 0:
                completed_frames = {}
                for i, f in enumerate(self._trj):
                    # Visit the analysis module for the heavy lifting
                    if i in frames_to_process:
                        if analysis.func.process(f, i):
                            completed_frames[i] = [dt.datetime.today()]

                    # Write the checkpoint (once you've completed enough frames)
                    if (len(completed_frames) > 0 and
                        len(completed_frames) % analysis.checkpoint == 0):
                        print " Reached checkpoint at frame: %i" % i
                        new_frame = pd.DataFrame.from_items(completed_frames.items(),
                                                            orient='index',
                                                            columns=['date_completed'])
                        self.sim.data.append(path+'/analysis_log', new_frame)
                        completed_frames = {}
                        self._write_state(path, analysis)

                # One final write for the final data
                print " Writing to %s" % path
                self._write_state(path, analysis)

    def _write_state(self, path, analysis):
        intresults = analysis.func.intresults()
        if not intresults.empty:
            self.sim.data[path+"/intermediate_data"] = intresults
        results = analysis.func.results()
        if not results.empty:
            self.sim.data[path] = results

    def run_allatonce_analysis(self, start, stop, skip=1):

        existing_data = {}
        for path, analysis in self._allatonce.items():
            print " Preparing %s" % path
            analysis.func.prepare(self._trj, u=self._u, ref=self._ref)
            existing_data[path] = self.sim.data.retrieve(path)
            self.sim.data[path+"/analysis_stats"] = np.array([start, stop, skip, -1, -1])

        print " Computing/Loading result data..."
        for path, analysis in self._allatonce.items():
            # If there is no data or synchronize is False
            if (existing_data[path] is None) or not analysis.sync:
                self.sim.data.add(path, analysis.func.results())
            elif analysis.func.results().shape[0] < existing_data[path].shape[0]:
                self.sim.data.append(path, analysis.func.results())
            else:
                print "Nothing done, file is complete!"
