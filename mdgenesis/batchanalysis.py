from collections import namedtuple, OrderedDict
from MDAnalysis import collection
import numpy as np
import pandas as pd
import mdsynthesis as mds

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

    def run(self, trj, u=None, ref=None, start=0, stop=-1, skip=1):
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

    def run_sequential_analysis(self, start, stop, skip=1):

        existing_data = {}
        existing_start_stops = {}

        for path, analysis in self._sequential.items():
            print " Preparing %s" % path
            existing_data[path] = self.sim.data.retrieve(path)

            # Initialize an analysis module from the start or resume it at a
            # later time if a checkpoint or synchronization flag has been set.
            if (analysis.sync and (path+"/analysis_stats" in self.sim.data)):
                cstart, cstop, cskip, ccheckpoint, cframecount = self.sim.data[path+"/analysis_stats"]
                print "Start/Stop/Skip/FC in analysis_stats: ", cstart, cstop, cskip, cframecount

                # Have there been changes in this execution versus the last time?
                if (cstart < start or cskip != skip or cstop > stop or
                    ccheckpoint != analysis.checkpoint):
                    print "Detected difference in start/stop/skip/checkpoint"
                    print "Gotta restart from t=0, hope you don't mind!"
                    analysis.func.prepare(self._trj, u=self._u, ref=self._ref,
                                          start=start, stop=stop)
                    existing_start_stops[path] = (start,stop)

                else:
                    # Load intermediate data and framedata if it exists!
                    if (path in self.sim.data):
                        framedata_arg = self.sim.data[path]
                    else:
                        framedata_arg = pd.DataFrame()
                    if (path+"/checkpoint" in self.sim.data):
                        intdata_arg = self.sim.data[path+"/checkpoint"]
                    else:
                        intdata_arg = pd.DataFrame()

                    print "Found checkpoint: ", self.sim.data[path+"/checkpoint"]
                    print "Starting at: ", start+int(cframecount)
                    #print "Found framedata: ", self.sim.data[path]
                    analysis.func.prepare(self._trj, u=self._u, ref=self._ref,
                                   start=start+int(cframecount)+1, stop=stop,
                                   framedata=framedata_arg,
                                   intdata=intdata_arg,
                                   frames_processed=int(cframecount))

                    existing_start_stops[path] = (start+int(cframecount)+1,stop)
            else:
                analysis.func.prepare(self._trj, u=self._u, ref=self._ref, start=start, stop=stop)
                existing_start_stops[path] = (start,stop)

        # We have to reslice the trajectory in case we are checkpointing/syncing
        # a previous run that had a different start/skip/stop time.
        for path, analysis in self._sequential.items():
            start, stop = existing_start_stops[path]
            print "Start, stop are: ", start, stop
            if stop != -1:
                frames = self._trj[start:stop]
            else:
                frames = self._trj[start:]

            if path+"/analysis_stats" not in self.sim.data:
                self.sim.data[path+"/analysis_stats"] = np.array([start, stop, skip,
                                                                  analysis.checkpoint,
                                                                  analysis.func.framecount()])

            #print " Processing %d frames..." % frames.numframes
            for i, f in enumerate(frames):
                current_frame = i + start
                print i, i+start, analysis.checkpoint, current_frame % analysis.checkpoint == 0
                # TODO: throw error if trajectory actually doesn't support checkpoints
                if analysis.checkpoint != 0:
                    if current_frame % analysis.checkpoint == 0 and i > 0:

                        # This stores how many frames we've analyzed, critical!
                        self.sim.data[path+"/analysis_stats"] = np.array([start, stop, skip, analysis.checkpoint,
                                                                 analysis.func.framecount()])

                        # The analysis.func routine may or may not update this.
                        results = analysis.func.results()
                        if len(results) > 0:
                            #print "Storing Results at", path, analysis.func.results(), analysis.func.framecount()
                            print "Storing Results at", path, analysis.func.framecount()
                            self.sim.data[path] = results

                        # The analysis.func routine may or may not update this.
                        intresults = analysis.func.intresults()
                        if len(intresults) > 0:
                            #print "Storing Checkpoint at", path+"/checkpoint", analysis.func.intresults(), analysis.func.framecount()
                            print "Storing Checkpoint at", path+"/checkpoint", analysis.func.framecount()
                            self.sim.data[path+"/checkpoint"] = intresults

                analysis.func.process(f)

        print " done."

        print " Loading result data..."
        for path, analysis in self._sequential.items():
            # If there is no data or synchronize is False
            if (existing_data[path] is None) or not analysis.sync:
                self.sim.data.add(path, analysis.func.results())
            else:
                print "Nothing done, file is complete!"

    def run_allatonce_analysis(self, start, stop, skip=1):

       existing_data = {}
       for path, analysis in self._allatonce.items():
           print " Preparing %s" % path
           analysis.func.prepare(self._trj, u=self._u, ref=self._ref, start=start, stop=stop)
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

