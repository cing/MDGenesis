from MDAnalysis import collection
#import numpy as np
#import pandas as pd
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
        self._sequential = {}
        self._dcd_timeseries = {}
        self._allatonce = {}

    def add_dcd_timeseries(self, path, timeseries):
        self._dcd_timeseries[path] = (timeseries)

    def add_to_sequence(self, path, processor, synchronize=False):
        self._sequential[path] = (processor, synchronize)

    def add_allatonce(self, path, processor, synchronize=False):
        self._allatonce[path] = (processor, synchronize)

    def run(self, trj, u=None, ref=None, start=0, stop=-1, skip=1, checkpoint=0):
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
            self.run_sequential_analysis(start, stop, skip=skip, checkpoint=checkpoint)
            print "Done sequential analysis."

        if len(self._allatonce) > 0:
            print "Running all at once analyses..."
            self.run_allatonce_analysis(start, stop, skip=skip)
            print "Done all-at-once analysis."

    def run_dcd_timeseries(self, start, stop, skip=1):
        collection.clear()
        for path, tpl in self._dcd_timeseries.items():
            print " Adding timeseries: %s" % path
            collection.addTimeseries(tpl)

        print " Computing..."
        collection.compute(self._u, start=start, stop=stop, skip=skip)
        print " Done computing."

        print "Loading data..."
        for i, path in enumerate(self._dcd_timeseries.keys()):
            print " loading table %s with %d values..." % (path, len(collection[i][0]))
            self.sim.data.add(path, collection[i][0])

    def run_sequential_analysis(self, start, stop, skip=1, checkpoint=0):

        existing_data = {}
        for path, tpl in self._sequential.items():
            print " Preparing %s" % path
            existing_data[path] = self.sim.data.retrieve(path)

            # Prepare the analysis module from a checkpoint or for the first time
            if (path+"/checkpoint" in self.sim.data) and (path+"/framecount" in self.sim.data):
                print "Found checkpoint: ", self.sim.data[path+"/checkpoint"]
                print "Found frames: ", self.sim.data[path+"/framecount"]
                tpl[0].prepare(self._trj, u=self._u, ref=self._ref, start=start, stop=stop,
                               intdata=self.sim.data[path+"/checkpoint"],
                               frames_processed=self.sim.data[path+"/framecount"])
            else:
                tpl[0].prepare(self._trj, u=self._u, ref=self._ref, start=start, stop=stop)

        if stop != -1:
            frames = self._trj[start:stop]
        else:
            frames = self._trj[start:]

        #print " Processing %d frames..." % frames.numframes
        for i, f in enumerate(frames):
            # TODO: throw error if trajectory actually doesn't support checkpoints
            if checkpoint != 0:
                if i % checkpoint == 0 and i > 0:
                    print ".",
                    print path+"/checkpoint", tpl[0].intresults(), tpl[0].framecount()
                    self.sim.data[path+"/checkpoint"] = tpl[0].intresults()
                    self.sim.data[path+"/framecount"] = tpl[0].framecount()

            for path, tpl in self._sequential.items():
                # if tpl[1] is True, check if frame needs to be analyzed
                # otherwise, just analyze
                if tpl[1] and existing_data[path] is not None:
                    if i >= existing_data[path].shape[0]:
                        tpl[0].process(f)
                else:
                    tpl[0].process(f)
        print " done."

        print " Loading result data..."
        for path, tpl in self._sequential.items():
            # If there is no data or synchronize is False
            if (existing_data[path] is None) or not tpl[1]:
                self.sim.data.add(path, tpl[0].results())
            elif tpl[0].results().shape[0] < existing_data[path].shape[0]:
            #elif tpl[0].results().shape[0] < existing_data[path].shape[0]+1: #Debug line to always append
                self.sim.data.append(path, tpl[0].results())
            else:
                print "Nothing done, file is complete!"

    def run_allatonce_analysis(self, start, stop, skip=1):

       existing_data = {}
       for path, tpl in self._allatonce.items():
           print " Preparing %s" % path
           tpl[0].prepare(self._trj, u=self._u, ref=self._ref, start=start, stop=stop)
           existing_data[path] = self.sim.data.retrieve(path)

       print " Computing/Loading result data..."
       for path, tpl in self._allatonce.items():
           # If there is no data or synchronize is False
           if (existing_data[path] is None) or not tpl[1]:
               self.sim.data.add(path, tpl[0].results())
           elif tpl[0].results().shape[0] < existing_data[path].shape[0]:
           #elif tpl[0].results().shape[0] < existing_data[path].shape[0]+1: #Debug line to always append
               #self.sim.data.append(path, pd.DataFrame(tpl[0].results()))
               self.sim.data.append(path, tpl[0].results())
           else:
               print "Nothing done, file is complete!"

