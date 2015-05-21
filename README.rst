==========================================================================
MDGenesis: a batch molecular dynamics analysis executor/wrapper
==========================================================================

Molecular dynamicists study the dynamic fluctuations of atoms, molecules,
and proteins through comptuational experiments. They build models,
they simulate them, they extract quantitative data,
and finally they make figures that result in mind-blowing papers.
Extracting data from simulation data sets is typically performed with
analysis software like `MDAnalysis <https://github.com/MDAnalysis/mdanalysis>`_, 
`LOOS <http://loos.sourceforge.net/>`_, `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_,
`MDTraj <https://github.com/mdtraj/mdtraj>`_, or other simulation package utilities.

It is the aim of this package to provide a higher level organizational framework
for executing molecular dynamics analysis that supports many analysis packages
(assuming they have a python interface) and stores resulting data in a unified
manner. It is not the aim of this package to be "yet another molecular dynamics
analysis library" (YAMDAP), but rather to provide a wrapper around any analysis
performed in post-processing (after the trajectory data has been output to disk).

Benefits of Wrapping your Analysis
==================================

MDGenesis is an analysis wrapper that aims to provide a higher level
interface for executing analysis modules of all varieties. While it currently
offers no benefits other than organizational ones, future functionality may
include the following:

* Parallelism using analysis batches (analysis #1 on CPU1, analysis #2 on CPU2)
* Parallelism using chunked analysis (frames 0-10 on CPU1, frames 10-20 on CPU2)
* Detection and analysis of new data (+10 frames in source, append analysis)
* Fault-tolerant analysis (recovery from crashed analysis)
* Auto-population of analysis across datasets (new analysis, update for all data)
* Persistent analysis (daemon-mode, MDGenesis sniffs out and analyzes new data)

Shout Outs
==========

Let's be honest here, below the thin veil of README jargon and empty promises,
MDGenesis is nothing more than a glorified wrapper around 
`MDSynthesis <https://github.com/Becksteinlab/MDSynthesis>`_.
Thanks to David Dotson for creating it!

The principle behind this package and several analysis routines are taken
directly from the library `Batcha <https://github.com/davecap/batcha>`__.
