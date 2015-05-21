==========================================================================
MDGenesis: an asynchronus analysis manager for molecular dynamics data
==========================================================================

Molecular dynamicists study the dynamic fluctuations of atoms, molecules,
and proteins through comptuational experiments. First they build models,
then they simulate them, then they extract data, and then they perform analysis,
ultimately leading to figures and tables that make up mind-blowing papers.
Extraction of data from simulation trajectories is typically performed with
analysis software like MDAnalysis, LOOS, VMD, MDTraj, or simulation package
specific utilities. Despite the existence of a wide selection of analysis
packages, many types of analysis scripts are routinely performed.

It is the aim of this package to provide a higher level organizational framework
for molecular dynamics analysis that encompasses various simulation packages
and stores data in a unified format. It is not the aim of this package to
provide yet another molecular dynamics analysis library (YAMDAP), but rather
to provide a suggested "best practice" when performing analysis that stands to
improve consistency and reproducibility with respect to the execution of
analysis. But let's be honest here, Below the surface, MDGenesis is nothing more
than a glorified wrapper around `MDSynthesis<https://github.com/Becksteinlab/MDSynthesis>`__ 
(but it might be something more someday!).

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

The principle behind this package and several analysis routines are taken
directly from the library `Batcha<https://github.com/davecap/batcha>`__.
