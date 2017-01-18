# H2hh2bbtautauAnalysis

N.B. Basis forked from DESY CMSSW School 2016

Quick tutorial on how to run the analysis:

## Skimming

Modify datasets.cfg and use `./crabControl.py --config datasets.cfg --submit` to submit your skim jobs to the grid.
To check the status of your jobs, use `./crabControl.py --config datasets.cfg --status`. You can find more information in CrabSubmit/README. To merge and retrieve your job output, use `./crabControl.py --config datasets.cfg --merge --outputfolder /nfs/dust/cms/user/USER/Skim`

## FWLite analysis

To run the main FWLite analysis on all skimmed EDM files, run `./MainAnalysis/runMainAnalysis.py --edmfolder /nfs/dust/cms/user/USER/Skim`.
Without additional options, the FWLite analysis output will be placed in the `fwoutput` subfolder.

## Plotting

Modify and run `./Plotter/stackplotter.py`.
