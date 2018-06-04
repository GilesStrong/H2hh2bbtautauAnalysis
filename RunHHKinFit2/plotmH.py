#! /usr/bin/env python

import ROOT
import sys
from DataFormats.FWLite import Events, Handle

# Make VarParsing object
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideAboutPythonConfigFile#VarParsing_Example
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.parseArguments()

# Events takes either
# - single file name
# - list of file names
# - VarParsing options

# use Varparsing object
events = Events (options)

# create handle outside of loop
handle_mH  = Handle ("double")


# Create histograms, etc.
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gROOT.SetStyle('Plain') # white background
h_mH = ROOT.TH1F ("mH", "KinFit mass;m_H;entries", 100,200,400)

# loop over events
for event in events:
    event.getByLabel ("kinfit", "mH", handle_mH)
    mH = handle_mH.product()[0]

    h_mH.Fill(mH)

# make a canvas, draw, and save it
c1 = ROOT.TCanvas()
h_mH.Draw()
c1.Print ("mH.png")
