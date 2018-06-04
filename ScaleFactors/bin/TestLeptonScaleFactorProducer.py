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
handle_muons  = Handle ("std::vector<pat::Muon>")
handle_effval  = Handle ("double")


# Create histograms, etc.
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gROOT.SetStyle('Plain') # white background
effHist = ROOT.TH1F ("muoneff", "Muon Efficiency", 100, 0, 1)

# loop over events
for event in events:
    event.getByLabel ("selectedMuons", handle_muons)
    muons = handle_muons.product()

    event.getByLabel ("muonScaleFactor", "EfficiencyValue", handle_effval)
    muoneff = handle_effval.product()[0]

    # use muons to make Z peak
    for muon in muons:
        print(muon.pt())

    print(muoneff)
    print("--------")

    effHist.Fill(muoneff)

# make a canvas, draw, and save it
c1 = ROOT.TCanvas()
effHist.Draw()
c1.Print ("eff_py.png")
