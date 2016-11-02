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

handle_muons  = Handle ("std::vector<pat::Muon>")
handle_taus   = Handle ("std::vector<pat::Tau>")
handle_jets   = Handle ("std::vector<pat::Jet>")
handle_mets   = Handle ("std::vector<pat::MET>")


# Create histograms, etc.
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gROOT.SetStyle('Plain') # white background
h_mH_kinfit = ROOT.TH1F ("mHkinfit", "reconstructed Higgs mass;m_{H,reco};entries", 100,200,400)
h_mH_4body  = ROOT.TH1F ("mH4body",  "reconstructed Higgs mass;m_{H,reco};entries", 100,200,400)

# loop over events
for event in events:
    event.getByLabel("kinfit", "mH", handle_mH)
    mH = handle_mH.product()[0]
    h_mH_kinfit.Fill(mH)

    event.getByLabel("selectedMuons", handle_muons)
    event.getByLabel("selectedTaus",  handle_taus)
    event.getByLabel("selectedJets",  handle_jets)
    event.getByLabel("slimmedMETs",   handle_mets)
    mu=handle_muons.product()[0]
    tau=handle_taus.product()[0]
    jet1=handle_jets.product()[0]
    jet2=handle_jets.product()[1]
    met=handle_mets.product()[0]

    mu4v = ROOT.TLorentzVector(mu.px(),mu.py(),mu.pz(),mu.energy())
    tau4v = ROOT.TLorentzVector(tau.px(),tau.py(),tau.pz(),tau.energy())
    jet14v = ROOT.TLorentzVector(jet1.px(),jet1.py(),jet1.pz(),jet1.energy())
    jet24v = ROOT.TLorentzVector(jet2.px(),jet2.py(),jet2.pz(),jet2.energy())
    met4v = ROOT.TLorentzVector(met.px(),met.py(),met.pz(),met.energy())

    h_mH_4body.Fill((mu4v+tau4v+jet14v+jet24v+met4v).M())


# make a canvas, draw, and save it
c1 = ROOT.TCanvas()
h_mH_kinfit.SetLineWidth(2)
h_mH_kinfit.SetLineColor(ROOT.kBlue)
h_mH_kinfit.Draw()

h_mH_4body.SetLineWidth(2)
h_mH_4body.SetLineColor(ROOT.kRed)
h_mH_4body.Draw("same")
c1.Print ("mH.png")