#!/usr/bin/env python
import ROOT, os
from libstackplotter import stackplotter
import collections
import copy

# path to FWlite output files:
path = "/afs/desy.de/user/s/school73/school/heavyHiggs/CMSSW_8_0_10/src/H2hh2bbtautauAnalysis/MainAnalysis"

# select signal region / QCD region
dataDrivenQCD = True
suffix = ""
#~ dataDrivenQCD = False
#~ suffix = "_LS_Iso_selection"
logy = True

# specify file locations
rootFiles = collections.OrderedDict({})
rootFiles["Data"] = ['d', "%s/allData%s.root" % (path,suffix), ROOT.kBlack]
rootFiles["W jets"] = ['b', "%s/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8%s.root" % (path, suffix), ROOT.kBlue]
rootFiles["t#bar{t}"] = ['b', "%s/TT_TuneCUETP8M1_13TeV-powheg-pythia8%s.root" % (path, suffix), ROOT.kOrange]

for mass in [(250,ROOT.kRed),(400,ROOT.kViolet),(900,ROOT.kBlue)]:
    rootFiles["gg #rightarrow H*_{M=%i GeV} #rightarrow 2b2#tau" % mass[0]] = ['s', "%s/GluGluToRadionToHHTo2B2Tau_M-%i%s.root" % (path, mass[0], suffix), mass[1]]

#lumi = 12.906 #/fb, corresponds to ICHEP16 dataset (golden JSON)
lumi = 12.876 #/fb, Skim40

histoname = "weighted/nVert"
outputfilename = "nVert_%s" % suffix
stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "number of vertices",
             legendposition="0.49,0.67,0.95,0.9", ymin=1e-1, ymax=1e5, logy=logy)

histoname = "noPU/nVert"
outputfilename = "nVert_noPU_%s" % suffix
stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "number of vertices",
             legendposition="0.49,0.67,0.95,0.9", ymin=1e-1, ymax=1e5, logy=logy)

