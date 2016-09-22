#!/usr/bin/env python                                                                                                                                                                                          
import ROOT, os
from libstackplotter import stackplotter
import collections

# path to FWlite output files:                                                                                                                                                                                 
path = "/afs/desy.de/user/s/school73/school/heavyHiggs/CMSSW_8_0_10/src/H2hh2bbtautauAnalysis/MainAnalysis/"

rootFiles = collections.OrderedDict({})
rootFiles["Data"] = ['d', "%s/allData.root" % (path), ROOT.kBlack]
rootFiles["W jets"] = ['b', "%s/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root" % (path), ROOT.kBlue]
rootFiles["t#bar{t}"] = ['b', "%s/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root" % (path), ROOT.kOrange]

lumi = 12.906 #/fb  ICHEP16 dataset                                                                                                                                                                            

histoname = "weighted/nVert"
outputfilename = "nVert"
stackplotter(rootFiles, os.getcwd(), histoname, outputfilename, lumi, "number of vertices",
             legendposition="0.49,0.67,0.95,0.9", binSize=1, ymin=1e-1, ymax=1e5)
