#!/usr/bin/env python
import ROOT, os
from libstackplotter import stackplotter
import collections
import copy

# path to FWlite output files:
path = "/nfs/dust/cms/user/kutznerv/Skim40/fwoutput"
logy = True

# specify file locations
rootFiles = collections.OrderedDict({})
rootFiles["Data"] = ['d', "%s/allData.root" % (path), ROOT.kBlack]
rootFiles["W jets"] = ['b', "%s/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root" % (path), ROOT.kBlue]

lumi = 12.876 #/fb, Skim40

histoname = "noPU/nVert"
outputfilename = "Plotter/nVert_noPU"
stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "number of vertices",
             legendposition="0.49,0.67,0.95,0.9", ymin=1e-1, ymax=1e5, logy=logy)


    
    
    
