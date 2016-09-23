#!/usr/bin/env python
import ROOT, os
from libstackplotter import stackplotter
import collections
import copy

# path to FWlite output files:
path = "/afs/desy.de/user/s/school73/school/heavyHiggs/CMSSW_8_0_10/src/H2hh2bbtautauAnalysis/MainAnalysis"

# select signal region / QCD region
for region in range(4):
    dataDrivenQCD = None
    suffix = None
#A
    if region == 3:
        dataDrivenQCD = False
        suffix = ""
#B
    if region == 2:
        dataDrivenQCD = False
        suffix = "_LS_Iso_selection"
#D
    if region == 1:
        dataDrivenQCD = False
        suffix = "_LS_InvIso_selection"
#C
    if region == 0:
        dataDrivenQCD = False
        suffix = "_OS_InvIso_selection"
        
    logy = False
    print region
# specify file locations
    rootFiles = collections.OrderedDict({})
    rootFiles["Data"] = ['d', "%s/allData%s.root" % (path,suffix), ROOT.kBlack]
    rootFiles["W jets"] = ['b', "%s/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8%s.root" % (path, suffix), ROOT.kBlue]
    rootFiles["t#bar{t}"] = ['b', "%s/TT_TuneCUETP8M1_13TeV-powheg-pythia8%s.root" % (path, suffix), ROOT.kOrange]
    for mass in [(250,ROOT.kRed),(400,ROOT.kViolet),(900,ROOT.kBlue)]:
       rootFiles["gg #rightarrow H*_{M=%i GeV} #rightarrow 2b2#tau" % mass[0]] = ['s', "%s/GluGluToRadionToHHTo2B2Tau_M-%i%s.root" % (path, mass[0], suffix), mass[1]]

    lumi = 12.906 #/fb, corresponds to ICHEP16 dataset (golden JSON)
#~ lumi = 12.876 #/fb, Skim40
#lumi = 1

#~ histoname = "weighted/h_monitor"
#~ outputfilename = "Plotter/monitor_%s" % suffix
#~ stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, 1e-3, "",
             #~ legendposition="0.49,0.67,0.95,0.9", ymin=1e-2, ymax=5e5,
             #~ dataDrivenQCD=False, normalize=True, logy=True)

    histoname = "weighted/muonPt"
    outputfilename = "muonpt_%s" % suffix
    stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "p_{T}^{#mu} (GeV)",
                 legendposition="0.49,0.67,0.95,0.9", binSize=10, ymin=1e-2, ymax=5e5,
                 dataDrivenQCD=dataDrivenQCD, logy=logy)

    histoname = "weighted/m_met"
    outputfilename = "met_%s" % suffix
    stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "E_{T}^{miss} (GeV)",
                 legendposition="0.49,0.67,0.95,0.9", binSize=10, ymin=1e-2, ymax=5e5,
                 dataDrivenQCD=dataDrivenQCD, logy=logy)

    histoname = "weighted/svFitMass"
    outputfilename = "svFitMass_%s" % suffix
    stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "m_{SV} (GeV)",
                 legendposition="0.49,0.67,0.95,0.9", binSize=10, ymin=1e-5, ymax=5e5,
                 dataDrivenQCD=dataDrivenQCD, logy=logy)

    histoname = "weighted/m_kinfit"
    outputfilename = "m_kinfit_%s" % suffix
    stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "m_{kinfit} (GeV)",
                 legendposition="0.49,0.67,0.95,0.9", binSize=10, ymin=1e-2, ymax=5e5,
                 dataDrivenQCD=dataDrivenQCD, logy=logy)

    histoname = "weighted/P_chi2"
    outputfilename = "P_chi2_%s" % suffix
    stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "p",
                 legendposition="0.49,0.67,0.95,0.9", ymin=1e-2, ymax=5e5,
                 dataDrivenQCD=dataDrivenQCD, logy=logy)

    histoname = "weighted/m_bbar"
    outputfilename = "m_bbar_%s" % suffix
    stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "m_{b#bar{b}} (GeV)",
                 legendposition="0.49,0.67,0.95,0.9", binSize=10, ymin=1e-2, ymax=5e5,
                 dataDrivenQCD=dataDrivenQCD, logy=logy)


    histoname = "weighted/m_tauvis"
    outputfilename = "m_tauvis_%s" % suffix
    stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "m_{#tau vis} (GeV)",
                 legendposition="0.49,0.67,0.95,0.9", binSize=10, ymin=1e-2, ymax=5e5,
                 dataDrivenQCD=dataDrivenQCD, logy=logy)

    histoname = "weighted/m_t_mu"
    outputfilename = "m_t_mu_%s" % suffix
    stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "m (GeV)",
                 legendposition="0.49,0.67,0.95,0.9", binSize=10, ymin=1e-2, ymax=5e5,
                 dataDrivenQCD=dataDrivenQCD, logy=logy)

    histoname = "weighted/nVert"
    outputfilename = "nVert_%s" % suffix
    stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "number of vertices",
                 legendposition="0.49,0.67,0.95,0.9", ymin=1e-1, ymax=1e5, logy=logy)

    histoname = "noPU/nVert"
    outputfilename = "nVert_noPU_%s" % suffix
    stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "number of vertices",
                 legendposition="0.49,0.67,0.95,0.9", ymin=1e-1, ymax=1e5, logy=logy)

# plot the covariance matrix
#~ histoname = "weighted/h_cov_xx"
#~ outputfilename = "Plotter/h_cov_xx_%s" % suffix
#~ stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "covm",
             #~ legendposition="0.49,0.67,0.95,0.9", ymin=1e-1, ymax=1e5, logy=logy)
#~ histoname = "weighted/h_cov_xy"
#~ outputfilename = "Plotter/h_cov_xy_%s" % suffix
#~ stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "covm",
             #~ legendposition="0.49,0.67,0.95,0.9", ymin=1e-1, ymax=1e5, logy=logy)
#~ histoname = "weighted/h_cov_yx"
#~ outputfilename = "Plotter/h_cov_yx_%s" % suffix
#~ stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "covm",
             #~ legendposition="0.49,0.67,0.95,0.9", ymin=1e-1, ymax=1e5, logy=logy)
#~ histoname = "weighted/h_cov_yy"
#~ outputfilename = "Plotter/h_cov_yy_%s" % suffix
#~ stackplotter(rootFiles.copy(), os.getcwd(), histoname, outputfilename, lumi, "covm",
             #~ legendposition="0.49,0.67,0.95,0.9", ymin=1e-1, ymax=1e5, logy=logy)





