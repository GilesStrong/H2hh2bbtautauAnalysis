#!/usr/bin/env python

# adapted script from Dinko Ferencek, modified by Viktor Kutzner

# produces p_T/eta efficiency map from BTaggingEffAnalyzer output

from __future__ import division
import os, sys
import ROOT
from glob import glob
from array import array
import optparse

ROOT.gROOT.SetBatch(1)
outputFile = 0

def produceEfficiencyMaps(filename, subdirectory, binning):

  try:

      inputFile = ROOT.TFile(filename, 'READ')
      dataset = filename.split(".root")[0].split("/")[-1].split("_monitor")[0]

      print "Doing %s" % dataset

      outputFile.cd()
      outputFile.mkdir(dataset)
      outputFile.cd(dataset)

      for partonFlavor in ['b', 'c', 'udsg']:

        denominatorHisto = subdirectory + '/h2_BTaggingEff_Denom_' + partonFlavor
        numeratorHisto = subdirectory + '/h2_BTaggingEff_Num_' + partonFlavor

        denominatorIn = inputFile.Get(denominatorHisto)
        numeratorIn = inputFile.Get(numeratorHisto)

        xShift = denominatorIn.GetXaxis().GetBinWidth(1)/2.
        yShift = denominatorIn.GetYaxis().GetBinWidth(1)/2.

        binsX = array('d', binning[partonFlavor][0])
        binsY = array('d', binning[partonFlavor][1])

        denominatorOut = ROOT.TH2D('denominator_' + partonFlavor, '', (len(binsX)-1), binsX, (len(binsY)-1), binsY)
        numeratorOut   = ROOT.TH2D('numerator_' + partonFlavor, '', (len(binsX)-1), binsX, (len(binsY)-1), binsY)
        efficiencyOut  = ROOT.TH2D('efficiency_' + partonFlavor, '', (len(binsX)-1), binsX, (len(binsY)-1), binsY)

        # loop over all bins
        for i in range(1,denominatorOut.GetXaxis().GetNbins()+1):
          for j in range(1,denominatorOut.GetYaxis().GetNbins()+1):

            binXMin = denominatorIn.GetXaxis().FindBin(denominatorOut.GetXaxis().GetBinLowEdge(i)+xShift)
            binXMax = denominatorIn.GetXaxis().FindBin(denominatorOut.GetXaxis().GetBinUpEdge(i)-xShift)
            binYMinPos = denominatorIn.GetYaxis().FindBin(denominatorOut.GetYaxis().GetBinLowEdge(j)+yShift)
            binYMaxPos = denominatorIn.GetYaxis().FindBin(denominatorOut.GetYaxis().GetBinUpEdge(j)-yShift)
            binYMinNeg = denominatorIn.GetYaxis().FindBin(-denominatorOut.GetYaxis().GetBinUpEdge(j)+yShift)
            binYMaxNeg = denominatorIn.GetYaxis().FindBin(-denominatorOut.GetYaxis().GetBinLowEdge(j)-yShift)

            denominator = denominatorIn.Integral(binXMin,binXMax,binYMinPos,binYMaxPos)
            denominator = denominator + denominatorIn.Integral(binXMin,binXMax,binYMinNeg,binYMaxNeg)
            numerator = numeratorIn.Integral(binXMin,binXMax,binYMinPos,binYMaxPos)
            numerator = numerator + numeratorIn.Integral(binXMin,binXMax,binYMinNeg,binYMaxNeg)

            denominatorOut.SetBinContent(i,j,denominator)
            numeratorOut.SetBinContent(i,j,numerator)
            if(denominator>0.): efficiencyOut.SetBinContent(i,j,numerator/denominator)

        # check if there are any bins with 0 or 100% efficiency
        for i in range(1,denominatorOut.GetXaxis().GetNbins()+1):
          for j in range(1,denominatorOut.GetYaxis().GetNbins()+1):
             efficiency = efficiencyOut.GetBinContent(i,j)
             if(efficiency==0. or efficiency==1.):
                print 'Warning! Bin(%i,%i) for %s jets has a b-tagging efficiency of %.3f'%(i,j,partonFlavor,efficiency)

        # set efficiencies in overflow bins
        for i in range(1,denominatorOut.GetXaxis().GetNbins()+1):
          efficiencyOut.SetBinContent(i, denominatorOut.GetYaxis().GetNbins()+1, efficiencyOut.GetBinContent(i, denominatorOut.GetYaxis().GetNbins()))

        for j in range(1,denominatorOut.GetYaxis().GetNbins()+2):
          efficiencyOut.SetBinContent(denominatorOut.GetXaxis().GetNbins()+1, j, efficiencyOut.GetBinContent(denominatorOut.GetXaxis().GetNbins(), j))

        denominatorOut.Write()
        numeratorOut.Write()
        efficiencyOut.Write()

      print "b-tagging efficiency map for %s successfully created." % (filename)

  except:
      print "Skipping %s" % filename


if __name__ == "__main__":

    parser = optparse.OptionParser()

    parser.add_option("--edmFiles", action="store", type="string", default="/nfs/dust/cms/user/kutznerv/Skim35/*monitor.root", dest="edmFiles")
    (options, pfade) = parser.parse_args()

    binning = {'b':    [[0., 40., 60., 80., 100., 150., 200., 300., 1000.],[0., 0.6, 1.2, 2.4]],
             'c':    [[0., 40., 60., 80., 100., 150., 200., 1000.],[0., 0.6, 1.2, 2.4]],
             'udsg': [[0., 40., 60., 80., 100., 150., 200., 1000.],[0., 0.6, 1.2, 2.4]]}

    outputFile = ROOT.TFile("EfficiencyMaps.root", 'RECREATE')
    for filename in glob(options.edmFiles):
        produceEfficiencyMaps(filename, "bTaggingEffAnalyzerAK8PF", binning)
    outputFile.Close()
    
    outputFile = ROOT.TFile("EfficiencyMaps_LS_Iso.root", 'RECREATE')
    for filename in glob(options.edmFiles):
        produceEfficiencyMaps(filename, "bTaggingEffAnalyzerAK8PF_LS_Iso", binning)
    outputFile.Close()
    
    outputFile = ROOT.TFile("EfficiencyMaps_LS_InvIso.root", 'RECREATE')
    for filename in glob(options.edmFiles):
        produceEfficiencyMaps(filename, "bTaggingEffAnalyzerAK8PF_LS_InvIso", binning)
    outputFile.Close()

    
    outputFile = ROOT.TFile("EfficiencyMaps_OS_InvIso.root", 'RECREATE')
    for filename in glob(options.edmFiles):
        produceEfficiencyMaps(filename, "bTaggingEffAnalyzerAK8PF_OS_InvIso", binning)
    outputFile.Close()

