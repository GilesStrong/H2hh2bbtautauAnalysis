#!/usr/bin/env python

# plotting script to plot data, signal and stacked MC backgrounds
#
# comments @ viktor.kutzner@cern.ch
#
# todo list: variable binning, ratio plots, automatic legend placement, TDR style

import os
from math import sqrt
import sys
import ROOT
import glob
import optparse
import collections

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetStatColor(ROOT.kWhite)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetPalette(1)
ROOT.TH1.SetDefaultSumw2(True)


def constructQCD(rootFiles, path, histoName, binSize, storage, lumi, noNegativeEntries=True):

    additionalFiles = ["_LS_Iso_selection"]

    hists = {}
    hists["default"] = collections.OrderedDict({})
    for additionalFile in additionalFiles:
        hists[additionalFile] = collections.OrderedDict({})

    # load additional input files needed for QCD determination:
    for label in rootFiles.keys():
        for additionalFile in additionalFiles:
            datasettype = rootFiles[label][0]
            filepath = rootFiles[label][1]
            rootFiles[label + additionalFile] = [datasettype, filepath.split(".root")[0] + additionalFile + ".root"]
 
    for label in rootFiles:

        fileName = rootFiles[label][1]
        print "@", fileName

        file=ROOT.TFile(fileName, "READ")
        storage.append(file)
        try:
            hist=file.Get(histoName)
            hist.SetName(hist.GetName()+"_"+label)
            currentBinWidth = int(hist.GetBinWidth(0))
            if binSize:
                hist.Rebin(binSize/currentBinWidth)
            else:
                binSize = currentBinWidth          

            default = True
            for additionalFile in additionalFiles:
                if additionalFile in fileName:
                    hists[additionalFile][label] = (hist)
                    default = False
            if default:
                hists["default"][label] = (hist)
            
        except:
            print "Couldn't load", fileName
            continue
            
    # added BG:
    for htype in hists:
        for label in hists[htype].keys():
            #~ if label == "addedData": continue
            if rootFiles[label][0] == 'b':
                if not "addedBG" in hists[htype][label]:
                    hists[htype]["addedBG"] = hists[htype][label].Clone()
                else:
                    hists[htype]["addedBG"].Add(hists[htype][label])

    # lumi weighting (luminosity assumed to be in fb, xsections in pb):
    #for addFile in hists:
        #hists[addFile]["addedBG"].Scale(lumi)
    
    # for like-sign leptons and isolation criteria (LS_Iso), get difference between data/added BG:
    ddQCD = hists["_LS_Iso_selection"]["Data_LS_Iso_selection"].Clone()
    bg = hists["_LS_Iso_selection"]["addedBG"].Clone()
    ddQCD.Add(bg, -1)
    
    # cut off negative entries:
    if noNegativeEntries:
        for iBin in range(ddQCD.GetNbinsX()):
            if ddQCD.GetBinContent(iBin) < 0:
                ddQCD.SetBinContent(iBin, 0)
    
    return ddQCD



def stackplotter(rootFiles, path, histoName, outputfilename, lumi, xtitle,
                 binSize=None, normalize=False, logx=False, logy=True, debug=False,
                 xmin=None, xmax=None, ymin=None, ymax=None, offsetx=1.2, offsety=1.2,
                 saveFormats="pdf", biglabel="CMS", smalllabel="private work", linewidth=2,
                 sqrts=13, legendposition="0.5,0.75,0.95,0.9", legendTextSize=0.6,
                 axisLabelSize=0.6, errorBarSize=4, dataDrivenQCD=False, qcdColor=ROOT.kRed):

    ROOT.gStyle.SetEndErrorSize(errorBarSize)
    
    storage = []
    hists = collections.OrderedDict({})
   
    for label in rootFiles:

        fileName = rootFiles[label][1]

        if "/" in fileName:
            file=ROOT.TFile(fileName, "READ")
        else:
            file=ROOT.TFile(path + "/" + fileName, "READ")
        storage.append(file)
        try:
            hist=file.Get(histoName)
            hist.SetName(hist.GetName()+"_"+label)
            currentBinWidth = int(hist.GetBinWidth(0))
            if binSize:
                hist.Rebin(binSize/currentBinWidth)
            else:
                binSize = currentBinWidth
            
            if normalize: hist.Scale( 1./hist.Integral() )
            
            hists[label] = (hist)
        except:
            continue

    # load data-driven QCD:
    if dataDrivenQCD:
        qcdstorage = []
        hists["QCD"] = constructQCD(rootFiles, path, histoName, binSize, qcdstorage, lumi)
        rootFiles["QCD"] = ['b', "", qcdColor]

    canvas = ROOT.TCanvas("canvas","canvas", 800, 800)

    pad1 = ROOT.TPad("pad1","pad1",0,0.2,1,1)
    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.2)
    pad1.Draw()
    pad2.Draw()
    
    ####################################################################
    pad1.cd()                               # main data/signal/MC plot #
    ####################################################################
    
    lpos = []
    for ipos in legendposition.split(","):
        lpos.append(float(ipos))
    legend = ROOT.TLegend(lpos[0],lpos[1],lpos[2],lpos[3])

    pad1.SetLogx(logx)
    pad1.SetLogy(logy)
    
    currentDrawOptions = ""

    # lumi weighting (luminosity assumed to be in fb, xsections in pb):
    #for label in hists:
        #if rootFiles[label][0] != 'd':
            #if label != 'QCD':
                #hists[label].Scale(lumi)
            
    # get data/MC ratio from added MC background histogram
    addedBG = None
    for label in hists:
        if rootFiles[label][0] == 'b':
            if not addedBG:
                addedBG = hists[label].Clone()
                continue
            addedBG.Add(hists[label])
    
    binmax = addedBG.GetMaximumBin()
    binmaxValue = addedBG.GetBinContent(binmax)
    
    if logy:
        binmaxValue = binmaxValue * 1e4
        ymin = 1e-8
    else:
        binmaxValue = binmaxValue + 50
        ymin = 0
    ymax = binmaxValue
    
    
    mcratio = None
    for label in hists:
        if rootFiles[label][0] == 'd':
            mcratio = hists[label].Clone()
            mcratio.Divide(addedBG)
    
    # plot backgrounds:
    mcstack = ROOT.THStack("mcstack","")
    for label in hists:
        if rootFiles[label][0] == 'b':
            if debug: print "adding background to stack:", label
            color = rootFiles[label][2]
            hists[label].SetFillColor(color)
            hists[label].SetLineColor(color)
            hists[label].SetLineWidth(0)
            hists[label].SetMarkerColor(color)
            mcstack.Add(hists[label])

    mcstack.SetTitle(";;Events/(%i GeV)" % binSize)
    mcstack.Draw("HIST " + currentDrawOptions)
    currentDrawOptions = "same"
    if debug:
        addedBG.SetLineColor(ROOT.kViolet)
        addedBG.SetLineWidth(linewidth)
        addedBG.Draw(currentDrawOptions)
    mcstack.GetYaxis().SetTitleOffset(offsety)
    mcstack.GetXaxis().SetTitleOffset(offsetx)
    if xmax != None and xmin != None:
        mcstack.GetXaxis().SetRangeUser(xmin,xmax)
    if ymin != None: mcstack.SetMinimum(ymin)
    if ymax != None: mcstack.SetMaximum(ymax)
    
    # plot data:
    for label in hists:
        if rootFiles[label][0] == 'd':
            hists[label].SetLineWidth(linewidth)
            hists[label].SetLineColor(rootFiles[label][2])
            hists[label].Sumw2()
            hists[label].Draw(currentDrawOptions + " E1")
            currentDrawOptions = "same"

    # plot signal:
    for label in hists:
        if rootFiles[label][0] == 's':
            hists[label].SetLineWidth(linewidth)
            hists[label].SetLineColor(rootFiles[label][2])
            hists[label].SetMarkerColor(rootFiles[label][2])
            hists[label].Draw(currentDrawOptions + " hist")
            currentDrawOptions = "same"

    for label in hists:
        legend.AddEntry(hists[label], label)
        if xmax != None and xmin != None:
            hists[label].GetXaxis().SetRangeUser(xmin,xmax)

    legend.UseCurrentStyle()
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.Draw()

    pad1.SetBottomMargin(0.012)
    pad1.SetLeftMargin(0.125)
    pad1.SetRightMargin(0.03)
    pad1.SetTopMargin(0.07)
    
    l = pad1.GetLeftMargin()
    t = pad1.GetTopMargin()
    r = pad1.GetRightMargin()
    b = pad1.GetBottomMargin()

    legend.SetTextSize(legendTextSize*t)
    mcstack.GetXaxis().SetTitleSize(0)
    mcstack.GetXaxis().SetLabelSize(0)
    mcstack.GetXaxis().SetLabelOffset(999)
    mcstack.GetYaxis().SetTitleSize(axisLabelSize*t)
    mcstack.GetYaxis().SetLabelSize(axisLabelSize*t)

    latex=ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)
    
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(axisLabelSize*t)

    latex.DrawLatex(1-r,1-t+0.2*t,"%.2f fb^{-1} (%i TeV)" % (lumi, sqrts))

    latex.SetTextAlign(11)
    latex.SetTextFont(61)
    latex.SetTextSize(0.75*t)
    latex.DrawLatex(0.2,0.83, biglabel)
    latex.SetTextSize(0.5*t)
    latex.SetTextFont(52)
    latex.DrawLatex(0.2,0.79, smalllabel)
    
    ####################################################################
    pad2.cd()                                  # Data/MC ratio subplot #
    ####################################################################
   
    pad2.SetGrid()
    pad2.SetLogx(logx)
    
    mcratio.Draw("E1")
    pad2.SetBottomMargin(0.45)
    pad2.SetLeftMargin(0.125)
    pad2.SetRightMargin(0.03)
    pad2.SetTopMargin(0)
    
    if xmax != None and xmin != None:
        mcratio.GetXaxis().SetRangeUser(xmin,xmax)
    
    t=t*4
    mcratio.GetXaxis().SetTitleSize(axisLabelSize*t)
    mcratio.GetXaxis().SetLabelSize(axisLabelSize*t)
    mcratio.GetYaxis().SetTitleSize(axisLabelSize*t)
    mcratio.GetYaxis().SetLabelSize(axisLabelSize*t)
    mcratio.GetYaxis().SetRangeUser(0.01,1.99)
    mcratio.GetYaxis().SetNdivisions(4)
    mcratio.GetYaxis().SetTitleOffset(0.32)
    
    mcratio.SetLineWidth(linewidth)
    mcratio.SetLineColor(ROOT.kBlack)
    
    mcratio.SetTitle(";" + xtitle + ";Data/MC")
        
    canvas.cd()
    
    canvas.Update()
    for saveFormat in saveFormats.split(","):
        canvas.SaveAs(outputfilename + "." + saveFormat)
        
    # print dd qcd:
    if "QCD" in hists:
        # hists["QCD"].Draw()
        # canvas.SaveAs(outputfilename + "_qcd.pdf")
        folder=histoName.split("/")[0]
        histo=histoName.split("/")[-1]
        outputfile = ROOT.TFile(outputfilename + "_qcd.root","UPDATE")
        outputfolder = outputfile.mkdir(folder)
        outputfolder.cd()
        hists["QCD"].SetName(histo)
        hists["QCD"].Write()
        outputfile.Close() 
                
    return 0
    
        
