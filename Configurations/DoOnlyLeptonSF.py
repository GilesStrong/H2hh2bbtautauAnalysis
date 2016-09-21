import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.PythonUtilities.LumiList as LumiList
import os
from Configuration.AlCa.autoCond import autoCond
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
from PhysicsTools.PatAlgos.tools.tauTools import *
#~ from RecoMET.METPUSubtraction.MVAMETConfiguration_cff import runMVAMET

options = VarParsing('analysis')

options.register ('runOnData', False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "run on data (True/False)")

options.register ('redoPuppi', False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "rebuild puppiMET (True/False)")

options.register ('globalTag', "80X_mcRun2_asymptotic_v5",
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'input global tag to be used')

options.register ('pileupGeneratedFile', "PUweightProducer/etc/pileupMC.root",
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "path to generated PU distribution")

options.register ('pileupDataFile', "PUweightProducer/etc/MyDataPileupHistogram_true.root",
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "path to data PU distribution")
                  
options.register ('muonIDScaleFactorFile', os.environ["CMSSW_BASE"]+"/src/HTT-utilities/LepEffInterface/data/Muon/Run2016BCD/Muon_IdIso0p15_eff.root",
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "muon ID scale factor input root file")

options.register ('muonTriggerScaleFactorFile', os.environ["CMSSW_BASE"]+"/src/HTT-utilities/LepEffInterface/data/Muon/Run2016BCD/Muon_IsoMu22_eff.root",
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "muon trigger scale factor input root file")

options.register ('svFitInputFile', "TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root",
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "SVFit input root file")

options.register ('CSVFile', "BTaggingSFProducer/etc/CSVv2_ichep.csv",
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "BTag scale factors")

options.register ('EffMapFile', "BTaggingSFProducer/plugins/EfficiencyMaps.root",
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "BTag efficiency map")

options.register ('EffMapFolder', "output_monitor",
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "BTag efficiency map folder name (same as crab shortname)")

options.parseArguments()

process = cms.Process("SKIM2")

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

#load for ReRunMet
process.load("RecoMET/METProducers.METSignificance_cfi")
process.load("RecoMET/METProducers.METSignificanceParams_cfi")

process.GlobalTag.globaltag = options.globalTag
jetCollection = "slimmedJets"

if options.runOnData:
    process.GlobalTag.globaltag = autoCond['run2_data']
    triggerResultsTag = "TriggerResults::HLT"
else:
    process.GlobalTag.globaltag = autoCond['run2_mc']
    triggerResultsTag = "TriggerResults::HLT2"

process.options  = cms.untracked.PSet()

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
            '/store/user/vkutzner/GluGluToRadionToHHTo2B2Tau_M-300_narrow_13TeV-madgraph/HeavyHiggsSkim14_GluGluToRadionToHHTo2B2Tau_M-300/160906_152802/0000/output_2.root'
        ))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(150) )

outputMonitorFile = options.outputFile.split(".root")[0] + "_monitor.root"

process.muonIDScaleFactor = cms.EDProducer('LeptonScaleFactorProducer',
                                           src       = cms.InputTag("selectedMuons"),
                                           inputFile = cms.string(options.muonIDScaleFactorFile),
                                           setEfficiencyDataToOne = cms.bool(options.runOnData),
                                           setEfficiencyMCToOne = cms.bool(options.runOnData),
                                           getErrors = cms.bool(False),
                                           debug = cms.bool(False),
                                           )

process.muonTriggerScaleFactor = cms.EDProducer('LeptonScaleFactorProducer',
                                                src       = cms.InputTag("selectedMuons"),
                                                inputFile = cms.string(options.muonTriggerScaleFactorFile),
                                                setEfficiencyDataToOne = cms.bool(options.runOnData),
                                                setEfficiencyMCToOne = cms.bool(True),
                                                getErrors = cms.bool(False),
                                                debug = cms.bool(False),
                                                )

process.out_objectselection = cms.OutputModule("PoolOutputModule",
                                               fileName = cms.untracked.string(options.outputFile),
        )

  
process.objectselection = cms.Path(
    process.muonIDScaleFactor
    *process.muonTriggerScaleFactor
    )


process.e = cms.EndPath(process.out_objectselection)
