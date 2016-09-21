import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("SCALE")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring16MiniAODv2/GluGluToRadionToHHTo2B2Tau_M-300_narrow_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/80000/02EE1330-5D3B-E611-B5AD-001CC4A7C0A4.root'
#'file:/pnfs/desy.de/cms/tier2/store/data/Run2016D/SingleMuon/MINIAOD/PromptReco-v2/000/276/317/00000/04D90ACF-0C45-E611-92CE-02163E0136AD.root'
    )
)

process.muonsWithIso = cms.EDProducer("PATMuonIsolationProducer",
    muons = cms.InputTag('slimmedMuons'),
    isoAlgorithm = cms.int32(7)
)

process.selectedMuons = cms.EDFilter("PATMuonCloneSelectorFilter",
    src = cms.InputTag('muonsWithIso'),
    vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
    ptmin = cms.double(20),
    etamax = cms.double(2.1),   
    vtxdxymax = cms.double(0.045),                 
    vtxdzmax = cms.double(0.2),
    takeMuonID = cms.string("medium"),
    iso = cms.double(0.15),
    isoAlgorithm = cms.int32(7),
    useIso = cms.bool(True),
    invertIso = cms.bool(False),
    debug = cms.bool(False),
    filter = cms.bool(True)
)

process.muonScaleFactor = cms.EDProducer('LeptonScaleFactorProducer',
     src       = cms.InputTag( "selectedMuons" ),
     inputFile = cms.string( os.environ["CMSSW_BASE"]+"/src/HTT-utilities/LepEffInterface/data/Muon/Run2016BCD/Muon_IsoMu22_eff.root" ),
#     inputFile = cms.string( os.environ["CMSSW_BASE"]+"/src/HTT-utilities/LepEffInterface/data/Muon/Run2016BCD/Muon_IdIso0p15_eff.root" ),
     setEfficiencyDataToOne = cms.bool( False ),
     setEfficiencyMCToOne = cms.bool( True ),
     getErrors = cms.bool( False ),
     debug = cms.bool( True ),
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_slimmedMuons_*_*',
        'keep *_*_*_SCALE'
    ),
    SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring("p")
    )
)

process.p = cms.Path(process.muonsWithIso*process.selectedMuons*process.muonScaleFactor)
process.e = cms.EndPath(process.out)
