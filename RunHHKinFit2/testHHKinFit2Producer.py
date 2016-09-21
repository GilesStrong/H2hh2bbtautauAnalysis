import FWCore.ParameterSet.Config as cms
import os



process = cms.Process("HHKINFIT2")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(2)


process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring16MiniAODv2/GluGluToRadionToHHTo2B2Tau_M-300_narrow_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/80000/02EE1330-5D3B-E611-B5AD-001CC4A7C0A4.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring16MiniAODv2/GluGluToRadionToHHTo2B2Tau_M-300_narrow_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/80000/2048A454-5D3B-E611-8D46-24BE05C3DB21.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring16MiniAODv2/GluGluToRadionToHHTo2B2Tau_M-300_narrow_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/80000/2E6C8F37-5D3B-E611-9C0A-24BE05C6C741.root'
#'file:/pnfs/desy.de/cms/tier2/store/data/Run2016D/SingleMuon/MINIAOD/PromptReco-v2/000/276/317/00000/04D90ACF-0C45-E611-92CE-02163E0136AD.root'
    )
)
process.muonsWithIso = cms.EDProducer("PATMuonIsolationProducer",
    muons = cms.InputTag('slimmedMuons'),
    isoAlgorithm = cms.int32(7)
)

process.selectedMuons = cms.EDFilter("PATMuonCloneSelectorFilter",
    src = cms.InputTag('muonsWithIso'),
    #~ src = cms.InputTag('muonsWithIso'),
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

process.selectedTaus = cms.EDFilter("PATTauCloneSelectorFilter",
    src = cms.InputTag('slimmedTaus'),
    ptmin = cms.double(20),
    etamax = cms.double(2.3),
    decayModeFinding = cms.double(0.5),
    leadTauCandDZ = cms.double(0.2),
    iso = cms.double(0.5),
    isoCrit = cms.string("byTightIsolationMVArun2v1DBoldDMwLT"),
    electronRejection = cms.string("againstElectronVLooseMVA6"),
    muonRejection = cms.string("againstMuonTight3"),
    useIso = cms.bool(True),
    useElectronRejection = cms.bool(True),
    useMuonRejection = cms.bool(True),
    debug = cms.bool(False),
    filter = cms.bool(True)
)

process.selectedJets = cms.EDFilter("PATJetCloneSelectorFilter",
    src = cms.InputTag('slimmedJets'),  
    ptmin = cms.double(20),
    etamax = cms.double(4.7),   

    neutralhadronfractionmax = cms.double(0.99),
    neutralEMfractionmax = cms.double(0.99),
    constituentsmin = cms.int32(1),
    muonfrationmax = cms.double(9999),
    chargedhadronfractionmin = cms.double(0.00),
    chargedmultiplicitymin = cms.int32(0),
    chargedEMfractionmax = cms.double(0.99),
    neutralEMfractionforwardmax = cms.double(0.90),
    neutralparticlesforwardmin = cms.int32(10),

    checkoverlap_muons = cms.bool(True),
    muons = cms.InputTag('selectedMuons'),
    dR_muon = cms.double(0.5),

    checkoverlap_taus = cms.bool(True),
    taus = cms.InputTag('selectedTaus'),
    dR_tau = cms.double(0.5),

    usebtag = cms.bool(True),
    btagptmin = cms.double(20),
    btagetamax = cms.double(2.4),
    btagdesc = cms.double(0.8),
    btagalgo =    cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    invertbtag = cms.bool(False),

    debug = cms.bool(False),
    filter = cms.bool(True)
)

process.selectedElectrons = cms.EDFilter("PATElectronCloneSelectorFilter",
     src = cms.InputTag('slimmedElectrons'),
     vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
     electronID = cms.string("mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
     ptmin = cms.double(24),
     etamax = cms.double(2.1),   
     vtxdxymax = cms.double(0.045),                 
     vtxdzmax = cms.double(0.2),
     iso = cms.double(0.1),
     useIso = cms.bool(True),
     invertIso = cms.bool(False),
     noMatchedConversions = cms.bool(True),
     maxMissingHits = cms.int32(1),
     debug = cms.bool(False),
     filter = cms.bool(False)
)

process.jetfilter = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedJets'),
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(1000),
)

process.jetfilter2 = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedJets'),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(1),
)

process.taufilter = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedTaus'),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1),
)

process.muonfilter = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedMuons'),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1),
)

process.electronveto = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedElectrons'),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(0),
)

process.oppositechargefilter = cms.EDFilter("PATMuonTauOppositeChargeFilter",
    muons = cms.InputTag('selectedMuons'),
    taus = cms.InputTag('selectedTaus'),
    debug = cms.bool(False),
    filter = cms.bool(True)
)

process.kinfit = cms.EDProducer('HHKinFit2Producer',
    muons = cms.InputTag('selectedMuons'),
    taus  = cms.InputTag('selectedTaus'),
    bjets = cms.InputTag('selectedJets'),
    met   = cms.InputTag('slimmedMETs'),
    metCovMat = cms.InputTag('METSignificance','METCovariance'),
    debug = cms.bool(False),
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("test.root"),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep GenRunInfoProduct_generator__SIM',
#        'keep *_slimmedMuons_*_*',
#        'keep *_slimmedJets_*_*',
#        'keep *_slimmedTaus_*_*',
        'keep *_selectedMuons_*_*',
        'keep *_selectedTaus_*_*',
        'keep *_selectedJets_*_*',
        'keep *_prunedGenParticles_*_*',
        'keep *_slimmedGenJets_*_*',
        'keep *_nEvents*_*_*',
        'keep *_*_*_HHKINFIT2',
        'keep *_TriggerResults_*_SKIM'
    ),
    SelectEvents = cms.untracked.PSet( 
       SelectEvents = cms.vstring("objectselection")
    )
)

process.objectselection = cms.Path(
                     process.muonsWithIso
                     *process.selectedMuons
                     *process.selectedTaus
                     *process.selectedJets
                     *process.muonfilter
                     *process.taufilter
                     *process.jetfilter
                     *process.oppositechargefilter
                     *process.selectedElectrons
                     *process.electronveto
                     *process.kinfit
                     )

process.e = cms.EndPath(process.out)
