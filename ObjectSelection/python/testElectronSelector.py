import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

#process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
#        'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring16MiniAODv2/ZH_ZToMM_HToInvisible_M150_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/78ACFC0C-C826-E611-BCE9-0025901E4A16.root'
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring16MiniAODv2/GluGluToRadionToHHTo2B2Tau_M-300_narrow_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/80000/02EE1330-5D3B-E611-B5AD-001CC4A7C0A4.root'
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeSelection.root')
)

process.muonsWithIso = cms.EDProducer("PATMuonIsolationProducer",
    muons = cms.InputTag('slimmedMuons'),
    isoAlgorithm = cms.int32(7)
)

process.unpackedSelectedPatTrigger = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
              patTriggerObjectsStandAlone = cms.InputTag( 'selectedPatTrigger' ),
              triggerResults              = cms.InputTag( 'TriggerResults::HLT2' )
)

process.matchedPatMuonTriggerMatchHLTIsoMu22 = cms.EDProducer("PATTriggerMatcherDRLessByR" , 
              src     = cms.InputTag( "muonsWithIso" ),
              matched = cms.InputTag( "unpackedSelectedPatTrigger" ),          
              matchedCuts = cms.string( 'type( "TriggerMuon" ) && path( "HLT_IsoMu22_v*" )' ),
              #maxDPtRel = cms.double( 0.5 ),                  # in case you want to use 'PATTriggerMatcherDRDPtLessByR'
              maxDeltaR = cms.double( 0.5 ),
              resolveAmbiguities    = cms.bool( True ),        # only one match per trigger object
              resolveByMatchQuality = cms.bool( True )         # take best match found per reco object
)

process.allMuons = cms.EDProducer("PATTriggerMatchMuonEmbedder",
              src = cms.InputTag( 'muonsWithIso' ),
              matches = cms.VInputTag(
              cms.InputTag( 'matchedPatMuonTriggerMatchHLTIsoMu22' )
              )
)

process.triggeredMuons = cms.EDFilter("PATMuonTriggerMatchCloneSelectorFilter",
              src = cms.InputTag( 'allMuons' ),
              triggerpath = cms.string('HLT_IsoMu22_v*'),
              filter = cms.bool(True)
              )

process.selectedMuons = cms.EDFilter("PATMuonCloneSelectorFilter",
    src = cms.InputTag('triggeredMuons'),
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

process.selected = cms.EDFilter("PATElectronCloneSelectorFilter",
    src = cms.InputTag('slimmedElectrons'),
    vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
    electronID = cms.string("mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
    # (miniaod2 contains HEEPV60, cut-based, MVA electronIDs)
    ptmin = cms.double(24),
    etamax = cms.double(2.1),   
    vtxdxymax = cms.double(0.045),                 
    vtxdzmax = cms.double(0.2),
    iso = cms.double(0.1),
    useIso = cms.bool(True),
    invertIso = cms.bool(False),
    noMatchedConversions = cms.bool(True),
    maxMissingHits = cms.int32(1),
    debug = cms.bool(True),
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

process.oppositechargefilter = cms.EDFilter("PATMuonTauOppositeChargeFilter",
    muons = cms.InputTag('selectedMuons'),
    taus = cms.InputTag('selectedTaus'),
    debug = cms.bool(False),
    filter = cms.bool(True)
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep GenRunInfoProduct_generator__SIM',
        'keep *_generator_*_*',
#        'keep *_slimmedMuons_*_*',
#        'keep *_slimmedJets_*_*',
#        'keep *_slimmedTaus_*_*',
        'keep *_selectedElectrons_*_*',
        'keep *_selectedMuons_*_*',
        'keep *_selectedTaus_*_*',
        'keep *_selectedJets_*_*',
        'keep *_prunedGenParticles_*_*',
        'keep *_slimmedGenJets_*_*'
    ),
    SelectEvents = cms.untracked.PSet( 
       SelectEvents = cms.vstring("objectselection","p","pEle")
    )
)

process.mon1 = cms.EDAnalyzer("SelectionAnalyzer",
              muons = cms.untracked.InputTag( 'slimmedMuons' ),
              taus  = cms.untracked.InputTag( 'slimmedTaus' ),
              jets  = cms.untracked.InputTag( 'slimmedJets' ),
)

process.mon2 = process.mon1.clone( muons = cms.untracked.InputTag( 'triggeredMuons' ))
process.mon3 = process.mon2.clone( muons = cms.untracked.InputTag( 'selectedMuons' ))
process.mon4 = process.mon3.clone( taus = cms.untracked.InputTag( 'selectedTaus' ))
process.mon5 = process.mon4.clone( jets = cms.untracked.InputTag( 'selectedJets' ))
process.mon6 = process.mon5.clone( )


process.objectselection = cms.Path(
                     process.mon1
                     *process.muonsWithIso
                     *process.unpackedSelectedPatTrigger
                     *process.matchedPatMuonTriggerMatchHLTIsoMu22
                     *process.allMuons
                     *process.triggeredMuons
                     *process.mon2
                     *process.selectedMuons
                     *process.mon3
                     *process.selectedTaus
                     *process.mon4
                     *process.selectedJets
                     *process.mon5
                     *process.muonfilter
                     *process.taufilter
                     *process.jetfilter
                     *process.oppositechargefilter
                     *process.mon6
                     )

process.p = cms.Path(
                     process.muonsWithIso
                     *process.unpackedSelectedPatTrigger
                     *process.matchedPatMuonTriggerMatchHLTIsoMu22
                     *process.allMuons
                     *process.triggeredMuons
                     *process.selectedMuons
                     *process.selectedTaus
                     *process.selectedJets
                     *process.muonfilter
                     *process.taufilter
                     *process.jetfilter2
                     *process.oppositechargefilter
                     )

process.pEle = cms.Path(process.selectedElectrons)

process.e = cms.EndPath(process.out)
