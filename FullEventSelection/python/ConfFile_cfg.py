import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.autoCond import autoCond
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
from PhysicsTools.PatAlgos.tools.tauTools import *

process = cms.Process("SKIM")

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

process.GlobalTag.globaltag = autoCond['run2_mc']

process.options  = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                       allowUnscheduled = cms.untracked.bool(True),
                                       numberOfThreads=cms.untracked.uint32(10)
                                       )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring16MiniAODv2/GluGluToRadionToHHTo2B2Tau_M-300_narrow_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/80000/02EE1330-5D3B-E611-B5AD-001CC4A7C0A4.root',
        'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring16MiniAODv2/GluGluToRadionToHHTo2B2Tau_M-300_narrow_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/80000/2048A454-5D3B-E611-8D46-24BE05C3DB21.root',
        'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring16MiniAODv2/GluGluToRadionToHHTo2B2Tau_M-300_narrow_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/80000/2E6C8F37-5D3B-E611-9C0A-24BE05C6C741.root',
        'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring16MiniAODv2/GluGluToRadionToHHTo2B2Tau_M-300_narrow_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/80000/6ABB144F-5D3B-E611-BF90-A0000420FE80.root'
        ))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("fullSelection_monitoring.root")
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
                                    invertIso = cms.bool(False),
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
                                    btagalgo = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
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
#Chargefilter
process.oppositechargefilter = cms.EDFilter("PATMuonTauOppositeChargeFilter",
                                            muons = cms.InputTag('selectedMuons'),
                                            taus = cms.InputTag('selectedTaus'),
                                            invertChargeFilter = cms.bool(False),
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

# rerun MET for process:
runMetCorAndUncFromMiniAOD(process,
                     isData=False,
                           )


process.out_objectselection = cms.OutputModule("PoolOutputModule",
                                               fileName = cms.untracked.string("fullSelection.root"),
                                               outputCommands = cms.untracked.vstring(
        'keep *_generator_*_*',
        'keep *_selectedMuons_*_*',
        'keep *_selectedTaus_*_*',
        'keep *_selectedJets_*_*',
        'keep *_slimmedMETs_*_*',
        'keep *_kinfit_*_*',
        'keep *_nEventsTotal_*_*',
        'keep *_nEventsSelected_*_*'
        ),
                                               SelectEvents = cms.untracked.PSet( 
      SelectEvents = cms.vstring("objectselection")
        )
                                               )
  

process.mon_all = cms.EDAnalyzer("SelectionAnalyzer",
                              muons = cms.untracked.InputTag( 'slimmedMuons' ),
                              taus  = cms.untracked.InputTag( 'slimmedTaus' ),
                              jets  = cms.untracked.InputTag( 'slimmedJets' ),
                              ) #2

process.mon_muons = process.mon_all.clone( muons = cms.untracked.InputTag( 'selectedMuons' )) #3
process.mon_taus = process.mon_muons.clone( taus = cms.untracked.InputTag( 'selectedTaus' )) #4
process.mon_jets = process.mon_taus.clone( jets = cms.untracked.InputTag( 'selectedJets' )) #5
process.mon_category = process.mon_jets.clone( ) #6

process.nEventsTotal = cms.EDProducer("EventCountProducer")
process.nEventsSelected = cms.EDProducer("EventCountProducer")

process.objectselection = cms.Path(
    process.nEventsTotal
    *process.mon_all
    *process.selectedMuons
    *process.mon_muons
    *process.selectedTaus
    *process.mon_taus
    *process.selectedJets
    *process.mon_jets
    *process.muonfilter
    *process.taufilter
    *process.jetfilter
    *process.oppositechargefilter
    *process.selectedElectrons
    *process.electronveto
    *process.mon_category
    *process.nEventsSelected
    *process.kinfit
    )


process.e = cms.EndPath(process.out_objectselection)