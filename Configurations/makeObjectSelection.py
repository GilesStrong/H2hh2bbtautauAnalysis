import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.PythonUtilities.LumiList as LumiList
import os
from Configuration.AlCa.autoCond import autoCond
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
from PhysicsTools.PatAlgos.tools.tauTools import *
from RecoMET.METPUSubtraction.MVAMETConfiguration_cff import runMVAMET
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

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

#~ options.register ('pileupGeneratedFile', "PUweightProducer/etc/pileupMC.root",
options.register ('pileupGeneratedFile', "default",
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "path to generated PU distribution")

#~ options.register ('pileupDataFile', "PUweightProducer/etc/MyDataPileupHistogram_true.root",
options.register ('pileupDataFile', "default",
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

process.GlobalTag.globaltag = options.globalTag

if options.runOnData:
    process.GlobalTag.globaltag = autoCond['run2_data']
    triggerResultsTag = "TriggerResults::HLT"
else:
    process.GlobalTag.globaltag = autoCond['run2_mc']
    triggerResultsTag = "TriggerResults::HLT"

process.options  = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                       allowUnscheduled = cms.untracked.bool(True),
                                       SkipEvent = cms.untracked.vstring( 'BTagCalibration' ),
                                       )

process.source = cms.Source("PoolSource",
    #~ skipEvents=cms.untracked.uint32(75850),
    fileNames = cms.untracked.vstring(
        'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring16MiniAODv2/GluGluToRadionToHHTo2B2Tau_M-300_narrow_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/80000/02EE1330-5D3B-E611-B5AD-001CC4A7C0A4.root'
        #~ '/store/mc/RunIISpring16MiniAODv2/GluGluToRadionToHHTo2B2Tau_M-600_narrow_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/D04ECFFC-CE3A-E611-BF65-001D09FDD6A5.root'
        #~ '/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/00000/0064B539-803A-E611-BDEA-002590D0B060.root'
        #~ '/store/mc/RunIISpring16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/20000/02652A75-C654-E611-9734-0CC47A4D7686.root'
        #~ '/store/data/Run2016D/SingleMuon/MINIAOD/PromptReco-v2/000/276/315/00000/168C3DE5-F444-E611-A012-02163E014230.root'
        ))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

outputMonitorFile = options.outputFile.split(".root")[0] + "_monitor.root"
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputMonitorFile)
                                   )

process.muonsWithIso = cms.EDProducer("PATMuonIsolationProducer",
                                      muons = cms.InputTag('slimmedMuons'),
                                      isoAlgorithm = cms.int32(7)
                                      )

process.unpackedSelectedPatTrigger = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
                                                    patTriggerObjectsStandAlone = cms.InputTag( 'selectedPatTrigger' ),
                                                    triggerResults              = cms.InputTag( triggerResultsTag )
                                                    )

process.matchedPatMuonTriggerMatchHLTIsoMu22 = cms.EDProducer("PATTriggerMatcherDRLessByR" , 
                                                              src     = cms.InputTag( "muonsWithIso" ),
                                                              matched = cms.InputTag( "unpackedSelectedPatTrigger" ),          
                                                              matchedCuts = cms.string( 'type( "TriggerMuon" ) && path( "HLT_IsoMu22_v* || HLT_IsoTkMu22_v*" )' ),
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
                                      triggerpath = cms.string('HLT_IsoMu22_v* || HLT_IsoTkMu22_v*'),
                                      ignoreTriggerMatch = cms.bool(not options.runOnData),    ########  just a HOT FIX #########
                                      filter = cms.bool(options.runOnData)                  ########  just a HOT FIX #########
                                      )

process.selectedIDMuons = cms.EDFilter("PATMuonCloneSelectorFilter",
                                     src = cms.InputTag('triggeredMuons'),
                                     vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                     ptmin = cms.double(0),
                                     etamax = cms.double(100),   
                                     vtxdxymax = cms.double(100),                 
                                     vtxdzmax = cms.double(100),
                                     takeMuonID = cms.string("medium"),
                                     iso = cms.double(0.15),
                                     isoAlgorithm = cms.int32(7),
                                     useIso = cms.bool(False),
                                     invertIso = cms.bool(False),
                                     debug = cms.bool(False),
                                     filter = cms.bool(True)
                                     )

process.selectedIsoMuons = cms.EDFilter("PATMuonCloneSelectorFilter",
                                     src = cms.InputTag('selectedIDMuons'),
                                     vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                     ptmin = cms.double(0),
                                     etamax = cms.double(100),   
                                     vtxdxymax = cms.double(100),                 
                                     vtxdzmax = cms.double(100),
                                     takeMuonID = cms.string("medium"),
                                     iso = cms.double(0.15),
                                     isoAlgorithm = cms.int32(7),
                                     useIso = cms.bool(True),
                                     invertIso = cms.bool(False),
                                     debug = cms.bool(False),
                                     filter = cms.bool(True)
                                     )

process.selectedVtxMuons = cms.EDFilter("PATMuonCloneSelectorFilter",
                                     src = cms.InputTag('selectedIsoMuons'),
                                     vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                     ptmin = cms.double(0),
                                     etamax = cms.double(100),   
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

process.selectedEtaMuons = cms.EDFilter("PATMuonCloneSelectorFilter",
                                     src = cms.InputTag('selectedVtxMuons'),
                                     vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                     ptmin = cms.double(0),
                                     etamax = cms.double(2.4),   
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

process.selectedMuons = cms.EDFilter("PATMuonCloneSelectorFilter",
                                     src = cms.InputTag('selectedEtaMuons'),
                                     #~ src = cms.InputTag('muonsWithIso'),
                                     vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                     ptmin = cms.double(22),
                                     etamax = cms.double(2.4),   
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
                                    etamax = cms.double(2.4),
                                    decayModeFinding = cms.double(0.5),
                                    leadTauCandDZ = cms.double(0.2),
                                    iso = cms.double(0.5),
                                    isoCrit = cms.string("byTightIsolationMVArun2v1DBdR03oldDMwLT"),
                                    electronRejection = cms.string("againstElectronVLooseMVA6"),
                                    muonRejection = cms.string("againstMuonTight3"),
                                    useIso = cms.bool(True),
                                    useElectronRejection = cms.bool(True),
                                    useMuonRejection = cms.bool(True),
                                    invertIso = cms.bool(False),
                                    debug = cms.bool(False),
                                    filter = cms.bool(True)
                                    )

process.selectedTausinvertedIso = process.selectedTaus.clone(invertIso = cms.bool(True))

# Update JEC:
# Do not forget 'L2L3Residual' on data!
if options.runOnData:
    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
        labelName = 'UpdatedJEC',
        jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']), 'None')  # Do not forget 'L2L3Residual' on data!
    )
else:
    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
        labelName = 'UpdatedJEC',
        jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')
    )
    
jetCollection = "updatedPatJetsUpdatedJEC"

process.goodJets = cms.EDFilter("PATJetCloneSelectorFilter",
                                src = cms.InputTag('updatedPatJetsUpdatedJEC'),  
                                ptmin = cms.double(20),
                                etamax = cms.double(2.7),   
                                
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
                                
                                usebtag = cms.bool(False),
                                btagptmin = cms.double(20),
                                btagetamax = cms.double(2.7),
                                btagdesc = cms.double(0.0),
                                btagalgo =    cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                                invertbtag = cms.bool(False),
                                
                                debug = cms.bool(False),
                                filter = cms.bool(True)
                               )

process.selectedJets = process.goodJets.clone(src = cms.InputTag('goodJets'), 
                                              usebtag = cms.bool(False))

process.goodJetsInvTauIso = process.goodJets.clone(taus = cms.InputTag('selectedTausinvertedIso'))

process.selectedJetsInvTauIso = process.goodJetsInvTauIso.clone(src = cms.InputTag('goodJetsInvTauIso'), 
                                              usebtag = cms.bool(False))
                     
process.selectedElectrons = cms.EDFilter("PATElectronCloneSelectorFilter",
                                         src = cms.InputTag('slimmedElectrons'),
                                         vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                         electronID = cms.string("mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                                         ptmin = cms.double(24),
                                         etamax = cms.double(2.4),   
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
                                 maxNumber = cms.uint32(9999),
                                 )

process.jetfilterInvTauIso = process.jetfilter.clone(src = cms.InputTag('selectedJetsInvTauIso'))


process.muonveto = cms.EDFilter("PATCandViewCountFilter",
                                src = cms.InputTag('selectedMuons'),
                                minNumber = cms.uint32(0),
                                maxNumber = cms.uint32(0),
                                )

process.taufilter = cms.EDFilter("PATCandViewCountFilter",
                                 src = cms.InputTag('selectedTaus'),
                                 minNumber = cms.uint32(1),
                                 maxNumber = cms.uint32(1),
                                 )

process.taufilterInvTauIso = process.taufilter.clone(src = cms.InputTag('selectedTausinvertedIso'))

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

process.oppositechargefilterInvTauIso = process.oppositechargefilter.clone(taus = cms.InputTag('selectedTausinvertedIso'))

process.samechargefilter = process.oppositechargefilter.clone(invertChargeFilter = cms.bool(True))

process.samechargefilterInvTauIso = process.oppositechargefilter.clone(invertChargeFilter = cms.bool(True),
                                                                       taus = cms.InputTag('selectedTausinvertedIso')
                                                                       )

process.PUWeightProducer = cms.EDProducer("PUWeightProducer",
                                          generatedFile = cms.string(options.pileupGeneratedFile),
                                          dataFile = cms.string(options.pileupDataFile),
                                          GenHistName = cms.string("pileup"),
                                          DataHistName = cms.string("pileup"),
                                          pileupSummaryInfoLabel = cms.string("slimmedAddPileupInfo"),
                                          runOnData = cms.bool(options.runOnData),
                                          debug = cms.bool(False),
                                          )

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

process.bTaggingEffAnalyzerAK8PF = cms.EDAnalyzer('BTaggingEffAnalyzer',
                                                 JetsTag            = cms.InputTag("goodJets"),
                                                 DiscriminatorTag   = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
                                                 DiscriminatorValue = cms.double(0.8),
                                                 PtNBins            = cms.int32(100),
                                                 PtMin              = cms.double(0.),
                                                 PtMax              = cms.double(1000.),
                                                 EtaNBins           = cms.int32(60),
                                                 EtaMin             = cms.double(-3.),
                                                 EtaMax             = cms.double(3.)
)

process.bTaggingEffAnalyzerAK8PF_LS_Iso = process.bTaggingEffAnalyzerAK8PF.clone()

process.bTaggingEffAnalyzerAK8PF_LS_InvIso = process.bTaggingEffAnalyzerAK8PF.clone(JetsTag = cms.InputTag("goodJetsInvTauIso"))

process.bTaggingEffAnalyzerAK8PF_OS_InvIso = process.bTaggingEffAnalyzerAK8PF.clone(JetsTag = cms.InputTag("goodJetsInvTauIso"))

process.bTaggingSF      = cms.EDProducer('BTaggingSFProducer',
                                         JetsTag            = cms.InputTag("selectedJets"),
                                         DiscriminatorTag   = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
                                         DiscriminatorValue = cms.double(0.8),
                                         runOnData          = cms.bool(options.runOnData),
                                         debug              = cms.bool(True),
                                         CSVFile            = cms.string(options.CSVFile),
                                         EffMapFile         = cms.string(options.EffMapFile),
                                         EffMapFolder       = cms.string(options.EffMapFolder),
)

# load QCD-region efficiency maps as well:

process.bTaggingSFLSIso = process.bTaggingSF.clone(
                                                    EffMapFile  = cms.string(options.EffMapFile.split(".root")[0] + "_LS_Iso.root")
                                                  )  
process.bTaggingSFLSInvIso = process.bTaggingSF.clone(
                                                    EffMapFile = cms.string(options.EffMapFile.split(".root")[0] + "_LS_InvIso.root"),
                                                    JetsTag = cms.InputTag("selectedJetsInvTauIso")
                                                  )
process.bTaggingSFOSInvIso = process.bTaggingSF.clone(
                                                    EffMapFile = cms.string(options.EffMapFile.split(".root")[0] + "_OS_InvIso.root"),
                                                    JetsTag = cms.InputTag("selectedJetsInvTauIso")
                                                  )

process.SVFit = cms.EDProducer('SVFitProducer',
                               muons = cms.InputTag( 'selectedMuons' ),
                               taus = cms.InputTag( 'selectedTaus' ),
                               mets = cms.InputTag( 'slimmedMETs' ),
                               # 1-integrateMarkovChain or 2-integrateVEGAS (default)
                               integrateAlgorithm = cms.int32(2),
                               dataInputFile = cms.string(options.svFitInputFile),
                               MetCovMat = cms.InputTag('METSignificance','METCovariance'),
                               )

process.SVFitInvTauIso = process.SVFit.clone(taus = cms.InputTag( 'selectedTausinvertedIso' ))


process.kinfit = cms.EDProducer('HHKinFit2Producer',
                                muons = cms.InputTag('selectedMuons'),
                                taus  = cms.InputTag('selectedTaus'),
                                bjets = cms.InputTag('selectedJets'),
                                met   = cms.InputTag('slimmedMETs'),
                                metCovMat = cms.InputTag('METSignificance','METCovariance'),
                                debug = cms.bool(False),
                                )

process.kinfitInvTauIso = process.kinfit.clone(taus  = cms.InputTag( 'selectedTausinvertedIso' ),
                                               bjets = cms.InputTag('selectedJetsInvTauIso')
                                               )

# rerun MET for process:
runMetCorAndUncFromMiniAOD(process,
                           isData=options.runOnData,
                           )

if options.redoPuppi:
    makePuppiesFromMiniAOD( process );

    runMetCorAndUncFromMiniAOD(process,
                               isData=options.runOnData,
                               pfCandColl=cms.InputTag("puppiForMET"),
                               recoMetFromPFCs=True,
                               reclusterJets=True,
                               jetFlavor="AK4PFPuppi",
                               postfix="Puppi"
                               )

# configure MVA MET for process:
runMVAMET( process, jetCollectionPF = jetCollection)

process.out_objectselection = cms.OutputModule("PoolOutputModule",
                                               fileName = cms.untracked.string(options.outputFile),
                                               outputCommands = cms.untracked.vstring(
        'drop *',
        #~ 'keep *_externalLHEProducer_*_*',
        'keep *_generator_*_*',
        #~ 'keep *_slimmedMuons_*_*',
        #~ 'keep *_slimmedJets_*_*',
        #~ 'keep *_slimmedTaus_*_*',
        'keep *_selectedMuons_*_*',
        'keep *_selectedTaus_*_*',
        'keep *_selectedJets_*_*',
        'keep *_updatedPatJetsUpdatedJEC_*_*',
        'keep *_prunedGenParticles_*_*',
        'keep *_slimmedGenJets_*_*',
        'keep *_nEvents*_*_*',
        'keep *_TriggerResults_*_SKIM',
        'keep *_PUWeightProducer_*_*',
        # ReRunMET:
        'keep *_slimmedMETs_*_*',
        'keep *_METSignificance_*_*',
        'keep *_slimmedMETsNoHF_*_*',
        #~ "keep *_slimmedMETsPuppi_*_*",
        'keep *_PUWeightProducer_*_*',
        'keep *_SVFit_*_*',
        'keep *_muonIDScaleFactor_*_*',
        'keep *_muonTriggerScaleFactor_*_*',
        'keep *_kinfit_*_*',
        'keep *_TriggerResults_*_SKIM',
        'keep *_bTaggingSF_*_*',
        'keep *_offlineSlimmedPrimaryVertices_*_*',
        ),
                                               SelectEvents = cms.untracked.PSet( 
        SelectEvents = cms.vstring("objectselection")
        ),
                                               dataset = cms.untracked.PSet(
                                                                           filterName = cms.untracked.string(''),
        ),
                                               )
  
process.out_LS_Iso_selection = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string(options.outputFile.split(".root")[0] + "_LS_Iso_selection.root"),
                               outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_generator_*_*',
        'keep *_selectedMuons_*_*',
        'keep *_selectedTaus_*_*',
        'keep *_selectedJets_*_*',
        'keep *_updatedPatJetsUpdatedJEC_*_*',
        'keep *_prunedGenParticles_*_*',
        'keep *_nEvents*_*_*',
        'keep *_TriggerResults_*_SKIM',
        'keep *_PUWeightProducer_*_*',
        'keep *_SVFit_*_*',
        'keep *_slimmedMETs_*_*',
        'keep *_METSignificance_*_*',
        'keep *_muonIDScaleFactor_*_*',
        'keep *_muonTriggerScaleFactor_*_*',
        'keep *_kinfit_*_*',
        'keep *_TriggerResults_*_SKIM',
        'keep *_bTaggingSFLSIso_*_*',
        'keep *_offlineSlimmedPrimaryVertices_*_*',
        ),
                               SelectEvents = cms.untracked.PSet( 
        SelectEvents = cms.vstring("LS_Iso_selection")
        ),
                                               dataset = cms.untracked.PSet(
                                                                           filterName = cms.untracked.string(''),
        ),
                               )


process.out_LS_InvIso_selection = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string(options.outputFile.split(".root")[0] + "_LS_InvIso_selection.root"),
                               outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_generator_*_*',
        'keep *_selectedMuons_*_*',
        'keep *_selectedTausinvertedIso_*_*',
        'keep *_selectedJetsInvTauIso_*_*',
        'keep *_prunedGenParticles_*_*',
        'keep *_nEvents*_*_*',
        'keep *_TriggerResults_*_SKIM',
        'keep *_PUWeightProducer_*_*',
        'keep *_SVFitInvTauIso_*_*',
        'keep *_slimmedMETs_*_*',
        'keep *_METSignificance_*_*',
        'keep *_muonIDScaleFactor_*_*',
        'keep *_muonTriggerScaleFactor_*_*',
        'keep *_kinfitInvTauIso_*_*',
        'keep *_TriggerResults_*_SKIM',
        'keep *_bTaggingSFLSInvIso_*_*',
        'keep *_offlineSlimmedPrimaryVertices_*_*',
        ),
                               SelectEvents = cms.untracked.PSet( 
        SelectEvents = cms.vstring("LS_InvIso_selection")
        ),
                                               dataset = cms.untracked.PSet(
                                                                           filterName = cms.untracked.string(''),
        ),
                               )

process.out_OS_InvIso_selection = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string(options.outputFile.split(".root")[0] + "_OS_InvIso_selection.root"),
                               outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_generator_*_*',
        'keep *_selectedMuons_*_*',
        'keep *_selectedTausinvertedIso_*_*',
        'keep *_selectedJetsInvTauIso_*_*',
        'keep *_prunedGenParticles_*_*',
        'keep *_slimmedGenJets_*_*',
        'keep *_nEvents*_*_*',
        'keep *_TriggerResults_*_SKIM',
        'keep *_PUWeightProducer_*_*',
        'keep *_SVFitInvTauIso_*_*',
        'keep *_slimmedMETs_*_*',
        'keep *_METSignificance_*_*',
        'keep *_muonIDScaleFactor_*_*',
        'keep *_muonTriggerScaleFactor_*_*',
        'keep *_kinfitInvTauIso_*_*',
        'keep *_TriggerResults_*_SKIM',
        'keep *_bTaggingSFOSInvIso_*_*',
        'keep *_offlineSlimmedPrimaryVertices_*_*',
        ),
                               SelectEvents = cms.untracked.PSet( 
        SelectEvents = cms.vstring("OS_InvIso_selection")
        ),
                                               dataset = cms.untracked.PSet(
                                                                           filterName = cms.untracked.string(''),
        ),
                               )

process.mon1 = cms.EDAnalyzer("SelectionAnalyzer",
                              muons = cms.untracked.InputTag( 'slimmedMuons' ),
                              taus  = cms.untracked.InputTag( 'slimmedTaus' ),
                              jets  = cms.untracked.InputTag( 'slimmedJets' ),
                              )

process.mon2 = process.mon1.clone( muons = cms.untracked.InputTag( 'allMuons' ))
process.mon3 = process.mon2.clone( muons = cms.untracked.InputTag( 'triggeredMuons' ))
process.mon4 = process.mon3.clone( muons = cms.untracked.InputTag( 'selectedIDMuons' ))
process.mon5 = process.mon4.clone( muons = cms.untracked.InputTag( 'selectedIsoMuons' ))
process.mon6 = process.mon5.clone( muons = cms.untracked.InputTag( 'selectedVtxMuons' ))
process.mon7 = process.mon6.clone( muons = cms.untracked.InputTag( 'selectedEtaMuons' ))
process.mon8 = process.mon7.clone( muons = cms.untracked.InputTag( 'selectedMuons' ))
process.mon9 = process.mon8.clone( taus = cms.untracked.InputTag( 'selectedTaus' ))
process.mon10 = process.mon9.clone( jets = cms.untracked.InputTag( 'selectedJets' ))
process.mon11 = process.mon10.clone( )

process.nEventsTotal = cms.EDProducer("EventCountProducer")
process.nEventsSelected = cms.EDProducer("EventCountProducer")

process.objectselection = cms.Path(
    process.nEventsTotal
    *process.mon1
    *process.muonsWithIso
    *process.unpackedSelectedPatTrigger
    *process.matchedPatMuonTriggerMatchHLTIsoMu22
    *process.allMuons
    *process.mon2
    *process.triggeredMuons
    *process.mon3
    *process.selectedIDMuons
    *process.mon4
    *process.selectedIsoMuons
    *process.mon5
    *process.selectedVtxMuons
    *process.mon6
    *process.selectedEtaMuons
    *process.mon7
    *process.selectedMuons
    *process.mon8
    *process.selectedTaus
    *process.mon9
    *process.goodJets
    *process.bTaggingEffAnalyzerAK8PF
    *process.selectedJets
    *process.mon10
    *process.muonfilter
    *process.taufilter
    *process.jetfilter
    *process.oppositechargefilter
    *process.selectedElectrons
    *process.electronveto
    *process.mon11
    *process.nEventsSelected
    *process.PUWeightProducer
    *process.SVFit
    *process.kinfit
    *process.muonIDScaleFactor
    *process.muonTriggerScaleFactor
    *process.bTaggingSF
    )

process.LS_Iso_selection = cms.Path(
    process.muonsWithIso
    *process.unpackedSelectedPatTrigger
    *process.matchedPatMuonTriggerMatchHLTIsoMu22
    *process.allMuons
    *process.triggeredMuons
    *process.selectedMuons
    *process.selectedTaus
    *process.goodJets
    *process.bTaggingEffAnalyzerAK8PF_LS_Iso
    *process.selectedJets
    *process.muonfilter
    *process.taufilter
    *process.jetfilter
    *process.samechargefilter
    *process.selectedElectrons
    *process.electronveto
    *process.PUWeightProducer
    *process.SVFit
    *process.kinfit
    *process.muonIDScaleFactor
    *process.muonTriggerScaleFactor
    *process.bTaggingSFLSIso
    )

process.LS_InvIso_selection = cms.Path(
    process.muonsWithIso
    *process.unpackedSelectedPatTrigger
    *process.matchedPatMuonTriggerMatchHLTIsoMu22
    *process.allMuons
    *process.triggeredMuons
    *process.selectedMuons
    *process.selectedTausinvertedIso
    *process.goodJetsInvTauIso
    *process.bTaggingEffAnalyzerAK8PF_LS_InvIso
    *process.selectedJetsInvTauIso
    *process.muonfilter
    *process.taufilterInvTauIso
    *process.jetfilterInvTauIso
    *process.samechargefilterInvTauIso
    *process.selectedElectrons
    *process.electronveto
    *process.PUWeightProducer
    *process.SVFitInvTauIso
    *process.kinfitInvTauIso
    *process.muonIDScaleFactor
    *process.muonTriggerScaleFactor
    *process.bTaggingSFLSInvIso
    )

process.OS_InvIso_selection = cms.Path(
    process.muonsWithIso
    *process.unpackedSelectedPatTrigger
    *process.matchedPatMuonTriggerMatchHLTIsoMu22
    *process.allMuons
    *process.triggeredMuons
    *process.selectedMuons
    *process.selectedTausinvertedIso
    *process.goodJetsInvTauIso
    *process.bTaggingEffAnalyzerAK8PF_OS_InvIso
    *process.selectedJetsInvTauIso
    *process.muonfilter
    *process.taufilterInvTauIso
    *process.jetfilterInvTauIso
    *process.oppositechargefilterInvTauIso
    *process.selectedElectrons
    *process.electronveto
    *process.PUWeightProducer
    *process.SVFitInvTauIso
    *process.kinfitInvTauIso
    *process.muonIDScaleFactor
    *process.muonTriggerScaleFactor
    *process.bTaggingSFOSInvIso
    )


process.e = cms.EndPath(process.out_objectselection*process.out_LS_Iso_selection*process.out_LS_InvIso_selection*process.out_OS_InvIso_selection)


