import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

#process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'file:///pnfs/desy.de/cms/tier2/store/mc/RunIISpring16MiniAODv2/GluGluToRadionToHHTo2B2Tau_M-300_narrow_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/80000/02EE1330-5D3B-E611-B5AD-001CC4A7C0A4.root'
        'file:///afs/desy.de/user/r/riegerjo/public/CMSSW_8_0_10/src/H2hh2bbtautauAnalysis/ReRunMET/corMETMiniAOD.root'
        )
                            )




process.myProducerLabel = cms.EDProducer('SVFitProducer',
                                         muons = cms.InputTag( 'slimmedMuons' ),
                                         taus = cms.InputTag( 'slimmedTaus' ),
                                         mets = cms.InputTag( 'slimmedMETs' ),
                                         MetCovMat = cms.InputTag('METSignificance','METCovariance'), 
                                         integrateAlgorithm = cms.int32(2)         ## 1-integrateMarkovChain or 2-integrateVEGAS (default)
                                         )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('myOutputFile.root')
                               )


process.p = cms.Path(process.myProducerLabel)

process.e = cms.EndPath(process.out)
