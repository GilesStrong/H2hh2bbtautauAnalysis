import FWCore.ParameterSet.Config as cms

muonScaleFactor = cms.EDProducer('LeptonScaleFactorProducer',
     src       = cms.InputTag( "selectedMuons" ),
     inputFile = cms.string( "../../../HTT-utilities/LepEffInterface/data/Muon/Run2016BCD/Muon_IdIso0p20_eff.root" ),
     getScaleFactor = cms.bool( False ),
     getEfficiency = cms.bool( True ),
     getErrors = cms.bool( False ),
     isData = cms.bool( False ),
     debug =  cms.bool( False ),
)
