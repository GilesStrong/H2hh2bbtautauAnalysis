// -*- C++ -*-
//
// Package:    H2hh2bbtautauAnalysis/PATMuonIsolationProducer
// Class:      PATMuonIsolationProducer
// 
/**\class PATMuonIsolationProducer PATMuonIsolationProducer.cc H2hh2bbtautauAnalysis/PATMuonIsolationProducer/plugins/PATMuonIsolationProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Benedikt Vormwald
//         Created:  Thu, 07 Jul 2016 15:55:50 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

class PATMuonIsolationProducer : public edm::stream::EDProducer<> {
public:
  explicit PATMuonIsolationProducer(const edm::ParameterSet&);
  ~PATMuonIsolationProducer() {};
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginStream(edm::StreamID) override {};
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override {};

  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  int isoAlgo_;
};

PATMuonIsolationProducer::PATMuonIsolationProducer(const edm::ParameterSet& iConfig){
  produces<pat::MuonCollection>(); 
  muonsToken_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  isoAlgo_ = iConfig.getParameter<int>("isoAlgorithm");
}

void
PATMuonIsolationProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonsToken_,muons);
  std::unique_ptr<pat::MuonCollection> muonsCopy(new pat::MuonCollection(*muons));

  for(std::vector<pat::Muon>::iterator mu=muonsCopy->begin(); mu!=muonsCopy->end(); ++mu){
    float isolation = (mu->pfIsolationR04().sumChargedHadronPt + std::max(0., mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5*mu->pfIsolationR04().sumPUPt))/mu->pt();
    mu->setIsolation(pat::IsolationKeys(isoAlgo_),isolation);
  }


  iEvent.put(std::move(muonsCopy));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PATMuonIsolationProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons");
  desc.add<int>("isoAlgorithm", 7)->setComment("used isolation value: TrackIso=0, EcalIso=1, HcalIso=2, PfAllParticleIso=3,PfChargedHadronIso=4, PfNeutralHadronIso=5, PfGammaIso=6, User1Iso=7, User2Iso=8, User3Iso=9, User4Iso=10, User5Iso=11, UserBaseIso=7, CaloIso=-1, PfPUChargedHadronIso=12, PfChargedAllIso=13 (defined in pat::IsolationKeys)");
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATMuonIsolationProducer);
