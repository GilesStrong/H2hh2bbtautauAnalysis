// Description: Thinning Muon Collection 
// Implementation: adopted version of edm::ThinningProducer
// Original Author:  Jan Oliver Rieger<oliver.rieger@cern.ch>
// Created:  Wed, 15 Jun 2016 08:16:49 GMT

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/stream/ThinningProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include <vector>
#include <string>
#include <math.h>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

namespace edm {
  class ParameterSetDescription;
}

namespace {
  class PATMuonSelector {
  public:
    PATMuonSelector(edm::ParameterSet const& pset, edm::ConsumesCollector& cc);
    //static void fillDescription(edm::ParameterSetDescription & desc);
    //void preChoose(edm::Handle<pat::MuonCollection> muons, edm::Event const& event, edm::EventSetup const& es);
    void preChoose(edm::Event const& event, edm::EventSetup const& es);
    bool choose(unsigned int iIndex, pat::Muon const& iItem);

  private:
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::Handle<reco::VertexCollection> vertices;

    double ptmin_;
    double etamax_;
    double vtxdxymax_;
    double vtxdzmax_;

    double iso_;
    int isoAlgo_;
    bool useIso_;
    bool invertIso_;

    std::string takeMuonID_;

    bool debug_;
  };

  PATMuonSelector::PATMuonSelector(edm::ParameterSet const& pset, edm::ConsumesCollector& cc){
    // Get the values of parameters that you need
    vtxToken_ = cc.consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertices"));
    ptmin_ = pset.getParameter<double>("ptmin");
    etamax_ = pset.getParameter<double>("etamax");
    vtxdxymax_ = pset.getParameter<double>("vtxdxymax");
    vtxdzmax_ = pset.getParameter<double>("vtxdzmax");
    takeMuonID_ = pset.getParameter<std::string>("takeMuonID");
    iso_ = pset.getParameter<double>("iso");
    isoAlgo_ = pset.getParameter<int>("isoAlgorithm");
    useIso_ = pset.getParameter<bool>("useIso");
    invertIso_ = pset.getParameter<bool>("invertIso");
    debug_ = pset.getParameter<bool>("debug");
  }

  // void PATMuonSelector::fillDescription(edm::ParameterSetDescription & desc) {
  //   // Add any additional parameters for the selector here to the ParameterSetDescription.
  //   // For example, if you want to read a track collection in the selector, maybe this    
  //   // desc.add<edm::InputTag>("trackTag");
  //   desc.add<edm::InputTag>("vertices")->setComment("input collection of primary vertices");
  //   desc.add<double>("ptmin", 19)->setComment("minimal pT of the particle");
  //   desc.add<double>("etamax", 2.1)->setComment("maximal abs(eta) of the particle");
  //   desc.add<double>("vtxdxymax", 0.045)->setComment("maximal distance of the particle from primary vertex in x-y plane");
  //   desc.add<double>("vtxdzmax", 0.2)->setComment("maximal distance of the particle from primary vertex in z direction");
  //   desc.add<std::string>("takeMuonID", "medium")->setComment("used particle ID working point [loose, medium, tight]");
  //   desc.add<double>("iso", 0.15)->setComment("isolation threshold");
  //   desc.add<int>("isoAlgorithm", 7)->setComment("used isolation value: TrackIso=0, EcalIso=1, HcalIso=2, PfAllParticleIso=3,PfChargedHadronIso=4, PfNeutralHadronIso=5, PfGammaIso=6, User1Iso=7, User2Iso=8, User3Iso=9, User4Iso=10, User5Iso=11, UserBaseIso=7, CaloIso=-1, PfPUChargedHadronIso=12, PfChargedAllIso=13 (defined in pat::IsolationKeys)");
  //   desc.add<bool>("useIso", true)->setComment("switch to enable/disable isolation criterion");
  //   desc.add<bool>("invertIso", false)->setComment("switch to invert isolation criterion");
  //   desc.add<bool>("debug", false)->setComment("switch to enable/disable debugging output");
  // }

  //void PATMuonSelector::preChoose(edm::Handle<pat::MuonCollection> muons, edm::Event const& event, edm::EventSetup const& es) {
  void PATMuonSelector::preChoose(edm::Event const& event, edm::EventSetup const& es) {
    //Get VertexCollection
    event.getByToken(vtxToken_, vertices);
  }

  bool PATMuonSelector::choose( unsigned int iIndex, pat::Muon const& muon) {
    // Return true if you select a particular element of the master collection
    // and want it saved in the thinned collection. Otherwise return false.

    bool muon_iso = false;
    bool muon_id = false;
    bool muon_pv = false;
    bool muon_acc = false;

    const reco::Vertex &PV = vertices->front();
    //if (vertices->empty()) return; // skip the event if no PV found
    
    //Muon ID "loose","medium" and "tight"
    if(takeMuonID_=="medium" && muon.isMediumMuon()) 
      muon_id = true;
    if(takeMuonID_=="loose" && muon.isLooseMuon()) 
      muon_id = true;
    if(takeMuonID_=="tight" && muon.isTightMuon(PV))
      muon_id = true; 
    if(takeMuonID_!="loose" && takeMuonID_!="medium" && takeMuonID_!="tight")
      std::cout << "Muon ID has to be 'loose', 'medium', or 'tight'" << std::endl;

    //Muon PV
    if(fabs(muon.muonBestTrack()->dxy(PV.position())) < vtxdxymax_ && fabs(muon.muonBestTrack()->dz(PV.position()))  < vtxdzmax_) 
      muon_pv = true; 

    //Muon acceptance
    if(muon.pt() > ptmin_ && fabs(muon.eta()) < etamax_)
      muon_acc = true;

    //Muon isolation
    if(useIso_){
      double isovalue = muon.userIsolation(pat::IsolationKeys(isoAlgo_));
      if(invertIso_  && isovalue > iso_)
          muon_iso=true;
      if(!invertIso_ && isovalue < iso_)
          muon_iso=true;
      if(isovalue < 0){
        std::cout << "muon isolation not set correctly!" << std::endl;
        muon_iso=false;
      }
    }
    else{
      muon_iso=true;
    }
   
    if(debug_){
      std::cout << muon_id  << " |ID:         muonID=" << muon.isMediumMuon() << std::endl;
      std::cout << muon_pv  << " |vertex:     dxy=" << fabs(muon.muonBestTrack()->dxy(PV.position())) << " dz=" <<fabs(muon.muonBestTrack()->dz(PV.position())) << std::endl;
      std::cout << muon_acc << " |acceptance: pt=" << muon.pt() << " etaz=" << fabs(muon.eta()) << std::endl;
      if(useIso_) std::cout << muon_iso << " |isolation:  iso="<< muon.userIsolation(pat::IsolationKeys(isoAlgo_)) << std::endl;
      std::cout << "-----------------------------------------------------------" << std::endl;
    }
    
    
    return(muon_id && muon_pv && muon_acc && muon_iso);
  }
}
