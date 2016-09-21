#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/stream/ThinningProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

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
  class PATElectronSelector {
  public:
    PATElectronSelector(edm::ParameterSet const& pset, edm::ConsumesCollector& cc);
    void preChoose(edm::Event const& event, edm::EventSetup const& es);
    bool choose(unsigned int iIndex, pat::Electron const& iItem);

  private:
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::Handle<reco::VertexCollection> vertices;

    std::string electronID_;

    double ptmin_;
    double etamax_;
    double vtxdxymax_;
    double vtxdzmax_;

    double iso_;
    bool useIso_;
    bool invertIso_;

    bool noMatchedConversions_;
    int maxMissingHits_;

    bool debug_;
  };

  PATElectronSelector::PATElectronSelector(edm::ParameterSet const& pset, edm::ConsumesCollector& cc){
    vtxToken_ = cc.consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertices"));
    electronID_ = pset.getParameter<std::string>("electronID");
    ptmin_ = pset.getParameter<double>("ptmin");
    etamax_ = pset.getParameter<double>("etamax");
    vtxdxymax_ = pset.getParameter<double>("vtxdxymax");
    vtxdzmax_ = pset.getParameter<double>("vtxdzmax");
    iso_ = pset.getParameter<double>("iso");
    useIso_ = pset.getParameter<bool>("useIso");
    invertIso_ = pset.getParameter<bool>("invertIso");
    debug_ = pset.getParameter<bool>("debug");
    noMatchedConversions_ = pset.getParameter<bool>("noMatchedConversions");
    maxMissingHits_ = pset.getParameter<int>("maxMissingHits");
  }


  void PATElectronSelector::preChoose(edm::Event const& event, edm::EventSetup const& es) {
    //Get VertexCollection
    event.getByToken(vtxToken_, vertices);
  }

  bool PATElectronSelector::choose( unsigned int iIndex, pat::Electron const& electron) {
    // Return true if you select a particular element of the master collection
    // and want it saved in the thinned collection. Otherwise return false.

    bool electron_pv = false;
    bool electron_acc = false;
    bool electron_iso = false;
    bool electron_id = false;
    bool electron_matchedConversions = false;
    bool electron_missingHits = false;

    const reco::Vertex &PV = vertices->front();
    if (vertices->empty())
      return false;

    if(electron.electronID(electronID_))
      electron_id = true;

    // electron PV
    if(fabs(electron.gsfTrack()->dxy(PV.position())) < vtxdxymax_ && fabs(electron.gsfTrack()->dz(PV.position()))  < vtxdzmax_) 
      electron_pv = true;

    // electron acceptance
    if(electron.pt() > ptmin_ && fabs(electron.eta()) < etamax_)
      electron_acc = true;

    // electron isolation
    double candIsolation;
    if(useIso_){
      candIsolation = (electron.pfIsolationVariables().sumChargedHadronPt + std::max(
           electron.pfIsolationVariables().sumNeutralHadronEt +
           electron.pfIsolationVariables().sumPhotonEt - 
           0.5 * electron.pfIsolationVariables().sumPUPt, 0.0)) / electron.pt();

      if (!invertIso_ && candIsolation < iso_)
         electron_iso = true;
      if (invertIso_ && candIsolation > iso_)
         electron_iso = true;
    }
    else{
         electron_iso = true;
    }
    
    // additional ID requirements:
    if (noMatchedConversions_ && electron.passConversionVeto())
      electron_matchedConversions = true;
    if (electron.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) <= maxMissingHits_)
      electron_missingHits = true;
   
    // debug output
    if(debug_){
      std::cout << electron_pv  << " |vertex:         dxy=" << fabs(electron.gsfTrack()->dxy(PV.position())) << " dz=" <<fabs(electron.gsfTrack()->dz(PV.position())) << std::endl;
      std::cout << electron_acc << " |acceptance:     pt=" << electron.pt() << " etaz=" << fabs(electron.eta()) << std::endl;
      std::cout << electron_iso << " |isolation:      " << candIsolation << std::endl;
      std::cout << electron_id << " |electronID      " << std::endl;
      std::cout << electron_matchedConversions << " |matched conv.   " << electron.passConversionVeto() << std::endl;
      std::cout << electron_missingHits << " |missing hits:   " << electron.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) << std::endl;
      std::cout << "-----------------------------------------------------------" << std::endl;
    }
    
    return electron_pv && electron_acc && electron_iso && electron_id
           && electron_matchedConversions  && electron_missingHits;
  }
}
