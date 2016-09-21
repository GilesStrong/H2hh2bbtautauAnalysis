// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Particle.h"


class PATMuonTriggerMatchSelector {
public:
  PATMuonTriggerMatchSelector (const edm::ParameterSet & cfg, edm::ConsumesCollector& cc)
    : path_(cfg.getParameter<std::string>("triggerpath")),
    ignoreTriggerMatch_(cfg.getParameter<bool>("ignoreTriggerMatch")){
  }

  void preChoose(const edm::Event & ed, const edm::EventSetup & es){
  }

  bool choose( unsigned int iIndex, const pat::Muon & m ) const { 
    if (ignoreTriggerMatch_)
      return (true);
    else
      return (m.triggerObjectMatchByPath(path_)!=0); 
  }

private:
  std::string path_;
  bool ignoreTriggerMatch_;
};

