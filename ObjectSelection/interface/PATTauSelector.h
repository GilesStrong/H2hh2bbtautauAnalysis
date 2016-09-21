// Description: Thinning Tau Collection 
// Implementation: adopted version of edm::ThinningProducer
// Original Author:  Jan Oliver Rieger<oliver.rieger@cern.ch>
// Created:  Wed, 15 Jun 2016 08:16:49 GMT

// selector for hadronic tau as defined by
// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#Baseline_mu_tau_h
//
// this selector implements the recommended tauID (Hadron Plus Strip, HPS combined isolation algorithm).

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/stream/ThinningProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

namespace edm {
  class ParameterSetDescription;
}

namespace {
  class PATTauSelector {
  public:
    PATTauSelector(edm::ParameterSet const& pset, edm::ConsumesCollector& cc);
    static void fillDescription(edm::ParameterSetDescription & desc);
    //void preChoose(edm::Handle<pat::TauCollection> taus, edm::Event const& event, edm::EventSetup const& es);
    void preChoose(edm::Event const& event, edm::EventSetup const& es);
    bool choose( unsigned int iIndex, pat::Tau const& iItem);

  private:
    double ptmin_;
    double etamax_;
    double decayModeFinding_;
    double leadTauCandDZ_;
    double iso_;
    std::string isoCrit_;
    std::string electronRejection_;
    std::string muonRejection_;
    
    bool useIso_;
    bool useMuonRejection_;
    bool useElectronRejection_;
    bool invertIso_;
    
    bool debug_;
  };

  PATTauSelector::PATTauSelector(edm::ParameterSet const& pset, edm::ConsumesCollector& cc){
    // Get the values of parameters that you need
    ptmin_ = pset.getParameter<double>("ptmin");
    etamax_ = pset.getParameter<double>("etamax");
    iso_ = pset.getParameter<double>("iso");
    isoCrit_ = pset.getParameter<std::string>("isoCrit");
    decayModeFinding_ = pset.getParameter<double>("decayModeFinding");
    leadTauCandDZ_ = pset.getParameter<double>("leadTauCandDZ");
    electronRejection_ = pset.getParameter<std::string>("electronRejection");
    muonRejection_ = pset.getParameter<std::string>("muonRejection");
    useIso_ = pset.getParameter<bool>("useIso");    
    useMuonRejection_ = pset.getParameter<bool>("useMuonRejection");
    useElectronRejection_ = pset.getParameter<bool>("useElectronRejection");     
    invertIso_ = pset.getParameter<bool>("invertIso");
    debug_ = pset.getParameter<bool>("debug");
  }

  // void PATTauSelector::fillDescription(edm::ParameterSetDescription & desc) {
  //   // Add any additional parameters for the selector here to the ParameterSetDescription.
  //   // For example, if you want to read a track collection in the selector, maybe this    
  //   // desc.add<edm::InputTag>("trackTag");
  //   desc.add<double>("ptmin", 20)->setComment("minimal pT of the particle");
  //   desc.add<double>("etamax", 2.3)->setComment("maximal abs(eta) of the particle");
  //   desc.add<double>("decayModeFinding", 0.5)->setComment("decayModeFinding");
  //   desc.add<double>("leadTauCandDZ", 0.2)->setComment("maximum distance in z direction wrt. the first PV");
  //   desc.add<double>("iso", 0.5)->setComment("isolation threshold");
  //   desc.add<std::string>("isoCrit", "byTightIsolationMVArun2v1DBoldDMwLT")->setComment("isolation criterion");
  //   desc.add<std::string>("electronRejection", "againstElectronVLooseMVA6")->setComment("electron rejection ID");
  //   desc.add<std::string>("muonRejection", "againstMuonTight3")->setComment("muon rejection ID");

  //   desc.add<bool>("useIso", true)->setComment("switch to enable/disable isolation criterion");
  //   desc.add<bool>("useElectronRejection", true)->setComment("switch to enable/disable electron rejection criterion");
  //   desc.add<bool>("useMuonRejection", true)->setComment("switch to enable/disable muon rejection criterion");

  //   desc.add<bool>("debug", false)->setComment("switch to enable/disable debugging output");
  // }

  //void PATTauSelector::preChoose(edm::Handle<pat::TauCollection> taus, edm::Event const& event, edm::EventSetup const& es) {
  void PATTauSelector::preChoose(edm::Event const& event, edm::EventSetup const& es) {
  }

  bool PATTauSelector::choose( unsigned int iIndex, pat::Tau const& tau) {
    // Return true if you select a particular element of the master collection
    // and want it saved in the thinned collection. Otherwise return false.

    bool tau_acc = false;
    bool tau_id = false;
    bool tau_charge = false;
    bool tau_iso = false;
    bool tau_electronrejection = false;
    bool tau_muonrejection = false;

    // tau acceptance
    if(tau.pt() > ptmin_ && fabs(tau.eta()) < etamax_)
      tau_acc = true;

    // tau isolation
    if(useIso_){
          if(invertIso_ && tau.tauID(isoCrit_) < iso_)
            {
            tau_iso = true;
            }
          if(!invertIso_ && tau.tauID(isoCrit_) > iso_)
            {
            tau_iso = true;
            }
    }
    else
      {
	tau_iso=true;
      }
   

    // electron rejection
    if(useElectronRejection_) {
      if(tau.tauID(electronRejection_) > 0.5)
        tau_electronrejection = true;
    } else {
      tau_electronrejection = true;
    }

    // muon rejection
    if(useMuonRejection_) {
      if(tau.tauID(muonRejection_) > 0.5)
        tau_muonrejection = true;
    } else {
      tau_muonrejection = true;
    }

    // tau ID

    // get leading charged had. candidate:
    pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
    if(tau.tauID("decayModeFinding") > decayModeFinding_) {
      if(packedLeadTauCand->dz() < leadTauCandDZ_)
        tau_id = true;
    }

    // require that tau charge is equal to plus or minus 1
    if(fabs(tau.charge() == 1))
      tau_charge = true;

    if(debug_){
      std::cout << tau_acc << " |acceptance: pt=" << tau.pt() << " etaz=" << fabs(tau.eta()) << std::endl;
      std::cout << tau_id << " |id          dmf=" << tau.tauID("decayModeFinding") << ",dz=" << packedLeadTauCand->dz() << std::endl;
      std::cout << tau_charge << " |charge      q=" << tau.charge() << std::endl;
      std::cout << tau_iso << " |isolation:  iso=" << tau.tauID(isoCrit_) << std::endl;
      std::cout << tau_electronrejection << " |el. reject: " << tau.tauID(muonRejection_) << std::endl;
      std::cout << tau_muonrejection << " |mu. reject: " << tau.tauID(muonRejection_) << std::endl;
      std::cout << "-----------------------------------------------------------" << std::endl;
    }

    return(tau_acc && tau_iso && tau_id && tau_charge && tau_electronrejection && tau_muonrejection);
  }
}
