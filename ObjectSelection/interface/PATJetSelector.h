// Description: Thinning Pat Collection 
// Implementation: adopted version of edm::ThinningProducer
// Original Author:  Benedikt Vormwald<benedikt.vormwald@cern.ch>
// Created:  Mon, 11 Jul 2016 11:11:11 GMT

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/stream/ThinningProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include <vector>
#include <string>
#include <math.h>

/////////user include////////////
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"


namespace edm {
  class ParameterSetDescription;
}

namespace {
  class PATJetSelector {
  public:
    //    PATJetSelector(edm::ParameterSet const& pset, edm::ConsumesCollector&& cc);
    PATJetSelector(edm::ParameterSet const& pset, edm::ConsumesCollector& cc);
    static void fillDescription(edm::ParameterSetDescription & desc);
    //    void preChoose(edm::Handle<pat::JetCollection> muons, edm::Event const& event, edm::EventSetup const& es);
    void preChoose(const edm::Event & ed, const edm::EventSetup & es);

    bool choose( unsigned int iIndex, pat::Jet const& iItem);

  private:
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;
    edm::EDGetTokenT<pat::TauCollection> tauToken_;
    edm::Handle<pat::MuonCollection> muons;
    edm::Handle<pat::TauCollection> taus;
    bool checkoverlap_muons_;
    bool checkoverlap_taus_;
    double dR_muons_;
    double dR_taus_;

    double ptmin_;
    double etamax_;

    double nhfmax_;
    double nemfmax_;
    int nconstmin_;
    double mufmax_;
    double chfmin_;
    int cmulmin_;
    double cemfmax_;
    double nemfforwmax_;
    int nneutralforwmin_;

    double btagptmin_;
    double btagetamax_;
    double btagdesc_;
    std::string btagalgo_;
    int usebtag_;
    bool invertbtag_;

    bool debug_;
  };

  //  PATJetSelector::PATJetSelector(edm::ParameterSet const& pset, edm::ConsumesCollector&& cc){
  PATJetSelector::PATJetSelector(edm::ParameterSet const& pset, edm::ConsumesCollector& cc){
    //overlap check
    muonToken_ = cc.consumes<pat::MuonCollection>(pset.getParameter<edm::InputTag>("muons"));
    tauToken_ = cc.consumes<pat::TauCollection>(pset.getParameter<edm::InputTag>("taus"));
    checkoverlap_muons_ = pset.getParameter<bool>("checkoverlap_muons");
    checkoverlap_taus_ = pset.getParameter<bool>("checkoverlap_taus");
    dR_muons_ = pset.getParameter<double>("dR_muon");
    dR_taus_ = pset.getParameter<double>("dR_tau");

    //acceptance
    ptmin_ = pset.getParameter<double>("ptmin");
    etamax_ = pset.getParameter<double>("etamax");

    //jet id
    nhfmax_ = pset.getParameter<double>("neutralhadronfractionmax");
    nemfmax_ = pset.getParameter<double>("neutralEMfractionmax");
    nconstmin_ = pset.getParameter<int>("constituentsmin");
    mufmax_ = pset.getParameter<double>("muonfrationmax");
    chfmin_ = pset.getParameter<double>("chargedhadronfractionmin");
    cmulmin_ = pset.getParameter<int>("chargedmultiplicitymin");
    cemfmax_ = pset.getParameter<double>("chargedEMfractionmax");
    nemfforwmax_ = pset.getParameter<double>("neutralEMfractionforwardmax");
    nneutralforwmin_ = pset.getParameter<int>("neutralparticlesforwardmin");

    //b-tagging
    btagptmin_ = pset.getParameter<double>("btagptmin");
    btagetamax_ = pset.getParameter<double>("btagetamax");
    btagdesc_ = pset.getParameter<double>("btagdesc");
    btagalgo_ = pset.getParameter<std::string>("btagalgo");
    usebtag_ = pset.getParameter<bool>("usebtag");
    invertbtag_ = pset.getParameter<bool>("invertbtag");

    //debugging info
    debug_ = pset.getParameter<bool>("debug");
  }

  // void PATJetSelector::fillDescription(edm::ParameterSetDescription & desc) {
  //   //overlap check
  //   desc.add<edm::InputTag>("muons")->setComment("input collection of muons");
  //   desc.add<edm::InputTag>("taus")->setComment("input collection of taus");
  //   desc.add<bool>("checkoverlap_muons", false)->setComment("switch to enable/disable check for muon overlap");
  //   desc.add<bool>("checkoverlap_taus", false)->setComment("switch to enable/disable check for tau overlap");
  //   desc.add<double>("dR_muon", 0.5)->setComment("dR to check for muon overlap");
  //   desc.add<double>("dR_tau", 0.5)->setComment("dR to check for tau overlap");

  //   //acceptance
  //   desc.add<double>("ptmin", 20)->setComment("minimal pT of the particle");
  //   desc.add<double>("etamax", 4.7)->setComment("maximal abs(eta) of the particle");

  //   //jet id
  //   desc.add<double>("neutralhadronfractionmax",0.99)->setComment("maximal neutral hadron fraction");
  //   desc.add<double>("neutralEMfractionmax",0.99)->setComment("maximal neutral EM fraction");
  //   desc.add<int>("constituentsmin",1)->setComment("minimal number of constituents");
  //   desc.add<double>("muonfrationmax",100)->setComment("maximal muon energy fraction");
  //   desc.add<double>("chargedhadronfractionmin",0)->setComment("minimal charged hadron fraction");
  //   desc.add<int>("chargedmultiplicitymin",0)->setComment("minimal charged multiplicity");
  //   desc.add<double>("chargedEMfractionmax",0.99)->setComment("maximal charged EM fraction in forward region");
  //   desc.add<double>("neutralEMfractionforwardmax",0.90)->setComment("maximal neutral EM fraction in forward region");
  //   desc.add<int>("neutralparticlesforwardmin",10)->setComment("minimal number of particles in forward region");

  //   //b-tagging
  //   desc.add<double>("btagdesc", 0.5)->setComment("minimal b-tag discriminator value");
  //   desc.add<std::string>("btagalgo", "pfCombinedInclusiveSecondaryVertexV2BJetTags")->setComment("used btag discriminator algorithm");
  //   desc.add<double>("btagptmin", 20)->setComment("minimal pT of b-tagged jets");
  //   desc.add<double>("btagetamax", 2.4)->setComment("maximal abs(eta) of b-tagged jets");
  //   desc.add<bool>("usebtag", true)->setComment("switch to enable/disable b-tag criterion");
  //   desc.add<bool>("invertbtag", false)->setComment("switch to invert b-tag criterion");

  //   //debugging info
  //   desc.add<bool>("debug", false)->setComment("switch to enable/disable debugging output");
  // }

  //  void PATJetSelector::preChoose(edm::Handle<pat::JetCollection> jets, edm::Event const& event, edm::EventSetup const& es) {
  void PATJetSelector::preChoose(edm::Event const& event, edm::EventSetup const& es) {
    event.getByToken(muonToken_, muons);
    event.getByToken(tauToken_, taus);
  }

  bool PATJetSelector::choose( unsigned int iIndex, pat::Jet const& jet) {
    bool jet_acc = false;
    bool jet_id = false;
    bool jet_btag = false;
    bool jet_nooverlap_muons = false;
    bool jet_nooverlap_taus = false;

    //jet acceptance
    if(jet.pt() > ptmin_ && fabs(jet.eta()) < etamax_)
      jet_acc = true;


    //jet id
    double NHF = jet.neutralHadronEnergyFraction();
    double NEMF = jet.neutralEmEnergyFraction();
    int NumConst = (jet.chargedMultiplicity()+jet.neutralMultiplicity());
    double MUF = jet.muonEnergyFraction();
    double CHF = jet.chargedHadronEnergyFraction();
    int CHM = jet.chargedMultiplicity();
    double CEMF = jet.chargedEmEnergyFraction();
    int NumNeutralParticle = jet.neutralMultiplicity();

    if(fabs(jet.eta())<=3)
      jet_id = (NHF<nhfmax_ && NEMF<nemfmax_ && NumConst>nconstmin_ && MUF<mufmax_) && ((fabs(jet.eta())<=2.4 && CHF>chfmin_ && CHM>cmulmin_ && CEMF<cemfmax_) || fabs(jet.eta())>2.4);
    else
      jet_id = (NEMF<nemfforwmax_ && NumNeutralParticle>nneutralforwmin_);

    //jet b-tagging
    if(usebtag_){
      double btagvalue = jet.bDiscriminator(btagalgo_);
      double btageta = fabs(jet.eta());
      if(jet.pt()>btagptmin_){
        if(invertbtag_){
          if ((btageta > btagetamax_) || (btageta < btagetamax_ && btagvalue < btagdesc_))
            jet_btag=true;
        }
        else{
          if(btageta < btagetamax_ && btagvalue > btagdesc_)
            jet_btag=true;
        }
      }
    }
    else
      jet_btag=true;
    
    //overlap check muons
    if(checkoverlap_muons_){
      jet_nooverlap_muons=true;
      for (pat::MuonCollection::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){
        double deta2 = pow(mu->eta()-jet.eta(),2);
        double dphi2 = pow(mu->phi()-jet.phi(),2);
        double dR = sqrt(deta2+dphi2);
        if (dR<dR_muons_){
          jet_nooverlap_muons=false;
          break;
        }
      }
    }
    else
      jet_nooverlap_muons = true;
    
    //overlap check taus
    if(checkoverlap_taus_){
      jet_nooverlap_taus=true;
      for (pat::TauCollection::const_iterator tau=taus->begin(); tau!=taus->end(); ++tau){
        double deta2 = pow(tau->eta()-jet.eta(),2);
        double dphi2 = pow(tau->phi()-jet.phi(),2);
        double dR = sqrt(deta2+dphi2);
        if (dR<dR_taus_){
          jet_nooverlap_taus=false;
          break;
        }
      }
    }
    else
      jet_nooverlap_taus = true;
  
    //debugging information
    if(debug_){
      std::cout << jet_acc << " |acceptance: pt=" << jet.pt() << " eta=" << fabs(jet.eta()) << std::endl;
      std::cout << jet_id  << " |id" << std::endl;
      if(usebtag_) std::cout << jet_btag << " |btag:  descr="<< jet.bDiscriminator(btagalgo_) << " eta=" << fabs(jet.eta()) << std::endl;
      if(checkoverlap_muons_) std::cout << jet_nooverlap_muons << " |no muon overlap" << std::endl;
      if(checkoverlap_taus_)  std::cout << jet_nooverlap_taus  << " |no tau overlap" << std::endl;
      std::cout << "-----------------------------------------------------------" << std::endl;
    }
        
    return(jet_acc && jet_id && jet_btag && jet_nooverlap_muons && jet_nooverlap_taus);
  }
}

