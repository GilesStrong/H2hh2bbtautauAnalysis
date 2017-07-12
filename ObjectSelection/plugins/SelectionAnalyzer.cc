// -*- C++ -*-
//
// Package:    H2hh2bbtautauAnalysis/SelectionAnalyzer
// Class:      SelectionAnalyzer
// 
/**\class SelectionAnalyzer SelectionAnalyzer.cc H2hh2bbtautauAnalysis/SelectionAnalyzer/plugins/SelectionAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Benedikt Vormwald
//         Created:  Tue, 19 Jul 2016 19:24:29 GMT
//
//


// system include files
#include <memory>
#include <map>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TH1.h"

class SelectionAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit SelectionAnalyzer(const edm::ParameterSet&);
      ~SelectionAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  /// check if histogram was booked
  bool booked(const std::string histName) const { return hists_.find(histName.c_str())!=hists_.end(); };
  /// fill histogram if it had been booked before
  void fill(const std::string histName, double value) const { if(booked(histName.c_str())) hists_.find(histName.c_str())->second->Fill(value); };

  // simple map to contain all histograms;
  // histograms are booked in the beginJob()
  // method
  std::map<std::string, TH1F*> hists_;

  // input tags
  edm::EDGetTokenT<edm::View<pat::Muon> > muonsToken_;
  edm::EDGetTokenT<edm::View<pat::Tau> > tausToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetsToken_;
};


SelectionAnalyzer::SelectionAnalyzer(const edm::ParameterSet& iConfig)
  : hists_(),
    muonsToken_(consumes<edm::View<pat::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muons"))),
    tausToken_(consumes<edm::View<pat::Tau> >(iConfig.getUntrackedParameter<edm::InputTag>("taus"))),
    jetsToken_ (consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("jets" )))
{
  usesResource("TFileService");
}


SelectionAnalyzer::~SelectionAnalyzer(){
}


void
SelectionAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::Handle<edm::View<pat::Muon> > muons;
  edm::Handle<edm::View<pat::Tau> > taus;
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(muonsToken_,muons);
  iEvent.getByToken(tausToken_,taus);
  iEvent.getByToken(jetsToken_,jets);


  fill("yield", 0.5);
  fill("muonMult", muons->size());
  fill("tauMult", taus->size());
  fill("jetMult", jets->size());
  if( jets->size()>0 ) fill("jet0Pt", (*jets)[0].pt());
  if( jets->size()>1 ) fill("jet1Pt", (*jets)[1].pt());
  if( jets->size()>2 ) fill("jet2Pt", (*jets)[2].pt());
  if( jets->size()>3 ) fill("jet3Pt", (*jets)[3].pt());
  //for(std::vector<pat::Muon>::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){
  for(size_t i = 0; i < muons->size(); ++i){
    fill("mu_pt", (*muons)[i]->pt());
    fill("mu_eta", (*muons)[i]->eta());
    fill("mu_iso", (*muons)[i]->iso());
    fill("mu_vtxdxy", (*muons)[i]->vtxdxy());
    fill("mu_vtxdz", (*muons)[i]->vtxdz());
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
SelectionAnalyzer::beginJob()
{
  // register to the TFileService
  edm::Service<TFileService> fs;

  // book histograms:
  hists_["yield"   ]=fs->make<TH1F>("yield"   , "event yield"          ,   1, 0.,   1.);
  hists_["muonMult"]=fs->make<TH1F>("muonMult", "muon multiplicity"    ,  10, 0.,  10.);
  hists_["tauMult" ]=fs->make<TH1F>("tauMult",  "tau multiplicity"     ,  10, 0.,  10.);
  hists_["jetMult" ]=fs->make<TH1F>("jetMult" , "jet multiplicity"     ,  15, 0.,  15.);
  hists_["jet0Pt"  ]=fs->make<TH1F>("jet0Pt"  , "1. leading jet pt"    ,  50, 0., 250.);
  hists_["jet1Pt"  ]=fs->make<TH1F>("jet1Pt"  , "1. leading jet pt"    ,  50, 0., 250.);
  hists_["jet2Pt"  ]=fs->make<TH1F>("jet2Pt"  , "1. leading jet pt"    ,  50, 0., 200.);
  hists_["jet3Pt"  ]=fs->make<TH1F>("jet3Pt"  , "1. leading jet pt"    ,  50, 0., 200.);
  hists_["mu_pt"   ]=fs->make<TH1F>("mu_pt"   , "muon pt;pt[GeV];muons",  50 , 0., 100.);
  hists_["mu_eta"  ]=s->make<TH1F>("mu_eta"   , "muon eta;eta;muons"   ,  50  , -5., 5.);
  hists_["mu_iso"  ]=fs->make<TH1F>("mu_iso"  , "muon iso;iso;muons"   ,  10  , 0., 1.);
  hists_["u_vtxdxy"]=fs->make<TH1F>("mu_vtxdxy","muon vtxdxy;vtxdxy;muons",20  , 0., 0.1);
  hists_["mu_vtxdz"]=fs->make<TH1F>("mu_vtxdz", "muon vtxdz;vtxdz;muons"  , 20  , 0., 0.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SelectionAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SelectionAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SelectionAnalyzer);
