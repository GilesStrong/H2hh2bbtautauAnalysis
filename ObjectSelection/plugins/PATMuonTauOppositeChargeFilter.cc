#include <memory>
#include <string>
#include <vector>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

#include "DataFormats/Candidate/interface/OverlapChecker.h"

class  PATMuonTauOppositeChargeFilter : public edm::stream::EDFilter<> {
   public:
      explicit  PATMuonTauOppositeChargeFilter(const edm::ParameterSet&);
      ~ PATMuonTauOppositeChargeFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      
  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  edm::EDGetTokenT<pat::TauCollection> tausToken_;
      
  bool debug_;
  bool invertChargeFilter_;
  
};

 PATMuonTauOppositeChargeFilter:: PATMuonTauOppositeChargeFilter(const edm::ParameterSet& iConfig)
{
  // constructor
   
   muonsToken_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
   tausToken_ = consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"));

   debug_ = iConfig.getParameter<bool>("debug");
   invertChargeFilter_ = iConfig.getParameter<bool>("invertChargeFilter");
}


 PATMuonTauOppositeChargeFilter::~ PATMuonTauOppositeChargeFilter()
{
}


bool  PATMuonTauOppositeChargeFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<pat::MuonCollection> muons;
   Handle<pat::TauCollection> taus;
   iEvent.getByToken(muonsToken_,muons);
   iEvent.getByToken(tausToken_,taus);

   bool oppositeCharged = false;
   bool sameCharged = false;
   bool overlapCheck = false;

   OverlapChecker overlap;

   if((int32_t) muons->size() == 1) 
     {
       if((int32_t) taus->size() == 1) 
	 {
	   pat::Muon selMu = muons->front();
	   pat::Tau selTau = taus->front();
	   if(selMu.charge()*selTau.charge() == -1) 
	     {
	       oppositeCharged = true;
	     }
	   else if(selMu.charge()*selTau.charge() == 1)
	     {
	       sameCharged = true;
	     }
	   if (!overlap(selMu, selTau)) 
	     {
	       overlapCheck = true;
	     }
	   
	 }
     }
   
   if(debug_){
	   std::cout << "-----------------------------------------------------------" << std::endl;
	   std::cout << "   number of selected muons: " << muons->size() << std::endl;
	   std::cout << "   number of selected taus: " << taus->size() << std::endl;
	   std::cout << "   opposite charged leptons: " << oppositeCharged << std::endl;
	   std::cout << "   same  charged leptons: " << sameCharged << std::endl;
	   std::cout << "-----------------------------------------------------------" << std::endl;
   }

   if(!invertChargeFilter_) return (oppositeCharged && overlapCheck);
   return (sameCharged && overlapCheck); 
}

void  PATMuonTauOppositeChargeFilter::beginStream(edm::StreamID) {
}

void  PATMuonTauOppositeChargeFilter::endStream() {
}

void  PATMuonTauOppositeChargeFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<bool>("invertChargeFilter", false);
  desc.add<edm::InputTag>("muons");
  desc.add<edm::InputTag>("taus");
  desc.add<bool>("debug", false);                      
  desc.add<bool>("filter", true);           
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATMuonTauOppositeChargeFilter);
