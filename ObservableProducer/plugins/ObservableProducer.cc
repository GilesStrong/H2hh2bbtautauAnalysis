// -*- C++ -*-
//
// Package:    H2hh2bbtautauAnalysis/ObservableProducer
// Class:      ObservableProducer
// 
/**\class ObservableProducer ObservableProducer.cc H2hh2bbtautauAnalysis/ObservableProducer/plugins/ObservableProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jan Oliver Rieger
//         Created:  Thu, 28 Jul 2016 09:02:07 GMT
//
//
// system include files
#include <memory>
#include <vector>
#include <math.h>
#include <exception>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"


class ObservableProducer : public edm::stream::EDProducer<> {
   public:
      explicit ObservableProducer(const edm::ParameterSet&);
      ~ObservableProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  edm::EDGetTokenT<pat::METCollection> metsToken_;

  std::vector<math::XYZTLorentzVector> muon_p4;
  std::vector<math::XYZTLorentzVector> met_p4;

  double transmass_mu;
      
};


ObservableProducer::ObservableProducer(const edm::ParameterSet& iConfig)
{
  produces<double>(); 
  muonsToken_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  metsToken_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
}


ObservableProducer::~ObservableProducer()
{
}

void
ObservableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonsToken_,muons);
  
  Handle<pat::METCollection> mets;
  iEvent.getByToken(metsToken_,mets);
  
  std::auto_ptr<double> output(new double(transmass_mu));


  ///////////////////////conversion to lorentzvector/////////////////////////////////////////////////
  //muon
  for (pat::MuonCollection::const_iterator muon = muons->begin(); muon != muons->end(); ++muon)
    {
     muon_p4.push_back(muon->p4());
    }
  //met
  pat::MET met = mets->at(0);
  met_p4.push_back(met.corP4());
  ///////////////////////////////////////////////////////////////////////////////////////////////////
 
  //transverse mass for leading muon
  transmass_mu = sqrt(2 * muon_p4[0].Pt() * met_p4[0].Pt() * (1-cos(muon_p4[0].Phi()-met_p4[0].Phi())) );

  iEvent.put(output);

  //std::cout<< "transverse mass: " << transmass_mu <<std::endl;

  muon_p4.clear();
  met_p4.clear();
}

void ObservableProducer::beginStream(edm::StreamID){}

void ObservableProducer::endStream(){}


 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ObservableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ObservableProducer);
