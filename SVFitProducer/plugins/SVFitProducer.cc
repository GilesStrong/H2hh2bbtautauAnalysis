// -*- C++ -*-
//
// Package:    H2hh2bbtautauAnalysis/SVFitProducer
// Class:      SVFitProducer
// 
// Description: [one line class summary]
//
// Original Author:  Jan Oliver Rieger
//         Created:  Tue, 02 Aug 2016 10:26:33 GMT
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
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/METReco/interface/MET.h"
#include "RecoMET/METAlgorithms/interface/METSignificance.h"

#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
#include "TMatrixD.h" 

class SVFitProducer : public edm::stream::EDProducer<> {
public:
  explicit SVFitProducer(const edm::ParameterSet&);
  ~SVFitProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;
  
  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  edm::EDGetTokenT<pat::TauCollection> tauonsToken_;
  edm::EDGetTokenT<pat::METCollection> metsToken_;
  edm::EDGetTokenT<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> > > metcovmatToken_;


  std::vector<math::XYZTLorentzVector> muon_p4;
  std::vector<math::XYZTLorentzVector> tau_p4;
  std::vector<math::XYZTLorentzVector> met_p4;

  std::vector<int> DecayMode;
  
  double sigmaxx;
  double sigmaxy;
  double sigmayy;
  
  double svFitMass;
  
  int integrateAlgo;
  std::string dataInputFile_;
};

SVFitProducer::SVFitProducer(const edm::ParameterSet& iConfig)
{
  produces<double>(); 
  muonsToken_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  tauonsToken_ = consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"));
  metsToken_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  metcovmatToken_ = consumes<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> > >(iConfig.getParameter<edm::InputTag>("MetCovMat")); 
  integrateAlgo = iConfig.getParameter<int>("integrateAlgorithm");
  dataInputFile_  = iConfig.getParameter<std::string>("dataInputFile");
}

SVFitProducer::~SVFitProducer(){}

void
SVFitProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonsToken_,muons);

  Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauonsToken_,taus);
    
  Handle<pat::METCollection> mets;
  iEvent.getByToken(metsToken_,mets);

  Handle<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> > > MetCovMats;
  iEvent.getByToken(metcovmatToken_, MetCovMats); 

  std::auto_ptr<double> output_mass(new double(svFitMass));
  
  ///////////////////////conversion to lorentzvector/////////////////////////////////////////////////
   
  for (pat::MuonCollection::const_iterator muon = muons->begin(); muon != muons->end(); ++muon)
    {
      muon_p4.push_back(muon->p4());
    }
  
  for (pat::TauCollection::const_iterator tau = taus->begin(); tau != taus->end(); ++tau)
    {
      tau_p4.push_back(tau->p4());
      DecayMode.push_back(tau->decayMode());
    }
  
  for (pat::METCollection::const_iterator met = mets->begin(); met != mets->end(); ++met)
    {
      met_p4.push_back(met->corP4());
      // sigmaxx= met->getSignificanceMatrix()(0,0);
      // sigmaxy= met->getSignificanceMatrix()(1,0);
      // sigmayy= met->getSignificanceMatrix()(1,1);
    }

  //Entries of CovMat
  sigmaxx = (*MetCovMats)(0,0);
  sigmaxy = (*MetCovMats)(1,0);
  sigmayy = (*MetCovMats)(1,1);
   
  //check
  if(met_p4.empty() || tau_p4.empty() || muon_p4.empty()) return;
  //std::cout << "covMET: sigmaxx=" << sigmaxx << ", sigmaxy=" << sigmaxy << ", sigmayy=" << sigmayy << std::endl;

  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //########################################Begin:SVFIT###################################################
  ////////////////////////////////////////////////////////////////////////////////////////////////////////  
  
  TMatrixD covMET(2, 2);
  covMET[0][0] = sigmaxx;
  covMET[1][0] = sigmaxy;
  covMET[0][1] = sigmaxy; //keep in mind: CovMET is a symmetric Matrix!
  covMET[1][1] = sigmayy;

  double metPx =  met_p4[0].Px();
  double metPy =  met_p4[0].Py();
  svFitStandalone::Vector MET(metPx, metPy, 0.);
    
  // define lepton four vectors
  //l1 is hadronic tau
  double l1Pt = tau_p4[0].Pt();
  double l1Eta = tau_p4[0].Eta();
  double l1Phi = tau_p4[0].Phi();
  double l1Px   =  tau_p4[0].Px();
  double l1Py   =  tau_p4[0].Py();
  double l1Pz   =  tau_p4[0].Pz();
  double l1Mass =  tau_p4[0].M();
  double l1En   = TMath::Sqrt(l1Px*l1Px + l1Py*l1Py + l1Pz*l1Pz + l1Mass*l1Mass);
  svFitStandalone::LorentzVector l1(l1Px, l1Py, l1Pz, l1En); // Tau_had
      
  //l2 is lepton from tau decay. here: mu
  double l2Pt = muon_p4[0].Pt();
  double l2Eta = muon_p4[0].Eta();
  double l2Phi = muon_p4[0].Phi();
  double l2Px   =   muon_p4[0].Px();
  double l2Py   =   muon_p4[0].Py();
  double l2Pz   =   muon_p4[0].Pz();
  double l2Mass =   105.7e-3;
  double l2En   = TMath::Sqrt(l2Px*l2Px + l2Py*l2Py + l2Pz*l2Pz + l2Mass*l2Mass);
  svFitStandalone::LorentzVector l2(l2Px, l2Py, l2Pz, l2En); // Tau_mu
    
  std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
   
  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay,  l1Pt, l1Eta, l1Phi, l1Mass, DecayMode[0])); // tau -> 1prong0pi0 hadronic decay (Pt, eta, phi, mass, pat::Tau.decayMode())
  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, l2Pt, l2Eta,  l2Phi,  l2Mass)); // tau -> muon decay (Pt, eta, phi, mass) --->>> kTauToMuDecay or kTauToElecDecay
  SVfitStandaloneAlgorithm algo(measuredTauLeptons, metPx, metPy, covMET, 0);
  algo.addLogM(false);

  edm::FileInPath inputFileName_visPtResolution(dataInputFile_);
  TH1::AddDirectory(false);  
  TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
  algo.shiftVisPt(true, inputFile_visPtResolution);
  
  //Choose IntegrateAlgorithm
  if(integrateAlgo==1)
    {
      algo.integrateMarkovChain();
    }
  else if(integrateAlgo==2)
    { 
      algo.integrateVEGAS();
    }

  svFitMass = algo.getMass();
  
  if ( algo.isValidSolution() ) {
    std::cout << "found mass = " << svFitMass << std::endl;
  } else {
    std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //########################################End:SVFIT#####################################################
  //////////////////////////////////////////////////////////////////////////////////////////////////////// 

  delete inputFile_visPtResolution;

  iEvent.put(output_mass);
  
  muon_p4.clear();
  tau_p4.clear();
  met_p4.clear();
  DecayMode.clear();
}

void SVFitProducer::beginStream(edm::StreamID){}
void SVFitProducer::endStream(){}

void
SVFitProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  
  //desc.setUnknown();
  desc.add<edm::InputTag>("muons");
  desc.add<edm::InputTag>("taus");
  desc.add<edm::InputTag>("mets");
  desc.add<edm::InputTag>("MetCovMat")->setComment("used CovarianceMatrix");
  desc.add<int>("integrateAlgorithm", 2)->setComment("Algorithm: integrateMarkovChain (1) / integrateVEGAS (2)");
  desc.add<std::string>("dataInputFile", "TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root")->setComment("data input file");
  descriptions.addDefault(desc);
  }

DEFINE_FWK_MODULE(SVFitProducer);
