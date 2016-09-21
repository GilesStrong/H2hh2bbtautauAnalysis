// -*- C++ -*-
//
// Package:    H2hh2bbtautauAnalysis/LeptonScaleFactor
// Class:      LeptonScaleFactor
// 
/**\class LeptonScaleFactorProducer LeptonScaleFactorProducer.cc H2hh2bbtautauAnalysis/ScaleFactors/plugins/LeptonScaleFactorProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Benedikt Vormwald
//         Created:  Fri, 22 Jul 2016 13:16:10 GMT
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"


class LeptonScaleFactorProducer : public edm::stream::EDProducer<> {
public:
  explicit LeptonScaleFactorProducer(const edm::ParameterSet&);
  ~LeptonScaleFactorProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;
  

  edm::EDGetTokenT< edm::View<reco::Candidate> > lepToken_;
  bool setEffDataToOne_;
  bool setEffMCToOne_;
  bool getErr_;
  bool isDebug_;
  ScaleFactor * lepSF_;
};

LeptonScaleFactorProducer::LeptonScaleFactorProducer(const edm::ParameterSet& iConfig)
{
  setEffDataToOne_ = iConfig.getParameter<bool>("setEfficiencyDataToOne");
  setEffMCToOne_   = iConfig.getParameter<bool>("setEfficiencyMCToOne");
  getErr_ = iConfig.getParameter<bool>("getErrors");
  isDebug_ = iConfig.getParameter<bool>("debug");

  produces<double>("ScaleFactorValue").setBranchAlias("ScaleFactorValue");
  produces<double>("EfficiencyDataValue").setBranchAlias("EfficiencyDataValue");
  produces<double>("EfficiencyMCValue").setBranchAlias("EfficiencyMCValue");
  if (getErr_) produces<double>("ScaleFactorError").setBranchAlias("ScaleFactorError");
  if (getErr_) produces<double>("EfficiencyDataError").setBranchAlias("EfficiencyDataError");
  if (getErr_) produces<double>("EfficiencyMCError").setBranchAlias("EfficiencyMCError");
  
  lepToken_ = consumes< edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("src"));
  
  lepSF_ = new ScaleFactor();
  lepSF_->init_ScaleFactor(iConfig.getParameter<std::string>("inputFile"));
}


LeptonScaleFactorProducer::~LeptonScaleFactorProducer(){
}


void
LeptonScaleFactorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< edm::View<reco::Candidate> > leptons;
  iEvent.getByToken(lepToken_,leptons);
  
  double sf=1.0;
  double sferr=0.0;
  double effdata=1.0;
  double effdataerr=0.0;
  double effmc=1.0;
  double effmcerr=0.0;

  int muonindex=0;
  
  for(edm::View<reco::Candidate>::const_iterator lep=leptons->begin(); lep!=leptons->end(); ++lep){
    double thissf=1.0;
    double thissferr=0;
    double thiseffdata=1.0;
    double thiseffdataerr=0;
    double thiseffmc=1.0;
    double thiseffmcerr=0;

    if(!setEffDataToOne_){
      thiseffdata    = lepSF_->get_EfficiencyData(lep->pt(), lep->eta());
      if (getErr_) thiseffdataerr = lepSF_->get_EfficiencyDataError(lep->pt(), lep->eta());
    }
    if(!setEffMCToOne_){
      thiseffmc    = lepSF_->get_EfficiencyMC(lep->pt(), lep->eta());
      if (getErr_) thiseffmcerr = lepSF_->get_EfficiencyMCError(lep->pt(), lep->eta());
    }

    thissf    =  thiseffdata/thiseffmc;
    if (getErr_) thissferr = sqrt(pow(thiseffdataerr/thiseffdata,2)+pow(thiseffmcerr/thiseffmc,2))*thissf;
    
    effdata *= thiseffdata;
    effdataerr += pow(thiseffdataerr/thiseffdata,2); //adding squared relative error
    effmc *= thiseffmc;
    effmcerr += pow(thiseffmcerr/thiseffmc,2); //adding squared relative error
    sf *= thissf;
    sferr += pow(thissferr/thissf,2); //adding squared relative error

    if (isDebug_)
      std::cout << "muon " << ++muonindex << ": pt=" << lep->pt() << " eta=" <<  lep->eta() 
                << ": eff(MC)=" << thiseffmc << "+-" << thiseffmcerr 
                << ", eff(data)=" << thiseffdata << "+-" << thiseffdataerr 
                << ", sf=" << thissf << "+-" << thissferr 
                << ", " << thiseffdata/thiseffmc << std::endl;
  }

  effdataerr=sqrt(effdataerr)*effdata; // absolute error
  effmcerr=sqrt(effmcerr)*effmc; // absolute error
  sferr=sqrt(sferr)*sf;  // absolute erroor

  if (isDebug_)
    std::cout << "  eff(MC)=" << effmc << "+-" << effmcerr 
              << "  eff(data)=" << effdata << "+-" << effdataerr
              << ", sf=" << sf << "+-" << sferr << std::endl;
  
  std::auto_ptr<double> pOutSf(new double(sf));
  std::auto_ptr<double> pOutSfErr(new double(sferr));
  std::auto_ptr<double> pOutEffData(new double(effdata));
  std::auto_ptr<double> pOutEffDataErr(new double(effdataerr));
  std::auto_ptr<double> pOutEffMC(new double(effmc));
  std::auto_ptr<double> pOutEffMCErr(new double(effmcerr));

  iEvent.put(std::move(pOutSf),"ScaleFactorValue");
  iEvent.put(std::move(pOutEffData),"EfficiencyDataValue");
  iEvent.put(std::move(pOutEffMC),"EfficiencyMCValue");
  if (getErr_) iEvent.put(std::move(pOutSfErr),"ScaleFactorError");
  if (getErr_) iEvent.put(std::move(pOutEffDataErr),"EfficiencyDataError");
  if (getErr_) iEvent.put(std::move(pOutEffMCErr),"EfficiencyMCError");
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
LeptonScaleFactorProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
LeptonScaleFactorProducer::endStream() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LeptonScaleFactorProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src")->setComment("lepton input collection");
  desc.add<std::string>("inputFile")->setComment("input root file with lepton scale factors and efficiencies");
  desc.add<bool>("setEfficiencyDataToOne",false)->setComment("flag to set data lepton efficiency to one");
  desc.add<bool>("setEfficiencyMCToOne",false)->setComment("flag to set MC lepton efficiency to one");
  desc.add<bool>("getErrors",false)->setComment("flag to enable output of uncertainties");
  desc.add<bool>("debug", false)->setComment("flag to switch on debug output");
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonScaleFactorProducer);
