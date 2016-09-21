// Original Author:  Dinko Ferencek
// modified by: Viktor Kutzner

// produces p_T/eta maps for all jets and jets which pass the discriminator value threshold set in the pyset configuration.

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH2D.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

typedef std::vector<pat::Jet> PatJetCollection;

class BTaggingEffAnalyzer : public edm::EDAnalyzer {
   public:
      explicit BTaggingEffAnalyzer(const edm::ParameterSet&);
      ~BTaggingEffAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      const edm::InputTag jetsTag;
      const std::string   discriminatorTag;
      const double  discriminatorValue;
      const int     ptNBins;
      const double  ptMin;
      const double  ptMax;
      const int     etaNBins;
      const double  etaMin;
      const double  etaMax;
      edm::Service<TFileService>  fs;
      TH2D *h2_BTaggingEff_Denom_b;
      TH2D *h2_BTaggingEff_Denom_c;
      TH2D *h2_BTaggingEff_Denom_udsg;
      TH2D *h2_BTaggingEff_Num_b;
      TH2D *h2_BTaggingEff_Num_c;
      TH2D *h2_BTaggingEff_Num_udsg;

      edm::EDGetTokenT<PatJetCollection> JetToken_;

};

BTaggingEffAnalyzer::BTaggingEffAnalyzer(const edm::ParameterSet& iConfig) :

  jetsTag(iConfig.getParameter<edm::InputTag>("JetsTag")),
  discriminatorTag(iConfig.getParameter<std::string>("DiscriminatorTag")),
  discriminatorValue(iConfig.getParameter<double>("DiscriminatorValue")),
  ptNBins(iConfig.getParameter<int>("PtNBins")),
  ptMin(iConfig.getParameter<double>("PtMin")),
  ptMax(iConfig.getParameter<double>("PtMax")),
  etaNBins(iConfig.getParameter<int>("EtaNBins")),
  etaMin(iConfig.getParameter<double>("EtaMin")),
  etaMax(iConfig.getParameter<double>("EtaMax"))

{

   JetToken_ = consumes<pat::JetCollection>(jetsTag);

   //now do what ever initialization is needed
   h2_BTaggingEff_Denom_b    = fs->make<TH2D>("h2_BTaggingEff_Denom_b", ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax);
   h2_BTaggingEff_Denom_c    = fs->make<TH2D>("h2_BTaggingEff_Denom_c", ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax);
   h2_BTaggingEff_Denom_udsg = fs->make<TH2D>("h2_BTaggingEff_Denom_udsg", ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax);
   h2_BTaggingEff_Num_b    = fs->make<TH2D>("h2_BTaggingEff_Num_b", ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax);
   h2_BTaggingEff_Num_c    = fs->make<TH2D>("h2_BTaggingEff_Num_c", ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax);
   h2_BTaggingEff_Num_udsg = fs->make<TH2D>("h2_BTaggingEff_Num_udsg", ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax);
}


BTaggingEffAnalyzer::~BTaggingEffAnalyzer()
{
}

void BTaggingEffAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<PatJetCollection> jets;
  iEvent.getByToken(JetToken_,jets);

  // loop over jets
  for(PatJetCollection::const_iterator it = jets->begin(); it != jets->end(); ++it)
  {
    int partonFlavor = it->partonFlavour();
    if( abs(partonFlavor)==5 )
    {
      h2_BTaggingEff_Denom_b->Fill(it->pt(), it->eta());
      if( it->bDiscriminator(discriminatorTag.c_str()) >= discriminatorValue ) h2_BTaggingEff_Num_b->Fill(it->pt(), it->eta());
    }
    else if( abs(partonFlavor)==4 )
    {
      h2_BTaggingEff_Denom_c->Fill(it->pt(), it->eta());
      if( it->bDiscriminator(discriminatorTag.c_str()) >= discriminatorValue ) h2_BTaggingEff_Num_c->Fill(it->pt(), it->eta());
    }
    else
    {
      h2_BTaggingEff_Denom_udsg->Fill(it->pt(), it->eta());
      if( it->bDiscriminator(discriminatorTag.c_str()) >= discriminatorValue ) h2_BTaggingEff_Num_udsg->Fill(it->pt(), it->eta());
    }
  }
}


void BTaggingEffAnalyzer::beginJob() {}

void BTaggingEffAnalyzer::endJob() {}

void BTaggingEffAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {}

void BTaggingEffAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {}

void BTaggingEffAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

void BTaggingEffAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

void BTaggingEffAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BTaggingEffAnalyzer);
