#include <memory>
#include <vector>
#include <math.h>
#include <exception>

#include "TROOT.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

class PUWeightProducer : public edm::stream::EDProducer<> {
   public:
      explicit PUWeightProducer(const edm::ParameterSet&);
      ~PUWeightProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      std::string generatedFile_;
      std::string dataFile_;
      std::string GenHistName_;
      std::string DataHistName_;
      std::string pileupSummaryInfoLabel_;
      
      double PUweight_;
      bool debug_;
      bool runOnData_;
      
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PUsummaryToken_;

};

PUWeightProducer::PUWeightProducer(const edm::ParameterSet& iConfig)
{
      produces<double>("PUWeight").setBranchAlias("PUWeight");
      PUsummaryToken_ = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<std::string>("pileupSummaryInfoLabel"));

      generatedFile_ = iConfig.getParameter<std::string>("generatedFile");
      dataFile_ = iConfig.getParameter<std::string>("dataFile");
      GenHistName_ = iConfig.getParameter<std::string>("GenHistName");
      DataHistName_ = iConfig.getParameter<std::string>("GenHistName");
      runOnData_ = iConfig.getParameter<bool>("runOnData");
      debug_ = iConfig.getParameter<bool>("debug");
}


PUWeightProducer::~PUWeightProducer(){}


void PUWeightProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // check if not running on data:
   if (runOnData_) return;
   
   edm::LumiReWeighting *LumiWeights;
   
   if (generatedFile_ == "default" || dataFile_ == "default") {
       
        // when using "default", use hardcoded PU reweighting histograms:

        std::vector<float> pileupMC(50);
        pileupMC =
              {0.0,
              0.000829312873542,
              0.00124276120498,
              0.00339329181587,
              0.00408224735376,
              0.00383036590008,
              0.00659159288946,
              0.00816022734493,
              0.00943640833116,
              0.0137777376066,
              0.017059392038,
              0.0213193035468,
              0.0247343174676,
              0.0280848773878,
              0.0323308476564,
              0.0370394341409,
              0.0456917721191,
              0.0558762890594,
              0.0576956187107,
              0.0625325287017,
              0.0591603758776,
              0.0656650815128,
              0.0678329011676,
              0.0625142146389,
              0.0548068448797,
              0.0503893295063,
              0.040209818868,
              0.0374446988111,
              0.0299661572042,
              0.0272024759921,
              0.0219328403791,
              0.0179586571619,
              0.0142926728247,
              0.00839941654725,
              0.00522366397213,
              0.00224457976761,
              0.000779274977993,
              0.000197066585944,
              7.16031761328e-05,
              0.0,
              0.0,
              0.0,
              0.0,
              0.0,
              0.0,
              0.0,
              0.0,
              0.0,
              0.0,
              0.0};
    
       std::vector<float> DataTrue(50);
       DataTrue =
              {0.000000,
               5046.461095,
               240725.407653,
               783167.975339,
               1735164.036137,
               2367877.136296,
               3413712.978956,
               6120882.509950,
               24306941.636897,
               67829843.914266,
               144549067.470500,
               256920079.998321,
               406072497.963521,
               562521112.978453,
               705656855.020352,
               840607940.954170,
               953662997.748009,
               1029757706.806727,
               1064909937.082736,
               1060910036.865547,
               1019587527.735262,
               946581187.455070,
               851484805.248377,
               741086543.942303,
               619352639.058317,
               492896076.542409,
               372173111.134789,
               266851915.290031,
               181846804.443556,
               117532951.416225,
               71775609.196965,
               41314474.310614,
               22420762.026694,
               11488854.385398,
               5568219.230968,
               2558909.680524,
               1119948.472484,
               470202.650510,
               191720.109574,
               77771.504554,
               33019.033352,
               16074.953335,
               9872.036962,
               7673.367502,
               6917.364819,
               6659.508026,
               6558.478631,
               6487.194493,
               6398.026661,
               6277.763350,
              };
       
       LumiWeights = new edm::LumiReWeighting(pileupMC, DataTrue);
      
   } else {
       LumiWeights = new edm::LumiReWeighting(generatedFile_, dataFile_, GenHistName_, DataHistName_);
   }
   
   // now use LumiReWeighting to get the event weight according to its pile-up content:
   // also see https://twiki.cern.ch/twiki/bin/view/CMS/PileupMCReweightingUtilities

   Handle<std::vector<PileupSummaryInfo> > pileUpInfo;
   iEvent.getByToken(PUsummaryToken_,pileUpInfo);
   std::vector<PileupSummaryInfo>::const_iterator PVI;
   
   float Tnpv = -1;
   int BX = -1;
   for(PVI = pileUpInfo->begin(); PVI != pileUpInfo->end(); ++PVI) {
       BX = PVI->getBunchCrossing();
   
       if(BX == 0) { 
           Tnpv = PVI->getTrueNumInteractions();
           continue;
       }
   }
   
   PUweight_ = LumiWeights->weight( Tnpv );

   if (debug_) {
       std::cout << "-----------------------------------------------------------\n";
       std::cout << " true number of interactions (mean value): " << Tnpv << "\n";
       std::cout << " bunch crossing BX:                        " << BX << "\n";
       std::cout << " calculated PU weight:                     " << PUweight_ << "\n";
       std::cout << "-----------------------------------------------------------\n";
   }
  
   // save PU weight in event:
   std::auto_ptr<double> output(new double(PUweight_));
   iEvent.put(output, "PUWeight");

}

void PUWeightProducer::beginStream(edm::StreamID){}

void PUWeightProducer::endStream(){}

 
void PUWeightProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PUWeightProducer);
