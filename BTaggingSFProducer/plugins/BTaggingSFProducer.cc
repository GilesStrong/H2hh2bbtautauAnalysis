// calculation of b-tagging SF using officially provided SF file in CSV format
// SF dependent on flavour, p_T and eta

#include <memory>
#include <vector>
#include <TLorentzVector.h>
#include <math.h>
#include <exception>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

typedef std::vector<pat::Jet> PatJetCollection;

class BTaggingSFProducer : public edm::stream::EDProducer<> {
			public:
								explicit BTaggingSFProducer(const edm::ParameterSet&);
								~BTaggingSFProducer();

								static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

			private:
								virtual void beginStream(edm::StreamID) override;
								virtual void produce(edm::Event&, const edm::EventSetup&) override;
								virtual void endStream() override;

								std::string CSVFile_;
								std::string EffMapFile_;
								std::string EffMapFolder_;

								TFile *f_EffMap_AK4 = 0;

								bool debug_;
								bool runOnData_;

								edm::Service<TFileService>  fs;
								edm::EDGetTokenT<PatJetCollection> JetToken_;

								std::string discriminatorTag_;
								double discriminatorValue_;

								TH2D *h2_EffMapB_AK4;
								TH2D *h2_EffMapC_AK4;
								TH2D *h2_EffMapUDSG_AK4;

};

BTaggingSFProducer::BTaggingSFProducer(const edm::ParameterSet& iConfig)
{
								produces<double>("BTaggingSF").setBranchAlias("BTaggingSF");
								JetToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("JetsTag"));        
								runOnData_ = iConfig.getParameter<bool>("runOnData");
								debug_ = iConfig.getParameter<bool>("debug");
								discriminatorTag_ = iConfig.getParameter<std::string>("DiscriminatorTag");
								discriminatorValue_ = iConfig.getParameter<double>("DiscriminatorValue");
								CSVFile_ = iConfig.getParameter<std::string>("CSVFile");
								EffMapFile_ = iConfig.getParameter<std::string>("EffMapFile");
								EffMapFolder_ = iConfig.getParameter<std::string>("EffMapFolder");
}


BTaggingSFProducer::~BTaggingSFProducer(){}


void BTaggingSFProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
			using namespace edm;
			
			if(debug_) std::cout << "starting BTSFP " << std::endl;
				
			double weight;
		

			Handle<PatJetCollection> jets;
			iEvent.getByToken(JetToken_,jets);
	
			// read BTagging SF from CSV file:
			BTagCalibration calib_csv("csv",CSVFile_);
			
			// we need two CSV reader instances, one for light and one for heavy jets
			BTagCalibrationReader readerHeavy(BTagEntry::OP_MEDIUM,          //operating point
																								"central");                                //systematics type
			readerHeavy.load(&calib_csv, //calibration instance
																									BTagEntry::FLAV_B, //Jet Flavour (Check this GS: 23-06-17)
																									"comb");              //measurement type                   
			BTagCalibrationReader readerLight(BTagEntry::OP_MEDIUM,                      //operating point
																									"central");                                //systematics type
			readerLight.load(&calib_csv, //calibration instance
																									BTagEntry::FLAV_UDSG,
																									"incl");              //measurement type 
			if(debug_) std::cout << "reader set up " << std::endl;
			if(debug_) std::cout << "Eff map file " << EffMapFile_.c_str() << std::endl;

			// import efficiency map:
			f_EffMap_AK4 = TFile::Open(EffMapFile_.c_str());
			
			std::string filename = EffMapFolder_ + "/efficiency_b";
			h2_EffMapB_AK4 = (TH2D*)f_EffMap_AK4->Get(filename.c_str());
			filename = EffMapFolder_ + "/efficiency_c";
			h2_EffMapC_AK4 = (TH2D*)f_EffMap_AK4->Get(filename.c_str());
			filename = EffMapFolder_ + "/efficiency_udsg";
			h2_EffMapUDSG_AK4 = (TH2D*)f_EffMap_AK4->Get(filename.c_str());

			double pMC = 1.0;
			double pData = 1.0;
			
			if(debug_) std::cout << "imported efficiency map " << std::endl;
			
			// loop over jets
			for(PatJetCollection::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {
							
								// do the following only for MC:
								if (!runOnData_){
																
												bool bTagged = false;

												double csv = jet->bDiscriminator(discriminatorTag_);
												if( jet->bDiscriminator(discriminatorTag_) >= discriminatorValue_ ) {
																bTagged = true;
												}
												
												if(debug_) std::cout << "btagged: " << bTagged << std::endl;

												int partonFlavor = abs(jet->partonFlavour());
												double eta = fabs(jet->eta());
												double pt = fabs(jet->pt());
												if( eta>2.4) continue;
												if( partonFlavor==0) continue; //for jets with flavor 0, we ignore.
												if( csv >= 1 || csv < 0 ) continue;

												// get efficiency from efficiency map:
												double eff = 1;
												if( partonFlavor==5 ) {
														///the pt/eta dependent efficiency for b-tag for "b jet"
														eff = h2_EffMapB_AK4->GetBinContent( h2_EffMapB_AK4->GetXaxis()->FindBin(pt), h2_EffMapB_AK4->GetYaxis()->FindBin(fabs(eta)) );
												}else if( partonFlavor==4){
														///the pt/eta dependent efficiency for b-tag for "c jet"
														eff = h2_EffMapC_AK4->GetBinContent( h2_EffMapC_AK4->GetXaxis()->FindBin(pt), h2_EffMapC_AK4->GetYaxis()->FindBin(fabs(eta)) );
												}else{
														///the pt/eta dependent efficiency for b-tag for "light jet"
														eff = h2_EffMapUDSG_AK4->GetBinContent( h2_EffMapUDSG_AK4->GetXaxis()->FindBin(pt), h2_EffMapUDSG_AK4->GetYaxis()->FindBin(fabs(eta)) );
												}
												if(debug_) std::cout << "b-tagging efficiency: " << eff << std::endl;


												// get pT/eta-dependent scale factor from official CSV file:
												double SF = 1.0;
												if (pt>800) pt = 800;
												if (pt<20) pt=20;
				
												if ( partonFlavor == 5 ) {
														SF = readerHeavy.eval(BTagEntry::FLAV_B, eta,  pt);
												}
												else if ( partonFlavor == 4 ) {
														SF = readerHeavy.eval(BTagEntry::FLAV_C, eta, pt);
															}
												else {
														SF = readerLight.eval(BTagEntry::FLAV_UDSG, eta, pt);
												}

												if(debug_) std::cout << "b-tagging SF: " << SF << std::endl;            
												
												if (bTagged) {
																pMC = pMC * eff;
																pData = pData * eff * SF;
												} else {
																pMC = pMC * (1 - eff);
																pData = pData * (1 - eff*SF);
												}
								}
			}
			
			if(pData != 0) {
							weight = pMC/pData;
			} else {
							weight = 0;
			}
			
			if(debug_) std::cout << "b-tagging final weight: " << weight << std::endl;     
			
			
			// save btagging SF in event:
			std::auto_ptr<double> output(new double(weight));
			iEvent.put(output, "BTaggingSF");
			
}

void BTaggingSFProducer::beginStream(edm::StreamID){}

void BTaggingSFProducer::endStream(){}

	
void BTaggingSFProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
		edm::ParameterSetDescription desc;
		desc.setUnknown();
		descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BTaggingSFProducer);
