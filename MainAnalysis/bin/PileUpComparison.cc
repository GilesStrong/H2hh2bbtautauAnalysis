#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <strstream>

#include <TH1F.h>
#include <TTree.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"

#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"

#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoProduct.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"

#include "DataFormats/Common/interface/MergeableCounter.h"


int main(int argc, char* argv[])
{
  using pat::Muon;

  gSystem->Load( "libFWCoreFWLite" );
  FWLiteEnabler::enable();

  // command line parser:
  optutl::CommandLineParser parser ("Analyze FWLite Histograms");
  parser.integerValue ("maxEvents"  ) =    -1;
  parser.integerValue ("outputEvery") = 10000;
  parser.stringValue  ("outputFile" ) = "analyzeFWLiteHistograms.root";
  parser.addOption ("runOnData",     optutl::CommandLineParser::kBool, "", false);
  parser.addOption ("debug",         optutl::CommandLineParser::kBool, "", false);

  parser.parseArguments (argc, argv);
  int maxEvents_ = parser.integerValue("maxEvents");
  unsigned int outputEvery_ = parser.integerValue("outputEvery");
  std::string outputFile_ = parser.stringValue("outputFile");
  std::vector<std::string> inputFiles_ = parser.stringVector("inputFiles");
  bool runOnData = parser.boolValue("runOnData");
  bool debug = parser.boolValue("debug");

  // book a set of histograms:
  fwlite::TFileService fs = fwlite::TFileService(outputFile_.c_str());
  TFileDirectory dir = fs.mkdir("noPU");
  TH1F* muonPt_noPU_  = dir.make<TH1F>("muonPt"  , "pt"  ,   100,   0., 300.);
  TH1F* muonEta_noPU_ = dir.make<TH1F>("muonEta" , "eta" ,   100,  -3.,   3.);
  TH1F* muonPhi_noPU_ = dir.make<TH1F>("muonPhi" , "phi" ,   100,  -5.,   5.);
  TH1F* nVert_noPU_   = dir.make<TH1F>("nVert" ,   "nVert" ,  50,    0,   50);
  TFileDirectory dir2 = fs.mkdir("weighted");
  TH1F* muonPt_  = dir2.make<TH1F>("muonPt"  , "pt"  ,   100,   0., 300.);
  TH1F* muonEta_ = dir2.make<TH1F>("muonEta" , "eta" ,   100,  -3.,   3.);
  TH1F* muonPhi_ = dir2.make<TH1F>("muonPhi" , "phi" ,   100,  -5.,   5.); 
  TH1F* nVert_   = dir2.make<TH1F>("nVert" ,   "nVert" ,  50,    0,   50);   
  TH1F* m_kinfit_  = dir2.make<TH1F>("m_kinfit" ,   "KinFit Mass" ,  100,    0,   500); 

  double xsec = -1;
  int nev_total = 0;
  int nev_selected = 0;
  double weight = 1;
  double weight_noPU = 1;
  double weight_lumi = 1;
  double weight_PU = 1;
  double weight_btag = 1;
  double weight_muon_IdScaleFactor = 1;
  double weight_muon_triggerScaleFactor = 1;

  double m_kinfit = 0;

  // loop over all input files:
  int ievt=0;
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
    std::string filename = inputFiles_[iFile];
    TFile* inFile = TFile::Open(filename.c_str());

    if(filename.find("Data") != std::string::npos || filename.find("Run201") != std::string::npos ) runOnData = true;

    if( inFile ){
      
      // get luminosity weights:
      if(!runOnData){
	// get cross section from GenInfo:
	TTree *runs = (TTree*)inFile->Get("Runs");
	runs->GetEntry();
	GenRunInfoProduct myRunInfoProduct;
	TBranch* genInfoBranch = runs->GetBranch("GenRunInfoProduct_generator__SIM.obj");
          
	if(genInfoBranch){
	  genInfoBranch->SetAddress(&myRunInfoProduct);
	  genInfoBranch->GetEntry(0);
	  xsec = (double)myRunInfoProduct.crossSection();
	}
          
	// loop over all luminosity blocks and get number of events:
	fwlite::LuminosityBlock lumiBlock(inFile);
	int ilumiBlock = 0;
	for(lumiBlock.toBegin(); !lumiBlock.atEnd(); ++lumiBlock, ++ilumiBlock){
	  edm::Handle<edm::MergeableCounter> hCnt;
	  lumiBlock.getByLabel(std::string("nEventsTotal"), hCnt);
	  nev_total = nev_total + hCnt->value;
            
	  edm::Handle<edm::MergeableCounter> hCntSel;
	  lumiBlock.getByLabel(std::string("nEventsSelected"), hCntSel);
	  nev_selected = nev_selected + hCntSel->value;
	}
	if(xsec>0) weight_lumi = xsec / nev_total;

	if(debug){
	  std::cout << " xsec = " << xsec << "\n";
	  std::cout << " nev_total = " << nev_total << "\n";
	  std::cout << " nev_selected = " << nev_selected << "\n";
	  std::cout << " weight_lumi = " << weight_lumi << "\n";
	}
      }

      // loop over all events:
      fwlite::Event ev(inFile);
      for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
	edm::EventBase const & event = ev;
        if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
        if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false)
	  std::cout << "  processing event: " << ievt << std::endl;

        // get all handles present in both data and MC:
	edm::Handle<std::vector<Muon>> selectedmuons;
        event.getByLabel(edm::InputTag("selectedMuons"), selectedmuons);
	edm::Handle<std::vector<pat::Tau>> selectedtaus;
        event.getByLabel(edm::InputTag("selectedTaus"), selectedtaus);
	edm::Handle<std::vector<pat::Jet>> selectedjets;
        event.getByLabel(edm::InputTag("selectedJets"), selectedjets);
	edm::Handle<std::vector<pat::MET>> slimmedMET;
        event.getByLabel(edm::InputTag("slimmedMETs","","SKIM"), slimmedMET);
	edm::Handle<edm::ValueMap<float>> offlineSlimmedPrimaryVertices;
        event.getByLabel(edm::InputTag("offlineSlimmedPrimaryVertices"), offlineSlimmedPrimaryVertices);

	edm::Handle<std::vector<double>> bTaggingSF;
        event.getByLabel(edm::InputTag("bTaggingSF"), bTaggingSF);
        
        // edm::Handle<double> muonIDScaleFactorVal;
        // event.getByLabel(edm::InputTag("muonIDScaleFactor", "ScaleFactorValue","SKIM2"), muonIDScaleFactorVal);
        
	edm::Handle<double> METSignificance;
        event.getByLabel(edm::InputTag("METSignificance"), METSignificance);

	edm::Handle<double> kinfit_p;
        event.getByLabel(edm::InputTag("kinfit", "P"), kinfit_p);
	edm::Handle<double> kinfit_chi2;
        event.getByLabel(edm::InputTag("kinfit", "chi2"), kinfit_chi2);
	edm::Handle<double> kinfit_mH;
        event.getByLabel(edm::InputTag("kinfit", "mH"), kinfit_mH);
	edm::Handle<double> kinfit_convergence;
        event.getByLabel(edm::InputTag("kinfit", "convergence"), kinfit_convergence);
        
        // get specific handles for MC:
        if(!runOnData){
	  edm::Handle<double> PUweight;
	  event.getByLabel(edm::InputTag("PUWeightProducer","PUWeight","SKIM"), PUweight);
	  weight_PU = *PUweight.product();
            
	  // get muon weights:
	  edm::Handle<double> muon_triggerSF;
	  event.getByLabel(edm::InputTag("muonTriggerScaleFactor", "ScaleFactorValue","SKIM2"), muon_triggerSF);
	  weight_muon_triggerScaleFactor = *muon_triggerSF.product();
	  edm::Handle<double> muon_IDSF;
	  event.getByLabel(edm::InputTag("muonIDScaleFactor", "ScaleFactorValue","SKIM2"), muon_IDSF);
	  weight_muon_IdScaleFactor = *muon_IDSF.product();    

	  if(debug) std::cout << " weight_PU = " << weight_PU << std::endl;
        }
        
        if(debug){ std::cout << " muon collection size: " << selectedmuons->size() << std::endl;
	  std::cout << " btagging weight vector size: " << bTaggingSF->size() << std::endl;
	}

        // loop over btag weights:
        for(std::vector<double>::const_iterator it = bTaggingSF->begin(); it != bTaggingSF->end(); ++it) {
	  weight_btag = weight_btag * *it;
        }
        
        if(debug) std::cout << " b-tagging weight = " << weight_btag << std::endl;
        if(debug) std::cout << " weight_muon_IdScaleFactor * weight_btag = " << weight_muon_IdScaleFactor * weight_btag << std::endl;
      
        // apply all weights:
        weight_noPU = weight_lumi * weight_muon_IdScaleFactor * weight_muon_triggerScaleFactor; //* weight_btag;
        weight =      weight_lumi * weight_PU * weight_muon_IdScaleFactor * weight_muon_triggerScaleFactor; //* weight_btag;
        if(debug) std::cout << " weight = " << weight << std::endl;
        
        // loop muon collection and fill histograms:
        for(std::vector<Muon>::const_iterator mu1=selectedmuons->begin(); mu1!=selectedmuons->end(); ++mu1){
	  muonPt_ ->Fill( mu1->pt (), weight );
	  muonEta_->Fill( mu1->eta(), weight );
	  muonPhi_->Fill( mu1->phi(), weight );
	  muonPt_noPU_ ->Fill( mu1->pt (), weight_noPU );
	  muonEta_noPU_->Fill( mu1->eta(), weight_noPU );
	  muonPhi_noPU_->Fill( mu1->phi(), weight_noPU );
        }

	//plot Kinfit mass
	m_kinfit = *kinfit_mH.product();
	m_kinfit_ -> Fill(m_kinfit, weight);

        //edm::ValueMap<float> vtx = *offlineSlimmedPrimaryVertices;
        //double nVert = vtx.size();
        double nVert = offlineSlimmedPrimaryVertices->size();

        // number of vertices:
        nVert_->Fill( nVert, weight );
        nVert_noPU_->Fill( nVert, weight_noPU );
      }
    
      inFile->Close();
    }
    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
  }
  return 0;
}
