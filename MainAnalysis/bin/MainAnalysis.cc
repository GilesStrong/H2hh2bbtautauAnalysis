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
#include <TLorentzVector.h>
#include <math.h>

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
  parser.addOption ("setdebug",      optutl::CommandLineParser::kBool, "", false);

  parser.parseArguments (argc, argv);
  int maxEvents_ = parser.integerValue("maxEvents");
  unsigned int outputEvery_ = parser.integerValue("outputEvery");
  std::string outputFile_ = parser.stringValue("outputFile");
  std::vector<std::string> inputFiles_ = parser.stringVector("inputFiles");
  bool runOnData = parser.boolValue("runOnData");
  bool debug = parser.boolValue("setdebug");

  // book a set of histograms:
  fwlite::TFileService fs = fwlite::TFileService(outputFile_.c_str());
  
  TFileDirectory dir_noPU = fs.mkdir("noPU");
  TFileDirectory dir_weighted = fs.mkdir("weighted");
  TH1F* hmuonPt      = dir_weighted.make<TH1F>("muonPt"  , "pt"  ,   100,   0., 300.);
  TH1F* h_nPrimaryVtx = dir_weighted.make<TH1F>("nPrimaryVtx"  , "Number of primary vertices"  ,   100,   0., 100.);
  TH1F* h_nPrimaryVtx_unweighted = dir_weighted.make<TH1F>("nPrimaryVtx_unweighted"  , "Number of primary vertices"  ,   100,   0., 100.);
  // ...

  double weight = 1;
  double weight_lumi = 1;
  double xsec = -1;
  int nev_total = 0;
  int nev_selected = 0;

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
        edm::Handle<vector<reco::Vertex> > vertices;
        ev.getByLabel(std::string("offlineSlimmedPrimaryVertices"), vertices);
        
        if(!runOnData){
            edm::Handle<double> PUweight;
            event.getByLabel(edm::InputTag("PUWeightProducer", "PUWeight"), PUweight);
            weight = *PUweight.product();
            std::cout << " weight_PU = " << weight << std::endl;
            weight *= weight_lumi;
        }
 
        h_nPrimaryVtx->Fill(vertices.size(), weight);
        h_nPrimaryVtx_unweighted->Fill(vertices.size(), weight_lumi);
      }
    
      inFile->Close();
    }
    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
  }
  return 0;
}
