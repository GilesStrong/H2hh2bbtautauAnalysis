#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"


int main(int argc, char* argv[]) 
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  FWLiteEnabler::enable();

  // initialize command line parser
  optutl::CommandLineParser parser ("Analyze ScaleFactors with FWLite");

  // set defaults
  parser.integerValue ("maxEvents"  ) = 1000;

  // parse arguments
  parser.parseArguments (argc, argv);
  int maxEvents_ = parser.integerValue("maxEvents");
  std::vector<std::string> inputFiles_ = parser.stringVector("inputFiles");

  // loop the events
  int ievt=0;  
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
    // open input file (can be located on castor)
    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    if( inFile ){
      fwlite::Event ev(inFile);
      for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
	edm::EventBase const & event = ev;
	// break loop if maximal number of events is reached 
	if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;

	// Handle to the muon collection
	edm::Handle<std::vector<pat::Muon> > cand;
	edm::Handle<double> muon_eff_val;
        edm::Handle<double> muon_eff_err;
	edm::Handle<double> muon_sf_val;
        edm::Handle<double> muon_sf_err;

	event.getByLabel(edm::InputTag("selectedMuons"), cand);
	event.getByLabel(edm::InputTag("muonScaleFactor","EfficiencyValue"), muon_eff_val);
        //event.getByLabel(edm::InputTag("muonScaleFactor","EfficiencyError"), muon_eff_err);
        event.getByLabel(edm::InputTag("muonScaleFactor","ScaleFactorValue"), muon_sf_val);
        //event.getByLabel(edm::InputTag("muonScaleFactor","ScaleFactorError"), muon_sf_err);

	// loop muon collection and fill histograms
        int mucount=0;
	for(std::vector<pat::Muon>::const_iterator mu = cand->begin(); mu != cand->end(); ++mu){
          std::cout << "muon " << ++mucount << "  :"
                    << "  pt=" << mu->pt() 
                    << ", eta= " << mu->eta() 
                    << std::endl;
        }
        std::cout << "  eff=" << *muon_eff_val //<< "+-" << *muon_eff_err
                  << ", sf=" << *muon_sf_val //<< "+-" << *muon_sf_err
                  << std::endl;
        std::cout << "-------------------------------" << std::endl;

      }  
      // close input file
      inFile->Close();
    }
    // break loop if maximal number of events is reached:
    // this has to be done twice to stop the file loop as well
    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
  }
  return 0;
}
