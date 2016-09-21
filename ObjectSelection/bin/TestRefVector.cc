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

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"


int main(int argc, char* argv[]) 
{
  // define what muon you are using; this is necessary as FWLite is not 
  // capable of reading edm::Views
  using pat::Muon;

  // ----------------------------------------------------------------------
  // First Part: 
  //
  //  * enable FWLite 
  //  * book the histograms of interest 
  //  * open the input file
  // ----------------------------------------------------------------------

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  FWLiteEnabler::enable();

  // initialize command line parser
  optutl::CommandLineParser parser ("Analyze FWLite Histograms");

  // set defaults
  parser.integerValue ("maxEvents"  ) = 1000;
  parser.integerValue ("outputEvery") =   10;
  parser.stringValue  ("outputFile" ) = "analyzeFWLiteHistograms.root";

  // parse arguments
  parser.parseArguments (argc, argv);
  int maxEvents_ = parser.integerValue("maxEvents");
  unsigned int outputEvery_ = parser.integerValue("outputEvery");
  std::string outputFile_ = parser.stringValue("outputFile");
  std::vector<std::string> inputFiles_ = parser.stringVector("inputFiles");

  // loop the events
  int ievt=0;  
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
    // open input file (can be located on castor)
    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    if( inFile ){
      // ----------------------------------------------------------------------
      // Second Part: 
      //
      //  * loop the events in the input file 
      //  * receive the collections of interest via fwlite::Handle
      //  * fill the histograms
      //  * after the loop close the input file
      // ----------------------------------------------------------------------      
      fwlite::Event ev(inFile);
      for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
	edm::EventBase const & event = ev;
	// break loop if maximal number of events is reached 
	if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
	// simple event counter
	if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false) 
	  std::cout << "  processing event: " << ievt << std::endl;

	// Handle to the muon collection
	edm::Handle<std::vector<Muon> > muons;
	event.getByLabel(std::string("slimmedMuons"), muons);
        edm::Handle< edm::RefVector<std::vector<Muon> > > selectedmuons1;
        event.getByLabel(std::string("testMuonsPredef"), selectedmuons1);
        edm::Handle< edm::RefVector<std::vector<Muon> > > selectedmuons2;
        event.getByLabel(std::string("testMuonsOwn"), selectedmuons2);

        std::cout << "######################################" << std::endl;	
        std::cout << "size of original muon collection:           " << muons->size() << std::endl;
        std::cout << "size of selected (predef)  muon collection: " << selectedmuons1->size() << std::endl;
        std::cout << "size of selected (own)  muon collection:    " << selectedmuons2->size() << std::endl;
        std::cout << "--------------------------------------" << std::endl;

	// loop muon collection and fill histograms
        for(unsigned int i = 0; i < muons->size(); i++)
          std::cout << "original collection: pt=" << muons->at(i).pt() << std::endl;
	for(unsigned int i = 0; i < selectedmuons1->size(); i++)
          std::cout << "selected (predef) collection: pt=" << selectedmuons1->at(i)->pt() << std::endl;
        for(unsigned int i = 0; i < selectedmuons2->size(); i++)
          std::cout << "selected (own) collection:    pt=" << selectedmuons2->at(i)->pt() << std::endl;
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
