#include <string>
#include <vector>
#include <iostream>

#include <TSystem.h>
#include <TFile.h>

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

int main(int argc, char* argv[])
{
  gSystem->Load( "libFWCoreFWLite" );
  FWLiteEnabler::enable();

  // command line parser:
  optutl::CommandLineParser parser ("Analyze FWLite Histograms");
  parser.parseArguments (argc, argv);
  std::vector<std::string> inputFiles_ = parser.stringVector("inputFiles");

  double xsec = -1;
  int nev_total = 0;
  int nev_selected = 0;
  double weight_lumi = 0;

  // loop over all input files:
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
    std::string filename = inputFiles_[iFile];
    TFile* inFile = TFile::Open(filename.c_str());

    if( inFile ){
      // loop over all runs and get cross section:
      int iRun=0;
      fwlite::Run run(inFile);
      for(run.toBegin(); !run.atEnd(); ++run, ++iRun){
        edm::Handle<GenRunInfoProduct> hGrip;
        run.getByLabel(std::string("generator"), hGrip);
        xsec=(double)hGrip->crossSection();
      }
      
      // loop over all luminosity blocks and get number of events:
      fwlite::LuminosityBlock lumiBlock(inFile);
      int ilumiBlock = 0;
      for(lumiBlock.toBegin(); !lumiBlock.atEnd(); ++lumiBlock, ++ilumiBlock){
        edm::Handle<edm::MergeableCounter> hCnt;
        lumiBlock.getByLabel(std::string("nEventsTotal"), hCnt);
        nev_total = nev_total + hCnt->value;
        
        edm::Handle<reco::VertexCollection> vertices;
        lumiBlock.getByLabel(std::string("nEventsSelected"), hCntSel);
        nev_selected = nev_selected + hCntSel->value;
      }
      
      if(xsec>0) weight_lumi = xsec / nev_total;
      
      std::cout << " xsec = " << xsec << "\n";
      std::cout << " nev_total = " << nev_total << "\n";
      std::cout << " nev_selected = " << nev_selected << "\n";
      std::cout << " weight_lumi = " << weight_lumi << "\n";
      std::cout << " accessed runs:        " << iRun+1 << "\n";
      std::cout << " accessed lumi blocks: " << ilumiBlock+1 << "\n";
    }
  }
  return 0;
}