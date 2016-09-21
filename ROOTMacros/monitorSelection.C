#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"

/// Fill the cutflow plot from 'yield' histograms
TH1F* cutflow(TFile* sample, const std::vector<std::string>& validSteps);

void monitorSelection(const std::string& histName="yield")
{
  gStyle->SetOptStat(0);

  // list of valid selection steps
  std::vector<std::string> validSteps_;
  validSteps_.push_back("1"   );
  validSteps_.push_back("2"   );
  validSteps_.push_back("3"  );
  validSteps_.push_back("4"   );
  validSteps_.push_back("5"   );
  validSteps_.push_back("6"   );
  
  //open files
  TFile* file = new TFile("../Configurations/analyzeSelection.root");

  // load histograms
  TH1F* hist = cutflow(file, validSteps_);

  // setup the canvas and draw the histograms
  TCanvas* canv0 = new TCanvas("canv0", "canv0", 600, 600);
  canv0->cd(0);
  canv0->SetLogy(1);
  hist->SetTitle("Selection Steps");
  hist->SetMinimum(1.);
  hist->SetFillColor(kYellow);
  hist->Draw();
}

TH1F* cutflow(TFile* sample, const std::vector<std::string>& validSteps)
{
  // book histogram
  TH1F* hist = new TH1F(sample->GetName(), sample->GetName(), validSteps.size(), 0., validSteps.size());
  // set labels
  for(unsigned int idx=0; idx<validSteps.size(); ++idx){
    hist->GetXaxis()->SetBinLabel( idx+1 , validSteps[idx].c_str());
  }
  hist->LabelsOption("h", "X"); //"h", "v", "u", "d"
  // fill histogram
  for(unsigned int idx=0; idx<validSteps.size(); ++idx){
    TH1F* buffer = (TH1F*)sample->Get((std::string("mon").append(validSteps[idx]).append("/yield")).c_str());
    hist->SetBinContent(idx+1, buffer->GetBinContent(1));
  }
  return hist;
}
