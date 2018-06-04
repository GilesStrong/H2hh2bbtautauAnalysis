#!/bin/bash

# get pile-up distributons for data / MC for reweighting

# PU in 2016 data:
# "true" and "observed" PU distributions: 
# see https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData#Pileup_JSON_Files_For_Run_II
# 
# "true" mode puts mean number of interactions per crossing into the histogram, "observed" uses properly-normalized poisson distribution with a mean corresponding to the expected number of interactions per crossing for each LumiSection.

pileupCalc.py -i Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt --inputLumiJSON pileup_latest.txt  --calcMode true --minBiasXsec 63000 --maxPileupBin 50 --numPileupBins 50  MyDataPileupHistogram_true.root

pileupCalc.py -i Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt --inputLumiJSON pileup_latest.txt  --calcMode observed --minBiasXsec 63000 --maxPileupBin 50 --numPileupBins 50  MyDataPileupHistogram_observed.root

# PU used in Spring16 MC campaign:
#
# see https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupInformation#Pileup_Information_in_MC_full_sim
# conditions available in /SimGeneral/MixingModule/python/mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU_cfi.py

root -q -b createPileUpMC2016.cxx
