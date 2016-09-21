#!/bin/bash

startingpath="$PWD"

cd $CMSSW_BASE/src/

pfad="$PWD"
echo $pfad

if [ -d "$pfad/HTT-utilities" ]
then
    echo "The directory HTT-utilities already exists!!"
else 
    echo "Getting Lepton Efficiency (interface)..." 
    git clone https://github.com/CMS-HTT/LeptonEff-interface $CMSSW_BASE/src/HTT-utilities/LepEffInterface
    echo "Getting Lepton Efficiency (data)..."
    git clone https://github.com/CMS-HTT/LeptonEfficiencies.git $CMSSW_BASE/src/HTT-utilities/LepEffInterface/data
fi

echo "You are now in:" 
cd $startingpath
echo $startingpath