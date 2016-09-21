#!/bin/bash

startingpath="$PWD"

cd $CMSSW_BASE/src/

pfad="$PWD"
echo $pfad

if [ -d "$pfad/TauAnalysis" ]
then
    echo "The directory TauAnalysis already exists!!"
else 
    echo "Getting SVfitStandalone..." 
    git clone https://github.com/veelken/SVfit_standalone.git TauAnalysis/SVfitStandalone
    echo "Getting SVFitTF (TransferFunctions)..."
    git clone https://github.com/veelken/SVfitTF.git TauAnalysis/SVfitTF
    echo "switch branch..."
    cd $CMSSW_BASE/src/TauAnalysis/SVfitStandalone
    git checkout HIG-16-006
fi

echo "You are now in:" 
cd $startingpath
echo $startingpath
