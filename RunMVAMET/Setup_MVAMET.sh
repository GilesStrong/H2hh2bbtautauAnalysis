#!/bin/bash

startingpath="$PWD"
cd $CMSSW_BASE/src/

if [ -d "$CMSSW_BASE/src/RecoMET/" ] || ["$CMSSW_BASE/src/DataFormats/" ]
then
    echo "Your environment is already set!"
fi

if [ -d "$CMSSW_BASE/src" ]
then
    echo "The /src directory already exists!!"
    mv $CMSSW_BASE/src $CMSSW_BASE/tmp_src
    mkdir $CMSSW_BASE/src
    cd $CMSSW_BASE/src
    cmsenv
    git cms-addpkg RecoMET/METPUSubtraction
    git cms-addpkg DataFormats/METReco
    git remote add -f mvamet https://github.com/rfriese/cmssw.git
    git checkout mvamet/mvamet80 -b mvamet
    mkdir RecoMET/METPUSubtraction/data
    cd RecoMET/METPUSubtraction/data
    wget https://github.com/rfriese/cmssw/raw/MVAMET2_beta_0.6/RecoMET/METPUSubtraction/data/weightfile.root
    cd $CMSSW_BASE/src
    mv $CMSSW_BASE/tmp_src/* .
    rm -r $CMSSW_BASE/tmp_src
    scram b -j12
else 
    cmsenv
    git cms-addpkg RecoMET/METPUSubtraction
    git cms-addpkg DataFormats/METReco
    git remote add -f mvamet https://github.com/rfriese/cmssw.git
    git checkout mvamet/mvamet80 -b mvamet
    mkdir RecoMET/METPUSubtraction/data
    cd RecoMET/METPUSubtraction/data
    wget https://github.com/rfriese/cmssw/raw/MVAMET2_beta_0.6/RecoMET/METPUSubtraction/data/weightfile.root
    cd $CMSSW_BASE/src
    scram b -j 10
fi

echo "You are now in:" 
cd $startingpath
echo $startingpath

