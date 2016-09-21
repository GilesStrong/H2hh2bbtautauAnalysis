#!/bin/bash

startingpath="$PWD"

cd $CMSSW_BASE/src/

pfad="$PWD"
echo $pfad

if [ -d "$pfad/HHKinFit2" ]
then
    echo "The directory HHKinFit2 already exists!!"
else 
    echo "Getting HHKinFit2..." 
    git clone -b CMSSWversion https://github.com/bvormwald/HHKinFit2.git $CMSSW_BASE/src/HHKinFit2
fi

echo "You are now in:" 
cd $startingpath
echo $startingpath