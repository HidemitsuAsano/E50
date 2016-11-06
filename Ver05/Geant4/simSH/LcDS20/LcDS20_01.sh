#!/bin/sh
ssh -f highpc1 ./simCharm/Geant4/simSH/LcDS20/simLcDS20_1.sh 1 1 
ssh -f highpc2 ./simCharm/Geant4/simSH/LcDS20/simLcDS20_1.sh 1 2 
ssh -f highpc3 ./simCharm/Geant4/simSH/LcDS20/simLcDS20_1.sh 2 1 
ssh -f highpc4 ./simCharm/Geant4/simSH/LcDS20/simLcDS20_1.sh 2 2 
ssh -f highpc5 ./simCharm/Geant4/simSH/LcDS20/simLcDS20_1.sh 3 1 
ssh -f highpc6 ./simCharm/Geant4/simSH/LcDS20/simLcDS20_1.sh 3 2 

