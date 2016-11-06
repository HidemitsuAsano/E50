#!/bin/sh
cd simCharm/Geant4;
nice time ~/simCharm/Geant4/bin/simSpec ~/simCharm/Geant4/mac/run2.mac ~/simCharm/Geant4/conf/LambdacDSCheck20/sim.conf.LambdacDSCheck20_$1_$2 ~/simCharm/Geant4/histo_Sim/test_LambdacDSCheck20_$1_$2.root
