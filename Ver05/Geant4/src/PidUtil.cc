/*
  PidUtil.cc

  2012/5  K.Shirotori
*/

#include "PidUtil.hh"

const G4String PName[] = {
  "geantino", // 0
  "pi-", // 1
  "pi+", // 2
  "kaon-", // 3
  "kaon+", // 4
  "mu-", // 5
  "mu+", // 6
  "e-", // 7
  "e+", // 8
  "gamma", //9 
  "proton", // 10
  "neutron", // 11 
  "pi0", // 12
  "nu_e", // 13
  "anti_proton", // 14
  "lambda", // 15
  "kaon0", // 16
  "anti_kaon0", // 17
  "kaon0S", // 18
  "kaon0L", // 19
  "sigma0", // 20
  "sigma-", // 21
  "sigma+", // 22
  "xi-", // 23
  "xi0", // 24
  "omega-", // 25
  "eta", // 26 
  "eta_prime", // 27
  "D-", // 28
  "D0", // 29
};

const int MaxNum = sizeof(PName)/sizeof(G4String);
const int DefaultNum = -1;  

int PIDFunc( const G4String & name )
{
  for( int i=0; i<MaxNum; ++i ){
    if( name==PName[i] ) return i;
  }

  if( name=="nu_e" || name=="nu_mu" ||
      name=="anti_nu_e" || name=="anti_nu_mu" )
    return 13;
  else if( name=="kaonP1" ) return 4;
  else if( name=="kaonM1" ) return 3;
  else if( name=="pionP1" ) return 2;
  else if( name=="pionM1" ) return 1;
  else if( name=="DmesonM1" ) return 28;
  else if( name=="DmesonZ1" ) return 29;
  else 
    return DefaultNum;
}

const G4String & PIDParticleName( int Pid )
{
  if(Pid==DefaultNum || Pid>=MaxNum )
    return PName[0];
  else
    return PName[Pid];
}
