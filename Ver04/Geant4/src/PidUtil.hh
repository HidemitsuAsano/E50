/*
  PidUtil.hh

  2012/5  K.Shirotori
*/

#ifndef PidUtil_h
#define PidUtil_h 1 

#include "G4String.hh"

int PIDFunc( const G4String & name );

const G4String & PIDParticleName( int Pid );

#endif
