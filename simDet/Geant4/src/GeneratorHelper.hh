/*
  GeneratorHelper.hh

  2016/4  K.Shirotori
*/

#ifndef GeneratorHelper_h 
#define GeneratorHelper_h 1

#include "G4ThreeVector.hh"

class PrimaryGeneratorAction;
class G4ParticleGun;
class AngDisGenerator;

G4ThreeVector UniformDirectionInUV( double u0, double v0, 
				    double hu=0, double hv=0 );

G4ThreeVector GaussDirectionInUV( double u0, double v0, 
				  double su, double sv );

G4ThreeVector UniformPosition( double hx, double hy, double hz ); 

G4ThreeVector GaussPosition( double sx, double sy, double hz );

struct BeamParam {
  BeamParam( PrimaryGeneratorAction *pact );
  int pid;
  double x0,y0,z0,u0,v0,p0;
  double dx,dy,dz,du,dv,dp;
  G4ParticleGun *gun;
};

#endif
