/*
  GeneratorHelper.hh

  2012/5  K.Shirotori
*/

#ifndef GeneratorHelper_h 
#define GeneratorHelper_h 1

#include "G4ThreeVector.hh"

class DetectorConstruction;
class PrimaryGeneratorAction;
class G4ParticleGun;
class AngDisGenerator;

G4ThreeVector UniformDirectionInUV( double u0, double v0, 
				    double hu=0, double hv=0 );

G4ThreeVector GaussDirectionInUV( double u0, double v0, 
				  double su, double sv );

G4ThreeVector UniformPosition( double hx, double hy, double hz ); 

G4ThreeVector GaussPosition( double sx, double sy, double hz );

G4ThreeVector UniformDirectionInThetaPhi( double cost1, double cost2 );

double BreitWigner( double mean, double gamma );

struct BeamParam {
  BeamParam( DetectorConstruction *det, PrimaryGeneratorAction *pact );
  int pid;
  double x0,y0,z0,u0,v0,p0;
  double dx,dy,dz,du,dv,dp, dztarg;
  double cost1, cost2;
  G4ThreeVector gTargetPos;
  double TargetRotAngH;
  double RecMass;
  double RecWidth;
  G4ParticleGun *gun;
};

bool Decay2Body( double Mini, double Mf1, double Mf2, 
		 const G4ThreeVector & Pini, 
		 G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
		 G4ThreeVector & Pf1CM,  G4ThreeVector & Pf2CM,
		 const AngDisGenerator & generator ); 

bool Scattering2Body( double Mi1, double Mi2, double Mf1, double Mf2,
		      const G4ThreeVector & Pini, 
		      G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
		      G4ThreeVector & Pf1CM,  G4ThreeVector & Pf2CM,
		      const AngDisGenerator & generator );  

bool Decay3BodyPhaseSpace( double Mini, double Mf1, double Mf2, double Mf3,
			   const G4ThreeVector & Pini,
			   G4ThreeVector & Pf1,  G4ThreeVector & Pf2, 
			   G4ThreeVector & Pf3 ); 

bool Decay3BodyPhaseSpaceCM( double Mini, double Mf1, double Mf2, double Mf3,
			     const G4ThreeVector & Pini,
			     G4ThreeVector & Pf1,  G4ThreeVector & Pf2, 
			     G4ThreeVector & Pf3,
			     G4ThreeVector & Pf1CM,  G4ThreeVector & Pf2CM, 
			     G4ThreeVector & Pf3CM ); 

bool Scattering3Body( double Mi1, double Mi2, 
		      double Mf1, double Mf2, double Mf3,
		      const G4ThreeVector & Pini, 
		      G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
		      G4ThreeVector & Pf3,
		      const AngDisGenerator & generator );  

bool Decay4BodyPhaseSpace( double Mini, 
			   double Mf1, double Mf2, 
			   double Mf3, double Mf4,
			   const G4ThreeVector & Pini,
			   G4ThreeVector & Pf1, G4ThreeVector & Pf2, 
			   G4ThreeVector & Pf3, G4ThreeVector & Pf4 ); 

bool Decay5BodyPhaseSpace( double Mini, 
			   double Mf1, double Mf2, 
			   double Mf3, double Mf4, 
			   double Mf5,
			   const G4ThreeVector & Pini,
			   G4ThreeVector & Pf1, G4ThreeVector & Pf2, 
			   G4ThreeVector & Pf3, G4ThreeVector & Pf4,
			   G4ThreeVector & Pf5 ); 

bool Decay6BodyPhaseSpace( double Mini, 
			   double Mf1, double Mf2, 
			   double Mf3, double Mf4, 
			   double Mf5, double Mf6,
			   const G4ThreeVector & Pini,
			   G4ThreeVector & Pf1, G4ThreeVector & Pf2, 
			   G4ThreeVector & Pf3, G4ThreeVector & Pf4,
			   G4ThreeVector & Pf5, G4ThreeVector & Pf6 ); 

G4double Pmx(G4double e, G4double p1, G4double p2);

#endif
