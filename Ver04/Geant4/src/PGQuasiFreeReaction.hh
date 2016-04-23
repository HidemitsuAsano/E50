/*
  PGQuasiFreeReaction.hh

  2010/11/12 T.Takahashi
*/

#ifndef PGQuasiFreeReaction_h
#define PGQuasiFreeReaction_h 1

#include "G4ThreeVector.hh"

class PrimaryGeneratorAction;
class G4Event;

class PGQuasiFreeReaction
{
public:
  PGQuasiFreeReaction( PrimaryGeneratorAction *PGAction, 
		       double FermiMomentum, 
		       double TargMass, double ResidualMass )
    : pGen(PGAction), kF_(FermiMomentum), TargM_(TargMass),
    ResM_(ResidualMass)
  {}
  ~PGQuasiFreeReaction(){}

  void SigmaMinusQuasiFreeProduction( G4Event *anEvent,
				      bool generateSigma=false );

private:
  PrimaryGeneratorAction *pGen;
  double kF_, TargM_, ResM_;

  G4ThreeVector generateNucleonMomentumFG();

};

#endif
