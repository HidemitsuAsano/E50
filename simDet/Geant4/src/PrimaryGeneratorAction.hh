/*
  PrimaryGeneratorAction.hh

  2016/4  K.Shirotori
*/

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include <vector>

class Analysis;
class G4ParticleGun;
class BeamParam;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction( Analysis *anaManager=0 );
  ~PrimaryGeneratorAction();

public:
  void GeneratePrimaries( G4Event *anEvent );

protected:
  Analysis *anaMan_;
  BeamParam *BP_;

  G4ParticleGun *chooseGun( int Pid );
  G4ParticleGun *chooseGun( const G4String & name );

  void MakeGuns();
  void DeleteGuns();

private:
  G4ParticleGun *gunDefault_;
  G4ParticleGun *gunPiM_, *gunPiP_, *gunPi0_;
  G4ParticleGun *gunPiM1_, *gunPiP1_;
  G4ParticleGun *gunKM_, *gunKP_;
  G4ParticleGun *gunKM1_, *gunKP1_;
  G4ParticleGun *gunK0_, *gunK0bar_, *gunK0S_, *gunK0L_;
  G4ParticleGun *gunGamma_, *gunEM_, *gunEP_;
  G4ParticleGun *gunMuM_, *gunMuP_;
  G4ParticleGun *gunProton_, *gunNeutron_;
  G4ParticleGun *gunProtonb_, *gunNeutronb_;
  G4ParticleGun *gunLambda_;
  G4ParticleGun *gunSigmaP_, *gunSigmaZ_, *gunSigmaM_;
  G4ParticleGun *gunXiM_, *gunXiZ_;

  G4ParticleGun *gunDMesonM_, *gunDMesonZ_;
  G4ParticleGun *gunDMesonstarM_;
  G4ParticleGun *gunDMesonM1_;
  G4ParticleGun *gunDMesonZ1_, *gunDMesonZ2_, *gunDMesonZ3_, *gunDMesonZ4_;
  G4ParticleGun *gunLambdacP_, *gunSigmacPP_, *gunSigmacZ_;

  void GPBeamGaussInUVGaussInXY( G4Event *anEvent );
  void GPBeamGaussInUVUniformInXY( G4Event *anEvent );
  void GPBeamUniformInUVUniformInXY( G4Event *anEvent );
  void GPBeamUniformInUVGaussInXY( G4Event *anEvent );

  void SetPrimaryInformation( const G4ThreeVector & pos,
			      double mass1, 
			      double mass2, 
			      const G4ThreeVector & mom1,
			      const G4ThreeVector & mom2,
			      const G4ThreeVector & mom1CM,
			      const G4ThreeVector & mom2CM,
			      const G4ThreeVector & beamMom,
			      const G4ThreeVector & pisMom );
  
  void GPBeamDemonstration( G4Event *anEvent );

  friend class BeamParam;
};

#endif
