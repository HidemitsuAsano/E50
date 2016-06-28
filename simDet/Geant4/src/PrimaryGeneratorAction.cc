/*
  PrimaryGeneratorAction.cc

  2016/1  K.Shirotori

  2016/5 modified by H. Asano
*/

#include "PrimaryGeneratorAction.hh"
#include "Analysis.hh"
#include "PrimaryInfo.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "PidUtil.hh"
#include "GeneratorHelper.hh"

#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "Randomize.hh"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <zlib.h>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

PrimaryGeneratorAction::
PrimaryGeneratorAction( Analysis *analysisManager )
  : G4VUserPrimaryGeneratorAction(), 
    anaMan_(analysisManager), BP_(0)
{
  MakeGuns();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  if(BP_) delete BP_;
  DeleteGuns();
}

void PrimaryGeneratorAction::GeneratePrimaries( G4Event *anEvent )
{
  ConfMan *confMan = ConfMan::GetConfManager();
  G4int reactionMode = confMan->ReactionMode();
  G4int nparticle = confMan->GetNumberOfParticles();

  if(!BP_) BP_ = new BeamParam( this );

  switch(reactionMode){
  case 1:
    for(int ipart=0;ipart<nparticle;ipart++){
      GPBeamGaussInUVGaussInXY( anEvent ); //break;
      // GPBeamGaussInUVGaussInXY( anEvent ); break;
    }
    break;
  case 2:
    GPBeamGaussInUVUniformInXY( anEvent ); break;
  case 3:
    GPBeamUniformInUVUniformInXY( anEvent ); break;
  case 4:
    GPBeamUniformInUVGaussInXY( anEvent ); break;
  default:
    GPBeamGaussInUVGaussInXY( anEvent ); break;
  }
}

void PrimaryGeneratorAction::MakeGuns()
{
  gunDefault_ = new G4ParticleGun( G4Geantino::Definition() );
  gunPiM_     = new G4ParticleGun( G4PionMinus::Definition() );
  gunPiP_     = new G4ParticleGun( G4PionPlus::Definition() );
  gunPi0_     = new G4ParticleGun( G4PionZero::Definition() );

  gunKM_      = new G4ParticleGun( G4KaonMinus::Definition() );
  gunKP_      = new G4ParticleGun( G4KaonPlus::Definition() );
  gunK0_      = new G4ParticleGun( G4KaonZero::Definition() );
  gunK0bar_   = new G4ParticleGun( G4AntiKaonZero::Definition() );
  gunK0S_     = new G4ParticleGun( G4KaonZeroShort::Definition() );
  gunK0L_     = new G4ParticleGun( G4KaonZeroLong::Definition() );

  gunGamma_   = new G4ParticleGun( G4Gamma::Definition() ); 
  gunEM_      = new G4ParticleGun( G4Electron::Definition() );
  gunEP_      = new G4ParticleGun( G4Positron::Definition() );

  gunMuM_     = new G4ParticleGun( G4MuonMinus::Definition() );
  gunMuP_     = new G4ParticleGun( G4MuonPlus::Definition() );

  gunProton_  = new G4ParticleGun( G4Proton::Definition() );
  gunNeutron_ = new G4ParticleGun( G4Neutron::Definition() );
  gunProtonb_ = new G4ParticleGun( G4AntiProton::Definition() );
  gunNeutronb_= new G4ParticleGun( G4AntiNeutron::Definition() );
  gunLambda_  = new G4ParticleGun( G4Lambda::Definition() );
  gunSigmaP_  = new G4ParticleGun( G4SigmaPlus::Definition() );
  gunSigmaZ_  = new G4ParticleGun( G4SigmaZero::Definition() );
  gunSigmaM_  = new G4ParticleGun( G4SigmaMinus::Definition() );
  gunXiM_     = new G4ParticleGun( G4XiMinus::Definition() );
  gunXiZ_     = new G4ParticleGun( G4XiZero::Definition() );

  gunDMesonM_  = new G4ParticleGun( G4DMesonMinus::Definition() );
  gunDMesonZ_  = new G4ParticleGun( G4DMesonZero::Definition() );

  gunLambdacP_ = new G4ParticleGun( G4LambdacPlus::Definition() );
  gunSigmacPP_ = new G4ParticleGun( G4SigmacPlusPlus::Definition() );
  gunSigmacZ_  = new G4ParticleGun( G4SigmacZero::Definition() );

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *DMinus1, *DZero1, *DZero2, *DZero3, *DZero4;
  G4ParticleDefinition *PiM1, *KM1, *PiP1, *KP1;
  DMinus1 = particleTable->FindParticle("DmesonM1");
  DZero4  = particleTable->FindParticle("DmesonZ4");
  DZero3  = particleTable->FindParticle("DmesonZ3");
  DZero2  = particleTable->FindParticle("DmesonZ2");
  DZero1  = particleTable->FindParticle("DmesonZ1");
  PiP1    = particleTable->FindParticle("pionP1");
  PiM1    = particleTable->FindParticle("pionM1");
  KP1     = particleTable->FindParticle("kaonP1");
  KM1     = particleTable->FindParticle("kaonM1");
  gunDMesonM1_  = new G4ParticleGun( DMinus1 );
  gunDMesonZ4_  = new G4ParticleGun( DZero4 );
  gunDMesonZ3_  = new G4ParticleGun( DZero3 );
  gunDMesonZ2_  = new G4ParticleGun( DZero2 );
  gunDMesonZ1_  = new G4ParticleGun( DZero1 );
  gunPiP1_      = new G4ParticleGun( PiP1 );
  gunPiM1_      = new G4ParticleGun( PiM1 );
  gunKP1_       = new G4ParticleGun( KP1 );
  gunKM1_       = new G4ParticleGun( KM1 );
}

void PrimaryGeneratorAction::DeleteGuns()
{
  delete gunSigmacZ_; gunSigmacZ_=0;
  delete gunSigmacPP_; gunSigmacPP_=0;
  delete gunLambdacP_; gunLambdacP_=0;

  delete gunDMesonZ4_; gunDMesonZ4_=0;
  delete gunDMesonZ3_; gunDMesonZ3_=0;
  delete gunDMesonZ2_; gunDMesonZ2_=0;
  delete gunDMesonZ1_; gunDMesonZ1_=0;
  delete gunDMesonM1_; gunDMesonM1_=0;

  delete gunDMesonZ_; gunDMesonZ_=0;
  delete gunDMesonM_; gunDMesonM_=0;

  delete gunXiZ_; gunXiZ_=0;
  delete gunXiM_; gunXiM_=0;
  delete gunSigmaM_; gunSigmaM_=0;
  delete gunSigmaZ_; gunSigmaZ_=0;
  delete gunSigmaP_; gunSigmaP_=0;
  delete gunLambda_; gunLambda_=0;
  delete gunNeutronb_; gunNeutronb_=0;
  delete gunProtonb_; gunProtonb_=0;
  delete gunNeutron_; gunNeutron_=0;
  delete gunProton_; gunProton_=0;

  delete gunMuP_; gunMuP_=0; 
  delete gunMuM_; gunMuM_=0; 

  delete gunEP_; gunEP_=0;
  delete gunEM_; gunEM_=0;
  delete gunGamma_; gunGamma_=0;

  delete gunK0L_; gunK0L_=0;
  delete gunK0S_; gunK0S_=0;
  delete gunK0bar_; gunK0bar_=0; 
  delete gunK0_; gunK0_=0; 

  delete gunKP1_; gunKP1_=0;
  delete gunKM1_; gunKM1_=0;
  delete gunPiP1_; gunPiP1_=0;
  delete gunPiM1_; gunPiM1_=0;

  delete gunKP_; gunKP_=0;
  delete gunKM_; gunKM_=0;

  delete gunPi0_; gunPi0_=0;
  delete gunPiP_; gunPiP_=0;
  delete gunPiM_; gunPiM_=0;
}

G4ParticleGun * PrimaryGeneratorAction::chooseGun( G4int Pid )
{
  G4String name=PIDParticleName(Pid);
  return chooseGun( name );
}

G4ParticleGun * PrimaryGeneratorAction::chooseGun( const G4String & name  )
{
  if( name=="pi-" )             return gunPiM_;
  else if( name=="pi+" )        return gunPiP_;
  else if( name=="pi0" )        return gunPi0_;
  else if( name=="kaon-" )      return gunKM_;
  else if( name=="kaon+" )      return gunKP_;
  else if( name=="kaon0" )      return gunK0_;
  else if( name=="anti_kaon0" ) return gunK0bar_;
  else if( name=="kaon0S" )     return gunK0S_;
  else if( name=="kaon0L" )     return gunK0L_;
  else if( name=="gamma" )      return gunGamma_;
  else if( name=="e-" )         return gunEM_;
  else if( name=="e+" )         return gunEP_;
  else if( name=="mu-" )        return gunMuM_;
  else if( name=="mu+" )        return gunMuP_;
  else if( name=="proton" )     return gunProton_;
  else if( name=="neutron" )    return gunNeutron_;
  else if( name=="anti_proton" )  return gunProtonb_;
  else if( name=="anti_neutron" ) return gunNeutronb_;
  else if( name=="lambda" )     return gunLambda_;
  else if( name=="sigma+" )     return gunSigmaP_;
  else if( name=="sigma0" )     return gunSigmaZ_;
  else if( name=="sigma-" )     return gunSigmaM_;
  else if( name=="xi-" )        return gunXiM_;
  else if( name=="xi0" )        return gunXiZ_;
  else                          return gunDefault_;
}

void PrimaryGeneratorAction::GPBeamGaussInUVGaussInXY( G4Event *anEvent )
{
  ConfMan *confMan = ConfMan::GetConfManager(); 

  if(!BP_) BP_ = new BeamParam( this );

  G4ThreeVector LPos;
  LPos = G4ThreeVector( BP_->x0, BP_->y0, BP_->z0 )
    + GaussPosition( BP_->dx, BP_->dy, BP_->dz );

  G4ThreeVector LDir = 
    GaussDirectionInUV( BP_->u0, BP_->v0, BP_->du, BP_->dv );
  G4ThreeVector gPos = LPos, gDir = LDir;

  G4double p;
  if( confMan->MomBeam() ){
    G4int nev=anEvent->GetEventID();
    switch( nev%30 ){
    case 0:
      p=0.2*GeV; break;
    case 1:
      p=0.3*GeV; break;
    case 2:
      p=0.4*GeV; break;
    case 3:
      p=0.5*GeV; break;
    case 4:
      p=0.6*GeV; break;
    case 5:
      p=0.7*GeV; break;
    case 6:
      p=0.8*GeV; break;
    case 7:
      p=0.9*GeV; break;
    case 8:
      p=1.0*GeV; break;
    case 9:
      p=1.2*GeV; break;
    case 10:
      p=1.4*GeV; break;
    case 11:
      p=1.6*GeV; break;
    case 12:
      p=1.8*GeV; break;
    case 13:
      p=2.0*GeV; break;
    case 14:
      p=2.4*GeV; break;
    case 15:
      p=2.8*GeV; break;
    case 16:
      p=3.0*GeV; break;
    case 17:
      p=3.5*GeV; break;
    case 18:
      p=4.0*GeV; break;
    case 19:
      p=5.0*GeV; break;
    case 20:
      p=6.0*GeV; break;
    case 21:
      p=7.0*GeV; break;
    case 22:
      p=8.0*GeV; break;
    case 23:
      p=9.0*GeV; break;
    case 24:
      p=10.0*GeV; break;
    case 25:
      p=12.0*GeV; break;
    case 26:
      p=14.0*GeV; break;
    case 27:
      p=16.0*GeV; break;
    case 28:
      p=18.0*GeV; break;
    case 29:
      p=20.0*GeV; break;
    default:
      p=5.0*GeV; break;
    }
  }
  else{
    p = BP_->p0;
    if(BP_->dp!=0.0) p += (G4UniformRand()-0.5)*BP_->dp;
  }

  BP_->gun->SetParticlePosition( gPos );
  BP_->gun->SetParticleMomentum( p*gDir );
  BP_->gun->GeneratePrimaryVertex( anEvent );
  
  if( anaMan_ ){
    SetPrimaryInformation( LPos, 0.,0.,0.,0.,0.,0., p*LDir, 0. );
  }
  // GPBeamGaussInUVGausInXY
}

void PrimaryGeneratorAction::GPBeamGaussInUVUniformInXY( G4Event *anEvent )
{
  G4ThreeVector LPos = G4ThreeVector( BP_->x0, BP_->y0, BP_->z0 )
    + UniformPosition( BP_->dx, BP_->dy, BP_->dz );
  G4ThreeVector LDir = 
    GaussDirectionInUV( BP_->u0, BP_->v0, BP_->du, BP_->dv );
  G4ThreeVector gPos = LPos, gDir = LDir;

  G4double p = BP_->p0;
  if(BP_->dp!=0.0) p += (G4UniformRand()-0.5)*BP_->dp;

  BP_->gun->SetParticlePosition( gPos );
  BP_->gun->SetParticleMomentum( p*gDir );
  BP_->gun->GeneratePrimaryVertex( anEvent );

  if( anaMan_ ){
    SetPrimaryInformation( LPos, 0.,0.,0.,0.,0.,0., p*LDir, 0. );
  }
  // GPBeamGaussInUVUniformInXY
}

void PrimaryGeneratorAction::
GPBeamUniformInUVUniformInXY( G4Event *anEvent )
{
  G4ThreeVector LPos = G4ThreeVector( BP_->x0, BP_->y0, BP_->z0 )
    + UniformPosition( BP_->dx, BP_->dy, BP_->dz );
  G4ThreeVector LDir = 
    UniformDirectionInUV( BP_->u0, BP_->v0, BP_->du, BP_->dv );
  G4ThreeVector gPos = LPos, gDir = LDir;

  G4double p = BP_->p0;
  if(BP_->dp!=0.0) p += (G4UniformRand()-0.5)*BP_->dp;

  //  BP_->gun->SetParticleEnergy(0.0);
  BP_->gun->SetParticlePosition( gPos );
  BP_->gun->SetParticleMomentum( p*gDir );
  BP_->gun->GeneratePrimaryVertex( anEvent );

  //ConfMan *confMan = ConfMan::GetConfManager(); 
  //  G4cout << "ParticleGenerator: gun:" 
  //	 << BP_->gun->GetParticleDefinition()->GetParticleName() << G4endl;

  if( anaMan_ ){
    SetPrimaryInformation( LPos, 0.,0.,0.,0.,0.,0., p*LDir, 0.);
  }
  // GPBeamUniformInUVUniformInXY
}

void PrimaryGeneratorAction::GPBeamUniformInUVGaussInXY( G4Event *anEvent )
{
  if(!BP_) BP_ = new BeamParam( this );

  G4ThreeVector LPos = G4ThreeVector( BP_->x0, BP_->y0, BP_->z0 )
    + UniformPosition( BP_->dx, BP_->dy, BP_->dz );
  G4ThreeVector LDir = 
    GaussDirectionInUV( BP_->u0, BP_->v0, BP_->du, BP_->dv );
  G4ThreeVector gPos = LPos, gDir = LDir;

  G4double p = BP_->p0;
  if(BP_->dp!=0.0) p += (G4UniformRand()-0.5)*BP_->dp;

  BP_->gun->SetParticlePosition( gPos );
  BP_->gun->SetParticleMomentum( p*gDir );
  BP_->gun->GeneratePrimaryVertex( anEvent );
  
  if( anaMan_ ){
    SetPrimaryInformation( LPos, 0.,0.,0.,0.,0.,0., p*LDir, 0. );
  }
  // GPBeamUniformInUVGausInXY
}

void PrimaryGeneratorAction::
SetPrimaryInformation( const G4ThreeVector & pos,
		       G4double mass1, 
		       G4double mass2, 
		       const G4ThreeVector & mom1,
		       const G4ThreeVector & mom2,
		       const G4ThreeVector & mom1CM,
		       const G4ThreeVector & mom2CM,
		       const G4ThreeVector & beamMom,
		       const G4ThreeVector & mom3 )
{
  PrimaryInfo pI;
  pI.x=pos.x(); pI.y=pos.y(); pI.z=pos.z();
  pI.m1 = mass1; pI.m2 = mass2;
  pI.p1=mom1.mag();
  pI.p2=mom2.mag();
  pI.p1cm=mom1CM.mag();
  pI.p2cm=mom2CM.mag();
  pI.beta1=1./sqrt(1.+(mass1/mom1.mag())*(mass1/mom1.mag()));
  pI.beta2=1./sqrt(1.+(mass2/mom2.mag())*(mass2/mom2.mag()));

  if(mom1.mag()!=0.0){
    pI.theta1=mom1.theta(); pI.phi1=mom1.phi();
    pI.u1=mom1.x()/mom1.z(); pI.v1=mom1.y()/mom1.z();
  }
  else {
    pI.theta1=0.0; pI.phi1=0.0; pI.u1=0.0; pI.v1=0.0;
  }
  if(mom2.mag()!=0.0){
    pI.theta2=mom2.theta(); pI.phi2=mom2.phi();
    pI.u2=mom2.x()/mom2.z(); pI.v2=mom2.y()/mom2.z();
  }
  else {
    pI.theta2=0.0; pI.phi2=0.0; pI.u2=0.0; pI.v2=0.0;
  }
  if(mom1CM.mag()!=0.0){
    pI.theta1cm=mom1CM.theta(); pI.phi1cm=mom1CM.phi();
    pI.u1cm=mom1CM.x()/mom1CM.z(); pI.v1cm=mom1CM.y()/mom1CM.z();
  }
  else {
    pI.theta1cm=0.0; pI.phi1cm=0.0; pI.u1cm=0.0; pI.v1cm=0.0;
  }
  if(mom2CM.mag()!=0.0){
    pI.theta2cm=mom2CM.theta(); pI.phi2cm=mom2CM.phi();
    pI.u2cm=mom2CM.x()/mom2CM.z(); pI.v2cm=mom2CM.y()/mom2CM.z();
  }
  else {
    pI.theta2cm=0.0; pI.phi2cm=0.0; pI.u2cm=0.0; pI.v2cm=0.0;
  }

  pI.pbeam=beamMom.mag();
  if(beamMom.z()!=0.0){
    pI.ubeam=beamMom.x()/beamMom.z(); pI.vbeam=beamMom.y()/beamMom.z();
  }
  else{
    pI.ubeam=0.0; pI.vbeam=0.0;
  }

  pI.p3=mom3.mag();
  if(mom3.z()!=0.0){
    pI.u3=mom3.x()/mom3.z(); pI.v3=mom3.y()/mom3.z();
  }
  else{
    pI.u3=0.0; pI.v3=0.0;
  }
  anaMan_->PrimaryGeneration( &pI );
}

