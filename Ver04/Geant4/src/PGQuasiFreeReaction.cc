/*
  PGQuasiFreeReaction.cc

  2010/11/12 T.Takahashi
*/

#include "PGQuasiFreeReaction.hh"

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "Analysis.hh"
#include "PrimaryInfo.hh"
#include "GeneratorHelper.hh"
#include "AngDisGenerator.hh"

#include "G4ios.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTypes.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "G4LorentzVector.hh"


const int MaxTry=100;

// Fermi Gas Model
G4ThreeVector PGQuasiFreeReaction::generateNucleonMomentumFG()
{
  double cost=2.*G4UniformRand()-1.;
  double sint=sqrt(1.-cost*cost);
  double phi=2.*acos(-1.)*G4UniformRand();
  double k3=kF_*kF_*kF_*G4UniformRand();
  double k=pow(k3,1./3.);
  
  return G4ThreeVector( k*sint*cos(phi), k*sint*sin(phi), k*cost );
}


// pi- Nucleus -> K+ Sigma- Nucleus'
void PGQuasiFreeReaction::
SigmaMinusQuasiFreeProduction( G4Event *anEvent, bool generateSigma )
{
  // Caution
  // Uniform angular distribution in CM now

  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4SigmaMinus::Definition()->GetPDGMass();

  BeamParam *BP = pGen->BP_;
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  AGUniform gen;


  G4ThreeVector LPi1=pb*LBeamDir;
  double LEi1=sqrt(pb*pb+Mi1*Mi1); 
  G4LorentzVector LVi1(LPi1,LEi1);

  G4ThreeVector LPi2=generateNucleonMomentumFG();
  double LEi2=sqrt(LPi2.mag2()+Mi2*Mi2);
  G4LorentzVector LVi2(LPi2,LEi2);
    
  G4ThreeVector beta1( LPi2/LEi2 );
  LVi1.boost(-beta1); LVi2.boost(-beta1);
  G4ThreeVector LmPi1=LVi1.vect(), LmPi2=LVi2.vect();
  //  std::cout << "LmPi1=" << LmPi1/GeV  << " GeV/c  LmPi2=" 
  //	    << LmPi2/GeV << " GeV/c" << std::endl;
  
  G4ThreeVector LmPf1, LmPf2;
  bool status1=Scattering2Body( Mi1, Mi2, Mf1, Mf2, LmPi1, 
				LmPf1, LmPf2, gen );
  int n1=0;
  while (!status1){
    if(++n1>MaxTry){
      G4Exception("SigmaMinusQuasiFreeProduction:");
    }
    LPi2=generateNucleonMomentumFG(); LEi2=sqrt(LPi2.mag2()+Mi2*Mi2);
    LVi1=G4LorentzVector(LPi1,LEi1);  LVi2=G4LorentzVector(LPi2,LEi2); 
    beta1=LPi2/LEi2;
    LVi1.boost(-beta1); LVi2.boost(-beta1);
    LmPi1=LVi1.vect(); LmPi2=LVi2.vect();
    //    std::cout << "LmPi1=" << LmPi1/GeV  << " GeV/c  LmPi2=" 
    //	      << LmPi2/GeV << " GeV/c" << std::endl;
    status1=Scattering2Body( Mi1, Mi2, Mf1, Mf2, LmPi1, 
			     LmPf1, LmPf2, gen );
  }
  G4LorentzVector LVf1(LmPf1,sqrt(LmPf1.mag2()+Mf1*Mf1)),
    LVf2(LmPf2,sqrt(LmPf2.mag2()+Mf2*Mf2));
  LVf1.boost(beta1); LVf2.boost(beta1);
  G4ThreeVector LPf1=LVf1.vect(), LPf2=LVf2.vect();
  double LEf1=LVf1.e(), LEf2=LVf2.e();
  double Q=LEi1+TargM_-LEf1-LEf2-ResM_;
  G4ThreeVector LPres=LPi1-LPf1-LPf2;
  double KEres=0.5*LPres.mag2()/ResM_;
  
  //  std::cout << "Q=" << Q/GeV << " GeV KEres=" << KEres/GeV 
  //	    << " GeV" << std::endl;


  G4ThreeVector gPos(LPos), gPf1(LPf1);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  pGen->gunKP_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunKP_->SetParticleMomentum( gPf1 );
  pGen->gunKP_->GeneratePrimaryVertex( anEvent );

  if( generateSigma ){
    G4ThreeVector gPf2(LPf2);
    gPf2.rotateX(90.0*degree); gPf2.rotateZ( BP->TargetRotAngH );
    pGen->gunSigmaM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunSigmaM_->SetParticleMomentum( gPf2 );
    pGen->gunSigmaM_->GeneratePrimaryVertex( anEvent );
  }

  if( pGen->anaMan_ ){
    pGen->SetPrimaryInformation( pGen->gunKP_, LPos, LPf1, pb*LBeamDir );
  }
  // SigmaMinusQuasiFreeProduction
  // pi- Nucleus -> K+ Sigma- Nucleus'
}
