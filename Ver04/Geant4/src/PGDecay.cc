/* 
  PGDecay.cc

  2012/5  K.Shirotori
*/

#include "PGDecay.hh"

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

// D- -> K+ pi- pi- 
void PGDecay::
DMesonMinus3BodyDecayKPPiMPiM( G4Event *anEvent, bool generatePiP )
{
  G4double Mini=G4DMesonMinus::Definition()->GetPDGMass();
  G4double Mf1 =G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2 =G4PionMinus::Definition()->GetPDGMass();
  G4double Mf3 =G4PionMinus::Definition()->GetPDGMass();

  BeamParam *BP = pGen->BP_;
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dz );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;
  G4ThreeVector LPf1, LPf2, LPf3;

  bool status=Decay3BodyPhaseSpace( Mini, Mf1, Mf2, Mf3, pb*LBeamDir,
				    LPf1, LPf2, LPf3 );
  if(!status){
    G4Exception("DMesonMinus3BodyDecayKPPiMPiM under Threshold");
  }
  G4ThreeVector gPos(LPos), gPf1(LPf1), gPf2(LPf2);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  gPf2.rotateX(90.0*degree); gPf2.rotateZ( BP->TargetRotAngH );

  pGen->gunKP_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunKP_->SetParticleMomentum( gPf1 );
  pGen->gunKP_->GeneratePrimaryVertex( anEvent );
  pGen->gunPiM_->SetParticleMomentum( gPf2 );
  pGen->gunPiM_->GeneratePrimaryVertex( anEvent );

  if( generatePiP ){
    G4ThreeVector gPf3(LPf3);
    gPf3.rotateX(90.0*degree); gPf3.rotateZ( BP->TargetRotAngH );
    pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM_->SetParticleMomentum( gPf3 );
    pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
  }
//   if( pGen->anaMan_ ){
//     pGen->SetPrimaryInformation( pGen->gunKP_, LPos, LPf1, pb*LBeamDir );
//     //pGen->SetPrimaryInformation( pGen->gunPiM_, LPos, LPf2, pb*LBeamDir );
//   }
  // DMesonMinus3BodyDecayKPPiMPiM
  // D- -> K+ pi- pi- 
}


// D0 -> K+ pi-
void PGDecay::
DMesonZero2BodyDecayKPPiM( G4Event *anEvent, bool generatePiP )
{
  G4double Mini=G4DMesonZero::Definition()->GetPDGMass();
  G4double Mf1 =G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2 =G4PionMinus::Definition()->GetPDGMass();
  BeamParam *BP = pGen->BP_;
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dz );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;
  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;
  AGUniform gen( BP->cost1, BP->cost2 );
  bool status = Decay2Body( Mini, Mf1, Mf2, pb*LBeamDir, LPf1, LPf2, LPf1CM, LPf2CM, gen );
  if(!status){
    G4Exception("DMesonZero2BodyDecayKPPiM: under Threshold");
  }
  G4ThreeVector gPos(LPos), gPf1(LPf1);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  pGen->gunKP_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunKP_->SetParticleMomentum( gPf1 );
  pGen->gunKP_->GeneratePrimaryVertex( anEvent );
  if( generatePiP ){
    G4ThreeVector gPf2(LPf2);
    gPf2.rotateX(90.0*degree); gPf2.rotateZ( BP->TargetRotAngH );
    pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM_->SetParticleMomentum( gPf2 );
    pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
  }    

//   if( pGen->anaMan_ ){
//     //pGen->SetPrimaryInformation( pGen->gunKP_, LPos, LPf1, pb*LBeamDir );
//     pGen->SetPrimaryInformation( pGen->gunPiM_, LPos, LPf2, pb*LBeamDir );
//   }
  // DMesonZero2BodyDecayKPPiM
  // D0 -> K+ pi-
}


// K- -> mu- anti_nu_mu
void PGDecay::
KMinus2BodyDecayMuMNu( G4Event *anEvent, bool generateNeutrino )
{
  G4double Mini=G4KaonMinus::Definition()->GetPDGMass();
  G4double Mf1 =G4MuonMinus::Definition()->GetPDGMass();
  G4double Mf2 =G4AntiNeutrinoMu::Definition()->GetPDGMass();
  BeamParam *BP = pGen->BP_;
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dz );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;
  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;
  AGUniform gen( BP->cost1, BP->cost2 );
  bool status = Decay2Body( Mini, Mf1, Mf2, pb*LBeamDir, LPf1, LPf2, LPf1CM, LPf2CM, gen );
  if(!status){
    G4Exception("KMinus2BodyDecayMuMNu: under Threshold");
  }
  G4ThreeVector gPos(LPos), gPf1(LPf1);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  pGen->gunMuM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunMuM_->SetParticleMomentum( gPf1 );
  pGen->gunMuM_->GeneratePrimaryVertex( anEvent );

//   if( pGen->anaMan_ ){
//     pGen->SetPrimaryInformation( pGen->gunMuM_, LPos, LPf1, pb*LBeamDir );
//   }
  // KMinus2BodyDecayMuMNu
  // K- -> mu- anti_nu_mu
}

// K- -> pi- pi0
void PGDecay::
KMinus2BodyDecayPiMPi0( G4Event *anEvent, bool generatePi0 )
{
  G4double Mini=G4KaonMinus::Definition()->GetPDGMass();
  G4double Mf1 =G4PionMinus::Definition()->GetPDGMass();
  G4double Mf2 =G4PionZero::Definition()->GetPDGMass();
  BeamParam *BP = pGen->BP_;
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dz );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;
  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;
  AGUniform gen( BP->cost1, BP->cost2 );
  bool status = Decay2Body( Mini, Mf1, Mf2, pb*LBeamDir, LPf1, LPf2, LPf1CM, LPf2CM, gen );
  if(!status){
    G4Exception("KMinus2BodyDecayMuMPi0: under Threshold");
  }
  G4ThreeVector gPos(LPos), gPf1(LPf1);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM_->SetParticleMomentum( gPf1 );
  pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
  if( generatePi0 ){
    G4ThreeVector gPf2(LPf2);
    gPf2.rotateX(90.0*degree); gPf2.rotateZ( BP->TargetRotAngH );
    pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPi0_->SetParticleMomentum( gPf2 );
    pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  }    


//   if( pGen->anaMan_ ){
//     pGen->SetPrimaryInformation( pGen->gunPiM_, LPos, LPf1, pb*LBeamDir );
//   }
  // KMinus2BodyDecayPiMPi0
  // K- -> pi- pi0
}

// K- -> pi- pi- pi+ 
void PGDecay::
KMinus3BodyDecayPiMPiMPiP( G4Event *anEvent, bool generatePiP )
{
  G4double Mini=G4KaonMinus::Definition()->GetPDGMass();
  G4double Mf1 =G4PionMinus::Definition()->GetPDGMass();
  G4double Mf2 =G4PionMinus::Definition()->GetPDGMass();
  G4double Mf3 =G4PionPlus::Definition()->GetPDGMass();

  BeamParam *BP = pGen->BP_;
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dz );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;
  G4ThreeVector LPf1, LPf2, LPf3;

  bool status=Decay3BodyPhaseSpace( Mini, Mf1, Mf2, Mf3, pb*LBeamDir,
				    LPf1, LPf2, LPf3 );
  if(!status){
    G4Exception("KMinus3BodyDecayPiMPiMPiP: under Threshold");
  }
  G4ThreeVector gPos(LPos), gPf1(LPf1), gPf2(LPf2);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  gPf2.rotateX(90.0*degree); gPf2.rotateZ( BP->TargetRotAngH );

  pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM_->SetParticleMomentum( gPf1 );
  pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
  pGen->gunPiM_->SetParticleMomentum( gPf2 );
  pGen->gunPiM_->GeneratePrimaryVertex( anEvent );

  if( generatePiP ){
    G4ThreeVector gPf3(LPf3);
    gPf3.rotateX(90.0*degree); gPf3.rotateZ( BP->TargetRotAngH );
    pGen->gunPiP_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiP_->SetParticleMomentum( gPf3 );
    pGen->gunPiP_->GeneratePrimaryVertex( anEvent );
  }
//   if( pGen->anaMan_ ){
//     pGen->SetPrimaryInformation( pGen->gunPiM_, LPos, LPf1, pb*LBeamDir );
//   }
  // KMinus3BodyDecayPiMPiMPiP
  // K- -> pi- pi- pi+
}
			
// K- -> e- pi0 anti_nu_e
void PGDecay::
KMinus3BodyDecayEMPi0Nu( G4Event *anEvent, 
			 bool generatePi0, bool genetateNeutrino )
{
  G4double Mini=G4KaonMinus::Definition()->GetPDGMass();
  G4double Mf1 =G4Electron::Definition()->GetPDGMass();
  G4double Mf2 =G4PionZero::Definition()->GetPDGMass();
  G4double Mf3 =G4AntiNeutrinoE::Definition()->GetPDGMass();

  BeamParam *BP = pGen->BP_;
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dz );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;
  G4ThreeVector LPf1, LPf2, LPf3;

  bool status=Decay3BodyPhaseSpace( Mini, Mf1, Mf2, Mf3, pb*LBeamDir,
				    LPf1, LPf2, LPf3 );
  if(!status){
    G4Exception("KMinus3BodyDecayEMPi0Nu: under Threshold");
  }
  G4ThreeVector gPos(LPos), gPf1(LPf1);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  pGen->gunEM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunEM_->SetParticleMomentum( gPf1 );
  pGen->gunEM_->GeneratePrimaryVertex( anEvent );

  if( generatePi0 ){
    G4ThreeVector gPf2(LPf2);
    gPf2.rotateX(90.0*degree); gPf2.rotateZ( BP->TargetRotAngH );
    pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPi0_->SetParticleMomentum( gPf2 );
    pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  }

//   if( pGen->anaMan_ ){
//     pGen->SetPrimaryInformation( pGen->gunEM_, LPos, LPf1, pb*LBeamDir );
//   }
  // KMinus3BodyDecayEMPi0Nu
  // K- -> e- pi0 anti_nu_e
}

// K- -> mu- pi0 anti_nu_mu
void PGDecay::
KMinus3BodyDecayMuMPi0Nu( G4Event *anEvent,
			  bool generatePi0, bool genetateNeutrino )
{    
  G4double Mini=G4KaonMinus::Definition()->GetPDGMass();
  G4double Mf1 =G4MuonMinus::Definition()->GetPDGMass();
  G4double Mf2 =G4PionZero::Definition()->GetPDGMass();
  G4double Mf3 =G4AntiNeutrinoMu::Definition()->GetPDGMass();

  BeamParam *BP = pGen->BP_;
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dz );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;
  G4ThreeVector LPf1, LPf2, LPf3;

  bool status=Decay3BodyPhaseSpace( Mini, Mf1, Mf2, Mf3, pb*LBeamDir,
				    LPf1, LPf2, LPf3 );
  if(!status){
    G4Exception("KMinus3BodyDecayMuMPi0Nu: under Threshold");
  }
  G4ThreeVector gPos(LPos), gPf1(LPf1);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  pGen->gunMuM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunMuM_->SetParticleMomentum( gPf1 );
  pGen->gunMuM_->GeneratePrimaryVertex( anEvent );

  if( generatePi0 ){
    G4ThreeVector gPf2(LPf2);
    gPf2.rotateX(90.0*degree); gPf2.rotateZ( BP->TargetRotAngH );
    pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPi0_->SetParticleMomentum( gPf2 );
    pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  }

//   if( pGen->anaMan_ ){
//     pGen->SetPrimaryInformation( pGen->gunMuM_, LPos, LPf1, pb*LBeamDir );
//   }
  // KMinus3BodyDecayMuMPi0Nu
  // K- -> mu- pi0 anti_nu_mu
}

// K- > pi- pi0 pi0
void PGDecay::
KMinus3BodyDecayPiMPi0Pi0( G4Event *anEvent, bool generatePi0 )
{
  G4double Mini=G4KaonMinus::Definition()->GetPDGMass();
  G4double Mf1 =G4PionMinus::Definition()->GetPDGMass();
  G4double Mf2 =G4PionZero::Definition()->GetPDGMass();
  G4double Mf3 =G4PionZero::Definition()->GetPDGMass();

  BeamParam *BP = pGen->BP_;
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dz );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;
  G4ThreeVector LPf1, LPf2, LPf3;

  bool status=Decay3BodyPhaseSpace( Mini, Mf1, Mf2, Mf3, pb*LBeamDir,
				    LPf1, LPf2, LPf3 );
  if(!status){
    G4Exception("KMinus3BodyDecayPiMPi0Pi0: under Threshold");
  }
  G4ThreeVector gPos(LPos), gPf1(LPf1);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM_->SetParticleMomentum( gPf1 );
  pGen->gunPiM_->GeneratePrimaryVertex( anEvent );

  if( generatePi0 ){
    G4ThreeVector gPf2(LPf2), gPf3(LPf3);
    gPf2.rotateX(90.0*degree); gPf2.rotateZ( BP->TargetRotAngH );
    gPf3.rotateX(90.0*degree); gPf3.rotateZ( BP->TargetRotAngH );
    pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPi0_->SetParticleMomentum( gPf2 );
    pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
    pGen->gunPi0_->SetParticleMomentum( gPf3 );
    pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  }

//   if( pGen->anaMan_ ){
//     pGen->SetPrimaryInformation( pGen->gunPiM_, LPos, LPf1, pb*LBeamDir );
//   }
  // KMinus3BodyDecayPiMPi0Pi0
  // K- > pi- pi0 pi0
}
