/*
  PGReaction5.cc 

  2013/12  K.Shirotori
*/

////Xi production
////Both Pion and Kaon beam conditions

#include "PGReaction5.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "Analysis.hh"
#include "PrimaryInfo.hh"
#include "GeneratorHelper.hh"
#include "AngDisGenerator.hh"
#include "ParticleID.hh"
#include "ConfMan.hh"

#include "G4ios.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTypes.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtStdlibRandomEngine.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenExternal/EvtExternalGenList.hh"

#include <iostream>
#include <string>
#include <list>

const G4int MaxTry=100;

////////////////////////////////////////
////////////////////
//50: pi- p -> Xi*- K*0 K+ : FS Xi- K+ K+ pi- 
void PGReaction5::XiKStarMesonFromPiMinusProton1( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  // 

  BeamParam *BP = pGen->BP_;

  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  G4double Mf1=0.8960*GeV;//K*0 mass
  G4double Mf2=BP->RecMass;//Xi- mass
  //G4double Mf2=G4XiMinus::Definition()->GetPDGMass();//Xi- mass
  G4double Mf3=G4KaonPlus::Definition()->GetPDGMass();//K+ mass

  G4double Wid = 0.05030*GeV;
  static EvtId KS1 = EvtPDL::getId(std::string("K*0"));  

  // G4double Mass=BP->RecMass;
  // G4double Width=BP->RecWidth;
  // G4double Mf2=BreitWigner( Mass, Width );//Xi-* mass

  G4double MF1=BreitWigner( Mf1, Wid );

   G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4double Ei=sqrt( Mi1*Mi1+Mi2*Mi2+2.*Mi2*sqrt(Mi1*Mi1+pb*pb) );
  G4double Mmin=(0.4940+0.1396)*GeV;

  while(1){
    MF1=BreitWigner( Mf1, Wid );
    if(MF1>Mmin) break;
  }

  //Decay 
  G4ThreeVector LPf1, LPf2, LPf3;
  G4ThreeVector LPf1CM, LPf2CM;

  G4bool status=
    Decay3BodyPhaseSpace( Ei, MF1, Mf2, Mf3, pb*LBeamDir, LPf1, LPf2, LPf3 );

  G4int n=0;
  while(!status){
    if(++n>MaxTry){
      G4Exception("KStarMesonXiMinusKaonPlusFromPiMinusProton1: under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    while( !(Ei>MF1+Mf2+Mf3) ){
      while(1){
	MF1=BreitWigner( Mf1, Wid );
	if(MF1>Mmin) break;
      }
    }
    status = 
      Decay3BodyPhaseSpace( Ei, MF1, Mf2, Mf3, pb*LBeamDir, LPf1, LPf2, LPf3 );
  }

  G4ThreeVector gPos(LPos), gPf2(LPf2), gPf3(LPf3);

  //Beam
  pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM1_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );

  pGen->gunXiM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunXiM_->SetParticleMomentum( gPf2 );
  pGen->gunXiM_->GeneratePrimaryVertex( anEvent );
  
  pGen->gunKP1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunKP1_->SetParticleMomentum( gPf3 );
  pGen->gunKP1_->GeneratePrimaryVertex( anEvent );

  ///////K*0 decay by EvtGen
  EvtParticle* parent1(0);

  G4LorentzVector LvKS1;
  LvKS1.setVect(LPf1/GeV);
  LvKS1.setE(sqrt(MF1*MF1+LPf1.mag2())/GeV);
  EvtVector4R pInit1( LvKS1.e(),
		      LvKS1.vect().x(), LvKS1.vect().y(), LvKS1.vect().z() );

  parent1 = EvtParticleFactory::particleFactory(KS1, pInit1); 
  evtGenerator->generateDecay(parent1);

  static EvtId kp=EvtPDL::getId(std::string("K+"));
  static EvtId km=EvtPDL::getId(std::string("K-"));
  static EvtId pip=EvtPDL::getId(std::string("pi+"));
  static EvtId pim=EvtPDL::getId(std::string("pi-"));
  static EvtId pr=EvtPDL::getId(std::string("p+"));
  static EvtId apr=EvtPDL::getId(std::string("anti-p-"));
  static EvtId pi0=EvtPDL::getId(std::string("pi0"));
  static EvtId ne=EvtPDL::getId(std::string("n0"));
  static EvtId ane=EvtPDL::getId(std::string("anti-n0"));
  static EvtId kz=EvtPDL::getId(std::string("K0"));
  static EvtId ks=EvtPDL::getId(std::string("K_S0"));
  static EvtId kl=EvtPDL::getId(std::string("K_L0"));
  static EvtId gm=EvtPDL::getId(std::string("gamma"));
  static EvtId ep=EvtPDL::getId(std::string("e+"));
  static EvtId em=EvtPDL::getId(std::string("e-"));
  static EvtId mup=EvtPDL::getId(std::string("mu+"));
  static EvtId mum=EvtPDL::getId(std::string("mu-"));

  G4int nkp=0, nkm=0;
  G4int npip=0, npim=0;
  G4int nproton=0, nprotonb=0;
  G4int npiz=0;
  G4int nneutron=0, nneutronb=0;
  G4int nkz=0, nkzs=0, nkzl=0;
  G4int ngamma=0;
  G4int nep=0, nem=0;
  G4int nmup=0, nmum=0;

  EvtParticle *p1;
  p1=parent1->getDaug(0);

  do{
    EvtId type=p1->getId();
    
    if (type==kp)  nkp++;
    if (type==km)  nkm++;
    if (type==pip) npip++;
    if (type==pim) npim++;
    if (type==pr)  nproton++;
    if (type==apr) nprotonb++;
    if (type==pi0) npiz++;
    if (type==ne)  nneutron++;
    if (type==ane) nneutronb++;
    if (type==kz)  nkz++;
    if (type==ks)  nkzs++;
    if (type==kl)  nkzl++;
    if (type==gm)  ngamma++;
    if (type==ep)  nep++;
    if (type==em)  nem++;
    if (type==mup) nmup++;
    if (type==mum) nmum++;
    
    p1=p1->nextIter();
  }while(p1!=0);

  const G4int nKP=nkp;
  const G4int nKM=nkm;
  const G4int nPiP=npip;
  const G4int nPiM=npim;
  const G4int nP=nproton;
  const G4int nPb=nprotonb;
  const G4int nPiZ=npiz;
  const G4int nN=nneutron;
  const G4int nNb=nneutronb;
  const G4int nKZ=nkz;
  const G4int nKZS=nkzs;
  const G4int nKZL=nkzl;
  const G4int nG=ngamma;
  const G4int nEP=nep;
  const G4int nEM=nem;
  const G4int nMuP=nmup;
  const G4int nMuM=nmum;

  G4ThreeVector gPfKP[nKP], gPfKM[nKM];
  G4ThreeVector gPfPiP[nPiP], gPfPiM[nPiM];
  G4ThreeVector gPfP[nP], gPfPb[nPb];
  G4ThreeVector gPfPiZ[nPiZ];
  G4ThreeVector gPfN[nN], gPfNb[nNb];
  G4ThreeVector gPfKZ[nKZ];
  G4ThreeVector gPfKZS[nKZS], gPfKZL[nKZL];
  G4ThreeVector gPfG[nG];
  G4ThreeVector gPfEP[nEP], gPfEM[nEM];
  G4ThreeVector gPfMuP[nMuP], gPfMuM[nMuM];

  G4int ikp=0, ikm=0;
  G4int ipip=0, ipim=0;
  G4int ip=0, ipb=0;
  G4int ipiz=0;
  G4int in=0, inb=0;
  G4int ikz=0;
  G4int ikzs=0, ikzl=0;
  G4int igamma=0;
  G4int iep=0, iem=0;
  G4int imup=0, imum=0;

  EvtParticle *q1;
  q1=parent1->getDaug(0);

  do{
    EvtId type=q1->getId();
    EvtVector4R P4Lab=q1->getP4Lab();

    if( type == kp ){
      gPfKP[ikp].setX(P4Lab.get(1)*GeV);
      gPfKP[ikp].setY(P4Lab.get(2)*GeV);
      gPfKP[ikp].setZ(P4Lab.get(3)*GeV);
      ikp++; 
    }
    if( type == km ){
      gPfKM[ikm].setX(P4Lab.get(1)*GeV);
      gPfKM[ikm].setY(P4Lab.get(2)*GeV);
      gPfKM[ikm].setZ(P4Lab.get(3)*GeV);
      ikm++;
    }
    if( type == pip ){
      gPfPiP[ipip].setX(P4Lab.get(1)*GeV);
      gPfPiP[ipip].setY(P4Lab.get(2)*GeV);
      gPfPiP[ipip].setZ(P4Lab.get(3)*GeV);
      ipip++;
    }
    if( type == pim ){
      gPfPiM[ipim].setX(P4Lab.get(1)*GeV);
      gPfPiM[ipim].setY(P4Lab.get(2)*GeV);
      gPfPiM[ipim].setZ(P4Lab.get(3)*GeV);
      ipim++;
    }
    if( type == pr ){
      gPfP[ip].setX(P4Lab.get(1)*GeV);
      gPfP[ip].setY(P4Lab.get(2)*GeV);
      gPfP[ip].setZ(P4Lab.get(3)*GeV);
      ip++;
    }
    if( type == apr ){
      gPfPb[ipb].setX(P4Lab.get(1)*GeV);
      gPfPb[ipb].setY(P4Lab.get(2)*GeV);
      gPfPb[ipb].setZ(P4Lab.get(3)*GeV);
      ipb++;
    }
    if( type == pi0 ){
      gPfPiZ[ipiz].setX(P4Lab.get(1)*GeV);
      gPfPiZ[ipiz].setY(P4Lab.get(2)*GeV);
      gPfPiZ[ipiz].setZ(P4Lab.get(3)*GeV);
      ipiz++;
    }
    if( type == ne ){
      gPfN[in].setX(P4Lab.get(1)*GeV);
      gPfN[in].setY(P4Lab.get(2)*GeV);
      gPfN[in].setZ(P4Lab.get(3)*GeV);
      in++;
    }
    if( type == ane ){
      gPfNb[inb].setX(P4Lab.get(1)*GeV);
      gPfNb[inb].setY(P4Lab.get(2)*GeV);
      gPfNb[inb].setZ(P4Lab.get(3)*GeV);
      inb++;
    }
    if( type == kz ){
      gPfKZ[ikz].setX(P4Lab.get(1)*GeV);
      gPfKZ[ikz].setY(P4Lab.get(2)*GeV);
      gPfKZ[ikz].setZ(P4Lab.get(3)*GeV);
      ikz++;
    }
    if( type == ks ){
      gPfKZS[ikzs].setX(P4Lab.get(1)*GeV);
      gPfKZS[ikzs].setY(P4Lab.get(2)*GeV);
      gPfKZS[ikzs].setZ(P4Lab.get(3)*GeV);
      ikzs++;
    }
    if( type == kl ){
      gPfKZL[ikzl].setX(P4Lab.get(1)*GeV);
      gPfKZL[ikzl].setY(P4Lab.get(2)*GeV);
      gPfKZL[ikzl].setZ(P4Lab.get(3)*GeV);
      ikzl++;
    }
    if( type == gm ){
      gPfG[igamma].setX(P4Lab.get(1)*GeV);
      gPfG[igamma].setY(P4Lab.get(2)*GeV);
      gPfG[igamma].setZ(P4Lab.get(3)*GeV);
      igamma++;
    }
    if( type == ep ){
      gPfEP[iep].setX(P4Lab.get(1)*GeV);
      gPfEP[iep].setY(P4Lab.get(2)*GeV);
      gPfEP[iep].setZ(P4Lab.get(3)*GeV);
      iep++;
    }
    if( type == em ){
      gPfEM[iem].setX(P4Lab.get(1)*GeV);
      gPfEM[iem].setY(P4Lab.get(2)*GeV);
      gPfEM[iem].setZ(P4Lab.get(3)*GeV);
      iem++;
    }
    if( type == mup ){
      gPfMuP[imup].setX(P4Lab.get(1)*GeV);
      gPfMuP[imup].setY(P4Lab.get(2)*GeV);
      gPfMuP[imup].setZ(P4Lab.get(3)*GeV);
      imup++;
    }
    if( type == mum ){
      gPfMuM[imum].setX(P4Lab.get(1)*GeV);
      gPfMuM[imum].setY(P4Lab.get(2)*GeV);
      gPfMuM[imum].setZ(P4Lab.get(3)*GeV);
      imum++;
    }

    q1=q1->nextIter();
  }while(q1!=0);

 // G4cout<< "**********************" << G4endl;
  // G4cout<< ikp << " " << ikm << G4endl;
  // for( G4int i=0; i<ikp; i++) G4cout<<  gPfKP[i] << G4endl;
  // for( G4int i=0; i<ikm; i++) G4cout<<  gPfKM[i] << G4endl;
  // G4cout<< ipip << " " << ipim << G4endl;
  // for( G4int i=0; i<ipip; i++) G4cout<<  gPfPiP[i] << G4endl;
  // for( G4int i=0; i<ipim; i++) G4cout<<  gPfPiM[i] << G4endl;
  // G4cout<< ip << " " << ipb <<G4endl;
  // for( G4int i=0; i<ip; i++)  G4cout<<  gPfP[i] << G4endl;
  // for( G4int i=0; i<ipb; i++) G4cout<<  gPfPb[i] << G4endl;

  //SetGuns
  for( G4int i=0; i<ikp; i++){
    pGen->gunKP1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKP1_->SetParticleMomentum( gPfKP[i] );
    pGen->gunKP1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikm; i++){
    pGen->gunKM1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKM1_->SetParticleMomentum( gPfKM[i] );
    pGen->gunKM1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipip; i++){
    pGen->gunPiP1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiP1_->SetParticleMomentum( gPfPiP[i] );
    pGen->gunPiP1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipim; i++){
    pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM1_->SetParticleMomentum( gPfPiM[i] );
    pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ip; i++){
    pGen->gunProton_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunProton_->SetParticleMomentum( gPfP[i] );
    pGen->gunProton_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipb; i++){
    pGen->gunProtonb_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunProtonb_->SetParticleMomentum( gPfPb[i] );
    pGen->gunProtonb_->GeneratePrimaryVertex( anEvent );
  }
  // for( G4int i=0; i<ipiz; i++){
  //   pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunPi0_->SetParticleMomentum( gPfPiZ[i] );
  //   pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<in; i++){
  //   pGen->gunNeutron_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunNeutron_->SetParticleMomentum( gPfN[i] );
  //   pGen->gunNeutron_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<inb; i++){
  //   pGen->gunNeutronb_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunNeutronb_->SetParticleMomentum( gPfNb[i] );
  //   pGen->gunNeutronb_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikz; i++){
  //   pGen->gunK0_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0_->SetParticleMomentum( gPfKZ[i] );
  //   pGen->gunK0_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikzs; i++){
  //   pGen->gunK0S_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0S_->SetParticleMomentum( gPfKZS[i] );
  //   pGen->gunK0S_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikzl; i++){
  //   pGen->gunK0L_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0L_->SetParticleMomentum( gPfKZL[i] );
  //   pGen->gunK0L_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<igamma; i++){
  //   pGen->gunGamma_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunGamma_->SetParticleMomentum( gPfG[i] );
  //   pGen->gunGamma_->GeneratePrimaryVertex( anEvent );
  // }
  for( G4int i=0; i<iep; i++){
    pGen->gunEP_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunEP_->SetParticleMomentum( gPfEP[i] );
    pGen->gunEP_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<iem; i++){
    pGen->gunEM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunEM_->SetParticleMomentum( gPfEM[i] );
    pGen->gunEM_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<imup; i++){
    pGen->gunMuP_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunMuP_->SetParticleMomentum( gPfMuP[i] );
    pGen->gunMuP_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<imum; i++){
    pGen->gunMuM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunMuM_->SetParticleMomentum( gPfMuM[i] );
    pGen->gunMuM_->GeneratePrimaryVertex( anEvent );
  }

  if( pGen->anaMan_ ){
    pGen->SetPrimaryInformation( LPos, MF1, Mf2, 
                                 LPf1, LPf2, LPf1CM, LPf2CM, pb*LBeamDir, gPf3 );
  } 

  parent1->deleteTree();

  //XiKStarMesonFromPiMinusProton1
  //50: pi- p -> Xi*- K*0 : FS Xi- K+ K+ pi-
}


////////////////////
//51: pi- p -> Xi*- K*0 K* from Y* production : FS Xi- K+ K+ pi-
void PGReaction5::XiKStarMesonFromPiMinusProton2( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  //
  BeamParam *BP = pGen->BP_;

  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  G4double Mf1=0.8960*GeV;//K*0 mass
  G4double Wid=0.05030*GeV;
  static EvtId KS1 = EvtPDL::getId(std::string("K*0"));  

  //Y* mass ranging 2.0-3.0 GeV/c  
  G4double Mm1=(2.5+(G4UniformRand()-0.5))*GeV;

  G4double Mf2=BP->RecMass;//Xi- mass
  //G4double Mf2=G4XiMinus::Definition()->GetPDGMass();//Xi- mass
  G4double Mf3=G4KaonPlus::Definition()->GetPDGMass();//K+ mass

  // G4double Mass=BP->RecMass;
  // G4double Width=BP->RecWidth;
  // G4double Mf2=BreitWigner( Mass, Width );//Xi-* mass
  
  G4double MF1=BreitWigner( Mf1, Wid );

  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4double Ei=sqrt( Mi1*Mi1+Mi2*Mi2+2.*Mi2*sqrt(Mi1*Mi1+pb*pb) );
  G4double Mmin=(0.4940+0.1396)*GeV;

  while(1){
    MF1=BreitWigner( Mf1, Wid );
    if(MF1>Mmin) break;
  }

  G4ThreeVector LPf1, LPm1;
  G4ThreeVector LPf1CM, LPm1CM;
  AGUniform gen1( BP->cost1, BP->cost2 ); 
  //AGForward2 gen1( BP->cost1, BP->cost2 );
  G4bool status1 = Scattering2Body( Mi1, Mi2, MF1, Mm1, pb*LBeamDir,
				    LPf1, LPm1, LPm1CM, LPm1CM, gen1 );
  G4int n1=0;
  while(!status1){
    if(++n1>MaxTry){
      G4Exception("LambdacPlusDStartMesonMinusFromPiMinusProton: under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;
    
    while( !(Ei>MF1+Mm1) ){
      while(1){
	Mm1=(2.5+(G4UniformRand()-0.5))*GeV;
	MF1=BreitWigner( Mf1, Wid );
	if(MF1>Mmin) break;
      }
    }
    status1 = Scattering2Body( Mi1, Mi2, MF1, Mm1, pb*LBeamDir,
			       LPf1, LPm1, LPf1CM, LPm1CM, gen1 );
  }

  //Y* decay 
  G4ThreeVector LPf2, LPf3;
  G4ThreeVector LPf2CM, LPf3CM;
  AGUniform gen2;
  G4bool status2 = Decay2Body( Mm1, Mf2, Mf3, LPm1, LPf2, LPf3, LPf2CM, LPf3CM, gen2 );
  if(!status2){
    G4Exception("YStarMinus2BodyDecayXiKaonPlus under Threshold");
  }

  G4ThreeVector gPos(LPos), gPf2(LPf2), gPf3(LPf3);

  //Beam
  pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM1_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );

  pGen->gunXiM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunXiM_->SetParticleMomentum( gPf2 );
  pGen->gunXiM_->GeneratePrimaryVertex( anEvent );
  
  pGen->gunKP1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunKP1_->SetParticleMomentum( gPf3 );
  pGen->gunKP1_->GeneratePrimaryVertex( anEvent );

  ///////K*0 decay by EvtGen
  EvtParticle* parent1(0);

  G4LorentzVector LvKS1;
  LvKS1.setVect(LPf1/GeV);
  LvKS1.setE(sqrt(MF1*MF1+LPf1.mag2())/GeV);
  EvtVector4R pInit1( LvKS1.e(),
		      LvKS1.vect().x(), LvKS1.vect().y(), LvKS1.vect().z() );

  parent1 = EvtParticleFactory::particleFactory(KS1, pInit1); 
  evtGenerator->generateDecay(parent1);

  static EvtId kp=EvtPDL::getId(std::string("K+"));
  static EvtId km=EvtPDL::getId(std::string("K-"));
  static EvtId pip=EvtPDL::getId(std::string("pi+"));
  static EvtId pim=EvtPDL::getId(std::string("pi-"));
  static EvtId pr=EvtPDL::getId(std::string("p+"));
  static EvtId apr=EvtPDL::getId(std::string("anti-p-"));
  static EvtId pi0=EvtPDL::getId(std::string("pi0"));
  static EvtId ne=EvtPDL::getId(std::string("n0"));
  static EvtId ane=EvtPDL::getId(std::string("anti-n0"));
  static EvtId kz=EvtPDL::getId(std::string("K0"));
  static EvtId ks=EvtPDL::getId(std::string("K_S0"));
  static EvtId kl=EvtPDL::getId(std::string("K_L0"));
  static EvtId gm=EvtPDL::getId(std::string("gamma"));
  static EvtId ep=EvtPDL::getId(std::string("e+"));
  static EvtId em=EvtPDL::getId(std::string("e-"));
  static EvtId mup=EvtPDL::getId(std::string("mu+"));
  static EvtId mum=EvtPDL::getId(std::string("mu-"));

  G4int nkp=0, nkm=0;
  G4int npip=0, npim=0;
  G4int nproton=0, nprotonb=0;
  G4int npiz=0;
  G4int nneutron=0, nneutronb=0;
  G4int nkz=0, nkzs=0, nkzl=0;
  G4int ngamma=0;
  G4int nep=0, nem=0;
  G4int nmup=0, nmum=0;

  EvtParticle *p1;
  p1=parent1->getDaug(0);

  do{
    EvtId type=p1->getId();
    
    if (type==kp)  nkp++;
    if (type==km)  nkm++;
    if (type==pip) npip++;
    if (type==pim) npim++;
    if (type==pr)  nproton++;
    if (type==apr) nprotonb++;
    if (type==pi0) npiz++;
    if (type==ne)  nneutron++;
    if (type==ane) nneutronb++;
    if (type==kz)  nkz++;
    if (type==ks)  nkzs++;
    if (type==kl)  nkzl++;
    if (type==gm)  ngamma++;
    if (type==ep)  nep++;
    if (type==em)  nem++;
    if (type==mup) nmup++;
    if (type==mum) nmum++;
    
    p1=p1->nextIter();
  }while(p1!=0);

  const G4int nKP=nkp;
  const G4int nKM=nkm;
  const G4int nPiP=npip;
  const G4int nPiM=npim;
  const G4int nP=nproton;
  const G4int nPb=nprotonb;
  const G4int nPiZ=npiz;
  const G4int nN=nneutron;
  const G4int nNb=nneutronb;
  const G4int nKZ=nkz;
  const G4int nKZS=nkzs;
  const G4int nKZL=nkzl;
  const G4int nG=ngamma;
  const G4int nEP=nep;
  const G4int nEM=nem;
  const G4int nMuP=nmup;
  const G4int nMuM=nmum;

  G4ThreeVector gPfKP[nKP], gPfKM[nKM];
  G4ThreeVector gPfPiP[nPiP], gPfPiM[nPiM];
  G4ThreeVector gPfP[nP], gPfPb[nPb];
  G4ThreeVector gPfPiZ[nPiZ];
  G4ThreeVector gPfN[nN], gPfNb[nNb];
  G4ThreeVector gPfKZ[nKZ];
  G4ThreeVector gPfKZS[nKZS], gPfKZL[nKZL];
  G4ThreeVector gPfG[nG];
  G4ThreeVector gPfEP[nEP], gPfEM[nEM];
  G4ThreeVector gPfMuP[nMuP], gPfMuM[nMuM];

  G4int ikp=0, ikm=0;
  G4int ipip=0, ipim=0;
  G4int ip=0, ipb=0;
  G4int ipiz=0;
  G4int in=0, inb=0;
  G4int ikz=0;
  G4int ikzs=0, ikzl=0;
  G4int igamma=0;
  G4int iep=0, iem=0;
  G4int imup=0, imum=0;

  EvtParticle *q1;
  q1=parent1->getDaug(0);

  do{
    EvtId type=q1->getId();
    EvtVector4R P4Lab=q1->getP4Lab();

    if( type == kp ){
      gPfKP[ikp].setX(P4Lab.get(1)*GeV);
      gPfKP[ikp].setY(P4Lab.get(2)*GeV);
      gPfKP[ikp].setZ(P4Lab.get(3)*GeV);
      ikp++; 
    }
    if( type == km ){
      gPfKM[ikm].setX(P4Lab.get(1)*GeV);
      gPfKM[ikm].setY(P4Lab.get(2)*GeV);
      gPfKM[ikm].setZ(P4Lab.get(3)*GeV);
      ikm++;
    }
    if( type == pip ){
      gPfPiP[ipip].setX(P4Lab.get(1)*GeV);
      gPfPiP[ipip].setY(P4Lab.get(2)*GeV);
      gPfPiP[ipip].setZ(P4Lab.get(3)*GeV);
      ipip++;
    }
    if( type == pim ){
      gPfPiM[ipim].setX(P4Lab.get(1)*GeV);
      gPfPiM[ipim].setY(P4Lab.get(2)*GeV);
      gPfPiM[ipim].setZ(P4Lab.get(3)*GeV);
      ipim++;
    }
    if( type == pr ){
      gPfP[ip].setX(P4Lab.get(1)*GeV);
      gPfP[ip].setY(P4Lab.get(2)*GeV);
      gPfP[ip].setZ(P4Lab.get(3)*GeV);
      ip++;
    }
    if( type == apr ){
      gPfPb[ipb].setX(P4Lab.get(1)*GeV);
      gPfPb[ipb].setY(P4Lab.get(2)*GeV);
      gPfPb[ipb].setZ(P4Lab.get(3)*GeV);
      ipb++;
    }
    if( type == pi0 ){
      gPfPiZ[ipiz].setX(P4Lab.get(1)*GeV);
      gPfPiZ[ipiz].setY(P4Lab.get(2)*GeV);
      gPfPiZ[ipiz].setZ(P4Lab.get(3)*GeV);
      ipiz++;
    }
    if( type == ne ){
      gPfN[in].setX(P4Lab.get(1)*GeV);
      gPfN[in].setY(P4Lab.get(2)*GeV);
      gPfN[in].setZ(P4Lab.get(3)*GeV);
      in++;
    }
    if( type == ane ){
      gPfNb[inb].setX(P4Lab.get(1)*GeV);
      gPfNb[inb].setY(P4Lab.get(2)*GeV);
      gPfNb[inb].setZ(P4Lab.get(3)*GeV);
      inb++;
    }
    if( type == kz ){
      gPfKZ[ikz].setX(P4Lab.get(1)*GeV);
      gPfKZ[ikz].setY(P4Lab.get(2)*GeV);
      gPfKZ[ikz].setZ(P4Lab.get(3)*GeV);
      ikz++;
    }
    if( type == ks ){
      gPfKZS[ikzs].setX(P4Lab.get(1)*GeV);
      gPfKZS[ikzs].setY(P4Lab.get(2)*GeV);
      gPfKZS[ikzs].setZ(P4Lab.get(3)*GeV);
      ikzs++;
    }
    if( type == kl ){
      gPfKZL[ikzl].setX(P4Lab.get(1)*GeV);
      gPfKZL[ikzl].setY(P4Lab.get(2)*GeV);
      gPfKZL[ikzl].setZ(P4Lab.get(3)*GeV);
      ikzl++;
    }
    if( type == gm ){
      gPfG[igamma].setX(P4Lab.get(1)*GeV);
      gPfG[igamma].setY(P4Lab.get(2)*GeV);
      gPfG[igamma].setZ(P4Lab.get(3)*GeV);
      igamma++;
    }
    if( type == ep ){
      gPfEP[iep].setX(P4Lab.get(1)*GeV);
      gPfEP[iep].setY(P4Lab.get(2)*GeV);
      gPfEP[iep].setZ(P4Lab.get(3)*GeV);
      iep++;
    }
    if( type == em ){
      gPfEM[iem].setX(P4Lab.get(1)*GeV);
      gPfEM[iem].setY(P4Lab.get(2)*GeV);
      gPfEM[iem].setZ(P4Lab.get(3)*GeV);
      iem++;
    }
    if( type == mup ){
      gPfMuP[imup].setX(P4Lab.get(1)*GeV);
      gPfMuP[imup].setY(P4Lab.get(2)*GeV);
      gPfMuP[imup].setZ(P4Lab.get(3)*GeV);
      imup++;
    }
    if( type == mum ){
      gPfMuM[imum].setX(P4Lab.get(1)*GeV);
      gPfMuM[imum].setY(P4Lab.get(2)*GeV);
      gPfMuM[imum].setZ(P4Lab.get(3)*GeV);
      imum++;
    }

    q1=q1->nextIter();
  }while(q1!=0);

 // G4cout<< "**********************" << G4endl;
  // G4cout<< ikp << " " << ikm << G4endl;
  // for( G4int i=0; i<ikp; i++) G4cout<<  gPfKP[i] << G4endl;
  // for( G4int i=0; i<ikm; i++) G4cout<<  gPfKM[i] << G4endl;
  // G4cout<< ipip << " " << ipim << G4endl;
  // for( G4int i=0; i<ipip; i++) G4cout<<  gPfPiP[i] << G4endl;
  // for( G4int i=0; i<ipim; i++) G4cout<<  gPfPiM[i] << G4endl;
  // G4cout<< ip << " " << ipb <<G4endl;
  // for( G4int i=0; i<ip; i++)  G4cout<<  gPfP[i] << G4endl;
  // for( G4int i=0; i<ipb; i++) G4cout<<  gPfPb[i] << G4endl;

  //SetGuns
  for( G4int i=0; i<ikp; i++){
    pGen->gunKP1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKP1_->SetParticleMomentum( gPfKP[i] );
    pGen->gunKP1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikm; i++){
    pGen->gunKM1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKM1_->SetParticleMomentum( gPfKM[i] );
    pGen->gunKM1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipip; i++){
    pGen->gunPiP1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiP1_->SetParticleMomentum( gPfPiP[i] );
    pGen->gunPiP1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipim; i++){
    pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM1_->SetParticleMomentum( gPfPiM[i] );
    pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ip; i++){
    pGen->gunProton_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunProton_->SetParticleMomentum( gPfP[i] );
    pGen->gunProton_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipb; i++){
    pGen->gunProtonb_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunProtonb_->SetParticleMomentum( gPfPb[i] );
    pGen->gunProtonb_->GeneratePrimaryVertex( anEvent );
  }
  // for( G4int i=0; i<ipiz; i++){
  //   pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunPi0_->SetParticleMomentum( gPfPiZ[i] );
  //   pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<in; i++){
  //   pGen->gunNeutron_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunNeutron_->SetParticleMomentum( gPfN[i] );
  //   pGen->gunNeutron_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<inb; i++){
  //   pGen->gunNeutronb_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunNeutronb_->SetParticleMomentum( gPfNb[i] );
  //   pGen->gunNeutronb_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikz; i++){
  //   pGen->gunK0_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0_->SetParticleMomentum( gPfKZ[i] );
  //   pGen->gunK0_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikzs; i++){
  //   pGen->gunK0S_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0S_->SetParticleMomentum( gPfKZS[i] );
  //   pGen->gunK0S_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikzl; i++){
  //   pGen->gunK0L_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0L_->SetParticleMomentum( gPfKZL[i] );
  //   pGen->gunK0L_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<igamma; i++){
  //   pGen->gunGamma_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunGamma_->SetParticleMomentum( gPfG[i] );
  //   pGen->gunGamma_->GeneratePrimaryVertex( anEvent );
  // }
  for( G4int i=0; i<iep; i++){
    pGen->gunEP_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunEP_->SetParticleMomentum( gPfEP[i] );
    pGen->gunEP_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<iem; i++){
    pGen->gunEM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunEM_->SetParticleMomentum( gPfEM[i] );
    pGen->gunEM_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<imup; i++){
    pGen->gunMuP_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunMuP_->SetParticleMomentum( gPfMuP[i] );
    pGen->gunMuP_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<imum; i++){
    pGen->gunMuM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunMuM_->SetParticleMomentum( gPfMuM[i] );
    pGen->gunMuM_->GeneratePrimaryVertex( anEvent );
  }

  if( pGen->anaMan_ ){
    pGen->SetPrimaryInformation( LPos, MF1, Mf2, 
                                 LPf1, LPf2, LPf1CM, LPf2CM, pb*LBeamDir, gPf3 );
  } 

  parent1->deleteTree();

  //XiKStarMesonFromPiMinusProton2
  //pi- p -> Xi*- K*0 from Y* production : FS Xi- K+ K+ pi-
}


////////////////////
//K- p -> Xi*- K+
void PGReaction5::XiKStarMesonFromKMinusProton1( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  // 

  BeamParam *BP = pGen->BP_;

  G4double Mi1=G4KaonMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  //K+
  G4double Mf1=0.493677*GeV;//K+ mass
  G4double Wid = 0.0*GeV;
  static EvtId KS1 = EvtPDL::getId(std::string("K+"));  
  G4double MF1=BreitWigner( Mf1, Wid );

  //Xi
  G4double Mass=BP->RecMass;
  G4double Width=BP->RecWidth;

  static EvtId Xi1;
  if( 1.320*GeV<Mass && Mass<1.330*GeV ) Xi1 = EvtPDL::getId(std::string("Xi-")); 
  else if( 1.530*GeV<Mass && Mass<1.540*GeV ) Xi1 = EvtPDL::getId(std::string("Xi*-")); 
  else if( 1.815*GeV<Mass && Mass<1.825*GeV ) Xi1 = EvtPDL::getId(std::string("Xi(1820)-")); 
  // else if( 1.595*GeV<Mass && Mass<1.605*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1600)0")); 
  // else if( 1.665*GeV<Mass && Mass<1.675*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1670)0")); 
  // else if( 1.685*GeV<Mass && Mass<1.695*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1690)0")); 
  // else if( 1.795*GeV<Mass && Mass<1.805*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1800)0")); 
  // else if( 1.805*GeV<Mass && Mass<1.815*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1810)0")); 
  // else if( 1.815*GeV<Mass && Mass<1.825*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1820)0")); 
  // else if( 1.825*GeV<Mass && Mass<1.835*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1830)0")); 
  // else if( 1.885*GeV<Mass && Mass<1.895*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1890)0")); 
  else{
    Mass = 1.32171*GeV;
    Width = 0.0;
    Xi1 = EvtPDL::getId(std::string("Xi-")); 
  }

  G4double MF2=BreitWigner( Mass, Width );

  //Beam
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4double Ei=sqrt( Mi1*Mi1+Mi2*Mi2+2.*Mi2*sqrt(Mi1*Mi1+pb*pb) );
 
  G4double MKmin=0.493;
  G4double MYmin=(1.322+0.140)*GeV;//Xi + pi

  if( Mass<1.322*GeV ){
    while(1){
      MF1=BreitWigner( Mf1, Wid );
      if(MF1>MKmin) break;
    }
  }
  else{
    while(1){
      MF1=BreitWigner( Mf1, Wid );
      MF2=BreitWigner( Mass, Width );
      if(MF1>MKmin && MF2>MYmin) break;
    }
  }

  //Reaction
  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;
  AGUniform gen1( BP->cost1, BP->cost2 ); 
  //AGForward2 gen1( BP->cost1, BP->cost2 );

  G4bool status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
				    LPf1, LPf2, LPf1CM, LPf2CM, gen1 );
  G4int n=0;
  while(!status1){
    if(++n>MaxTry){
      G4Exception("XiKStarMesonFromKMinusProton1: under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    while( !(Ei>MF1+MF2) ){
      if( Mass<1.322*GeV ){
	while(1){
	  MF1=BreitWigner( Mf1, Wid );
	  if(MF1>MKmin) break;
	}
      }
      else{
	while(1){
	  MF1=BreitWigner( Mf1, Wid );
	  MF2=BreitWigner( Mass, Width );
	  if(MF1>MKmin && MF2>MYmin) break;
	}
      }
    }
    status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, gen1 );
  }

  G4ThreeVector gPos(LPos), gPf1(LPf1);

  //Beam
  pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM1_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );

  pGen->gunKP1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunKP1_->SetParticleMomentum( gPf1 );
  pGen->gunKP1_->GeneratePrimaryVertex( anEvent );

  ///////K+ && Xi* decay by EvtGen
  //EvtParticle* parent1(0);
  EvtParticle* parent2(0);
  
  // G4LorentzVector LvKS;
  // LvKS.setVect(LPf1/GeV);
  // LvKS.setE(sqrt(MF1*MF1+LPf1.mag2())/GeV);
  
  G4LorentzVector LvXi;
  LvXi.setVect(LPf2/GeV);
  LvXi.setE(sqrt(MF2*MF2+LPf2.mag2())/GeV);
  
  // EvtVector4R pInit1( LvKS.e(),
  // 		      LvKS.vect().x(), LvKS.vect().y(), LvKS.vect().z() );
  EvtVector4R pInit2( LvXi.e(),
		      LvXi.vect().x(), LvXi.vect().y(), LvXi.vect().z() );

  // parent1 = EvtParticleFactory::particleFactory(KS1, pInit1); 
  // evtGenerator->generateDecay(parent1);
  parent2 = EvtParticleFactory::particleFactory(Xi1, pInit2); 
  evtGenerator->generateDecay(parent2);

  //std::cout<< "******" << std::endl;

  static EvtId kp=EvtPDL::getId(std::string("K+"));
  static EvtId km=EvtPDL::getId(std::string("K-"));
  static EvtId pip=EvtPDL::getId(std::string("pi+"));
  static EvtId pim=EvtPDL::getId(std::string("pi-"));
  static EvtId pr=EvtPDL::getId(std::string("p+"));
  static EvtId apr=EvtPDL::getId(std::string("anti-p-"));
  static EvtId pi0=EvtPDL::getId(std::string("pi0"));
  static EvtId ne=EvtPDL::getId(std::string("n0"));
  static EvtId ane=EvtPDL::getId(std::string("anti-n0"));
  static EvtId kz=EvtPDL::getId(std::string("K0"));
  static EvtId ks=EvtPDL::getId(std::string("K_S0"));
  static EvtId kl=EvtPDL::getId(std::string("K_L0"));
  static EvtId gm=EvtPDL::getId(std::string("gamma"));
  static EvtId ep=EvtPDL::getId(std::string("e+"));
  static EvtId em=EvtPDL::getId(std::string("e-"));
  static EvtId mup=EvtPDL::getId(std::string("mu+"));
  static EvtId mum=EvtPDL::getId(std::string("mu-"));

  G4int nkp=0, nkm=0;
  G4int npip=0, npim=0;
  G4int nproton=0, nprotonb=0;
  G4int npiz=0;
  G4int nneutron=0, nneutronb=0;
  G4int nkz=0, nkzs=0, nkzl=0;
  G4int ngamma=0;
  G4int nep=0, nem=0;
  G4int nmup=0, nmum=0;

  //EvtParticle *p1;
  EvtParticle *p2;
  //p1=parent1->getDaug(0);
  p2=parent2->getDaug(0);

  // do{
  //   EvtId type=p1->getId();
    
  //   if (type==kp)  nkp++;
  //   if (type==km)  nkm++;
  //   if (type==pip) npip++;
  //   if (type==pim) npim++;
  //   if (type==pr)  nproton++;
  //   if (type==apr) nprotonb++;
  //   if (type==pi0) npiz++;
  //   if (type==ne)  nneutron++;
  //   if (type==ane) nneutronb++;
  //   if (type==kz)  nkz++;
  //   if (type==ks)  nkzs++;
  //   if (type==kl)  nkzl++;
  //   if (type==gm)  ngamma++;
  //   if (type==ep)  nep++;
  //   if (type==em)  nem++;
  //   if (type==mup) nmup++;
  //   if (type==mum) nmum++;
    
  //   p1=p1->nextIter();
  // }while(p1!=0);

  do{
    EvtId type=p2->getId();
    
    if (type==kp)  nkp++;
    if (type==km)  nkm++;
    if (type==pip) npip++;
    if (type==pim) npim++;
    if (type==pr)  nproton++;
    if (type==apr) nprotonb++;
    if (type==pi0) npiz++;
    if (type==ne)  nneutron++;
    if (type==ane) nneutronb++;
    if (type==kz)  nkz++;
    if (type==ks)  nkzs++;
    if (type==kl)  nkzl++;
    if (type==gm)  ngamma++;
    if (type==ep)  nep++;
    if (type==em)  nem++;
    if (type==mup) nmup++;
    if (type==mum) nmum++;
    
    p2=p2->nextIter();
  }while(p2!=0);

  const G4int nKP=nkp;
  const G4int nKM=nkm;
  const G4int nPiP=npip;
  const G4int nPiM=npim;
  const G4int nP=nproton;
  const G4int nPb=nprotonb;
  const G4int nPiZ=npiz;
  const G4int nN=nneutron;
  const G4int nNb=nneutronb;
  const G4int nKZ=nkz;
  const G4int nKZS=nkzs;
  const G4int nKZL=nkzl;
  const G4int nG=ngamma;
  const G4int nEP=nep;
  const G4int nEM=nem;
  const G4int nMuP=nmup;
  const G4int nMuM=nmum;

  G4ThreeVector gPfKP[nKP], gPfKM[nKM];
  G4ThreeVector gPfPiP[nPiP], gPfPiM[nPiM];
  G4ThreeVector gPfP[nP], gPfPb[nPb];
  G4ThreeVector gPfPiZ[nPiZ];
  G4ThreeVector gPfN[nN], gPfNb[nNb];
  G4ThreeVector gPfKZ[nKZ];
  G4ThreeVector gPfKZS[nKZS], gPfKZL[nKZL];
  G4ThreeVector gPfG[nG];
  G4ThreeVector gPfEP[nEP], gPfEM[nEM];
  G4ThreeVector gPfMuP[nMuP], gPfMuM[nMuM];

  G4int ikp=0, ikm=0;
  G4int ipip=0, ipim=0;
  G4int ip=0, ipb=0;
  G4int ipiz=0;
  G4int in=0, inb=0;
  G4int ikz=0;
  G4int ikzs=0, ikzl=0;
  G4int igamma=0;
  G4int iep=0, iem=0;
  G4int imup=0, imum=0;

  //EvtParticle *q1;
  EvtParticle *q2;
  //q1=parent1->getDaug(0);
  q2=parent2->getDaug(0);

  // do{
  //   EvtId type=q1->getId();
  //   EvtVector4R P4Lab=q1->getP4Lab();

  //   if( type == kp ){
  //     gPfKP[ikp].setX(P4Lab.get(1)*GeV);
  //     gPfKP[ikp].setY(P4Lab.get(2)*GeV);
  //     gPfKP[ikp].setZ(P4Lab.get(3)*GeV);
  //     ikp++; 
  //   }
  //   if( type == km ){
  //     gPfKM[ikm].setX(P4Lab.get(1)*GeV);
  //     gPfKM[ikm].setY(P4Lab.get(2)*GeV);
  //     gPfKM[ikm].setZ(P4Lab.get(3)*GeV);
  //     ikm++;
  //   }
  //   if( type == pip ){
  //     gPfPiP[ipip].setX(P4Lab.get(1)*GeV);
  //     gPfPiP[ipip].setY(P4Lab.get(2)*GeV);
  //     gPfPiP[ipip].setZ(P4Lab.get(3)*GeV);
  //     ipip++;
  //   }
  //   if( type == pim ){
  //     gPfPiM[ipim].setX(P4Lab.get(1)*GeV);
  //     gPfPiM[ipim].setY(P4Lab.get(2)*GeV);
  //     gPfPiM[ipim].setZ(P4Lab.get(3)*GeV);
  //     ipim++;
  //   }
  //   if( type == pr ){
  //     gPfP[ip].setX(P4Lab.get(1)*GeV);
  //     gPfP[ip].setY(P4Lab.get(2)*GeV);
  //     gPfP[ip].setZ(P4Lab.get(3)*GeV);
  //     ip++;
  //   }
  //   if( type == apr ){
  //     gPfPb[ipb].setX(P4Lab.get(1)*GeV);
  //     gPfPb[ipb].setY(P4Lab.get(2)*GeV);
  //     gPfPb[ipb].setZ(P4Lab.get(3)*GeV);
  //     ipb++;
  //   }
  //   if( type == pi0 ){
  //     gPfPiZ[ipiz].setX(P4Lab.get(1)*GeV);
  //     gPfPiZ[ipiz].setY(P4Lab.get(2)*GeV);
  //     gPfPiZ[ipiz].setZ(P4Lab.get(3)*GeV);
  //     ipiz++;
  //   }
  //   if( type == ne ){
  //     gPfN[in].setX(P4Lab.get(1)*GeV);
  //     gPfN[in].setY(P4Lab.get(2)*GeV);
  //     gPfN[in].setZ(P4Lab.get(3)*GeV);
  //     in++;
  //   }
  //   if( type == ane ){
  //     gPfNb[inb].setX(P4Lab.get(1)*GeV);
  //     gPfNb[inb].setY(P4Lab.get(2)*GeV);
  //     gPfNb[inb].setZ(P4Lab.get(3)*GeV);
  //     inb++;
  //   }
  //   if( type == kz ){
  //     gPfKZ[ikz].setX(P4Lab.get(1)*GeV);
  //     gPfKZ[ikz].setY(P4Lab.get(2)*GeV);
  //     gPfKZ[ikz].setZ(P4Lab.get(3)*GeV);
  //     ikz++;
  //   }
  //   if( type == ks ){
  //     gPfKZS[ikzs].setX(P4Lab.get(1)*GeV);
  //     gPfKZS[ikzs].setY(P4Lab.get(2)*GeV);
  //     gPfKZS[ikzs].setZ(P4Lab.get(3)*GeV);
  //     ikzs++;
  //   }
  //   if( type == kl ){
  //     gPfKZL[ikzl].setX(P4Lab.get(1)*GeV);
  //     gPfKZL[ikzl].setY(P4Lab.get(2)*GeV);
  //     gPfKZL[ikzl].setZ(P4Lab.get(3)*GeV);
  //     ikzl++;
  //   }
  //   if( type == gm ){
  //     gPfG[igamma].setX(P4Lab.get(1)*GeV);
  //     gPfG[igamma].setY(P4Lab.get(2)*GeV);
  //     gPfG[igamma].setZ(P4Lab.get(3)*GeV);
  //     igamma++;
  //   }
  //   if( type == ep ){
  //     gPfEP[iep].setX(P4Lab.get(1)*GeV);
  //     gPfEP[iep].setY(P4Lab.get(2)*GeV);
  //     gPfEP[iep].setZ(P4Lab.get(3)*GeV);
  //     iep++;
  //   }
  //   if( type == em ){
  //     gPfEM[iem].setX(P4Lab.get(1)*GeV);
  //     gPfEM[iem].setY(P4Lab.get(2)*GeV);
  //     gPfEM[iem].setZ(P4Lab.get(3)*GeV);
  //     iem++;
  //   }
  //   if( type == mup ){
  //     gPfMuP[imup].setX(P4Lab.get(1)*GeV);
  //     gPfMuP[imup].setY(P4Lab.get(2)*GeV);
  //     gPfMuP[imup].setZ(P4Lab.get(3)*GeV);
  //     imup++;
  //   }
  //   if( type == mum ){
  //     gPfMuM[imum].setX(P4Lab.get(1)*GeV);
  //     gPfMuM[imum].setY(P4Lab.get(2)*GeV);
  //     gPfMuM[imum].setZ(P4Lab.get(3)*GeV);
  //     imum++;
  //   }

  //   q1=q1->nextIter();
  // }while(q1!=0);

  do{
    EvtId type=q2->getId();
    EvtVector4R P4Lab=q2->getP4Lab();

    if( type == kp ){
      gPfKP[ikp].setX(P4Lab.get(1)*GeV);
      gPfKP[ikp].setY(P4Lab.get(2)*GeV);
      gPfKP[ikp].setZ(P4Lab.get(3)*GeV);
      ikp++; 
    }
    if( type == km ){
      gPfKM[ikm].setX(P4Lab.get(1)*GeV);
      gPfKM[ikm].setY(P4Lab.get(2)*GeV);
      gPfKM[ikm].setZ(P4Lab.get(3)*GeV);
      ikm++;
    }
    if( type == pip ){
      gPfPiP[ipip].setX(P4Lab.get(1)*GeV);
      gPfPiP[ipip].setY(P4Lab.get(2)*GeV);
      gPfPiP[ipip].setZ(P4Lab.get(3)*GeV);
      ipip++;
    }
    if( type == pim ){
      gPfPiM[ipim].setX(P4Lab.get(1)*GeV);
      gPfPiM[ipim].setY(P4Lab.get(2)*GeV);
      gPfPiM[ipim].setZ(P4Lab.get(3)*GeV);
      ipim++;
    }
    if( type == pr ){
      gPfP[ip].setX(P4Lab.get(1)*GeV);
      gPfP[ip].setY(P4Lab.get(2)*GeV);
      gPfP[ip].setZ(P4Lab.get(3)*GeV);
      ip++;
    }
    if( type == apr ){
      gPfPb[ipb].setX(P4Lab.get(1)*GeV);
      gPfPb[ipb].setY(P4Lab.get(2)*GeV);
      gPfPb[ipb].setZ(P4Lab.get(3)*GeV);
      ipb++;
    }
    if( type == pi0 ){
      gPfPiZ[ipiz].setX(P4Lab.get(1)*GeV);
      gPfPiZ[ipiz].setY(P4Lab.get(2)*GeV);
      gPfPiZ[ipiz].setZ(P4Lab.get(3)*GeV);
      ipiz++;
    }
    if( type == ne ){
      gPfN[in].setX(P4Lab.get(1)*GeV);
      gPfN[in].setY(P4Lab.get(2)*GeV);
      gPfN[in].setZ(P4Lab.get(3)*GeV);
      in++;
    }
    if( type == ane ){
      gPfNb[inb].setX(P4Lab.get(1)*GeV);
      gPfNb[inb].setY(P4Lab.get(2)*GeV);
      gPfNb[inb].setZ(P4Lab.get(3)*GeV);
      inb++;
    }
    if( type == kz ){
      gPfKZ[ikz].setX(P4Lab.get(1)*GeV);
      gPfKZ[ikz].setY(P4Lab.get(2)*GeV);
      gPfKZ[ikz].setZ(P4Lab.get(3)*GeV);
      ikz++;
    }
    if( type == ks ){
      gPfKZS[ikzs].setX(P4Lab.get(1)*GeV);
      gPfKZS[ikzs].setY(P4Lab.get(2)*GeV);
      gPfKZS[ikzs].setZ(P4Lab.get(3)*GeV);
      ikzs++;
    }
    if( type == kl ){
      gPfKZL[ikzl].setX(P4Lab.get(1)*GeV);
      gPfKZL[ikzl].setY(P4Lab.get(2)*GeV);
      gPfKZL[ikzl].setZ(P4Lab.get(3)*GeV);
      ikzl++;
    }
    if( type == gm ){
      gPfG[igamma].setX(P4Lab.get(1)*GeV);
      gPfG[igamma].setY(P4Lab.get(2)*GeV);
      gPfG[igamma].setZ(P4Lab.get(3)*GeV);
      igamma++;
    }
    if( type == ep ){
      gPfEP[iep].setX(P4Lab.get(1)*GeV);
      gPfEP[iep].setY(P4Lab.get(2)*GeV);
      gPfEP[iep].setZ(P4Lab.get(3)*GeV);
      iep++;
    }
    if( type == em ){
      gPfEM[iem].setX(P4Lab.get(1)*GeV);
      gPfEM[iem].setY(P4Lab.get(2)*GeV);
      gPfEM[iem].setZ(P4Lab.get(3)*GeV);
      iem++;
    }
    if( type == mup ){
      gPfMuP[imup].setX(P4Lab.get(1)*GeV);
      gPfMuP[imup].setY(P4Lab.get(2)*GeV);
      gPfMuP[imup].setZ(P4Lab.get(3)*GeV);
      imup++;
    }
    if( type == mum ){
      gPfMuM[imum].setX(P4Lab.get(1)*GeV);
      gPfMuM[imum].setY(P4Lab.get(2)*GeV);
      gPfMuM[imum].setZ(P4Lab.get(3)*GeV);
      imum++;
    }

    q2=q2->nextIter();
  }while(q2!=0);

  // G4cout<< "**********************" << G4endl;
  // G4cout<< ikp << " " << ikm << G4endl;
  // for( G4int i=0; i<ikp; i++) G4cout<<  gPfKP[i] << G4endl;
  // for( G4int i=0; i<ikm; i++) G4cout<<  gPfKM[i] << G4endl;
  // G4cout<< ipip << " " << ipim << G4endl;
  // for( G4int i=0; i<ipip; i++) G4cout<<  gPfPiP[i] << G4endl;
  // for( G4int i=0; i<ipim; i++) G4cout<<  gPfPiM[i] << G4endl;
  // G4cout<< ip << " " << ipb <<G4endl;
  // for( G4int i=0; i<ip; i++)  G4cout<<  gPfP[i] << G4endl;
  // for( G4int i=0; i<ipb; i++) G4cout<<  gPfPb[i] << G4endl;

  //SetGuns
  for( G4int i=0; i<ikp; i++){
    pGen->gunKP1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKP1_->SetParticleMomentum( gPfKP[i] );
    pGen->gunKP1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikm; i++){
    pGen->gunKM1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKM1_->SetParticleMomentum( gPfKM[i] );
    pGen->gunKM1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipip; i++){
    pGen->gunPiP1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiP1_->SetParticleMomentum( gPfPiP[i] );
    pGen->gunPiP1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipim; i++){
    pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM1_->SetParticleMomentum( gPfPiM[i] );
    pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ip; i++){
    pGen->gunProton_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunProton_->SetParticleMomentum( gPfP[i] );
    pGen->gunProton_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipb; i++){
    pGen->gunProtonb_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunProtonb_->SetParticleMomentum( gPfPb[i] );
    pGen->gunProtonb_->GeneratePrimaryVertex( anEvent );
  }
  // for( G4int i=0; i<ipiz; i++){
  //   pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunPi0_->SetParticleMomentum( gPfPiZ[i] );
  //   pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<in; i++){
  //   pGen->gunNeutron_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunNeutron_->SetParticleMomentum( gPfN[i] );
  //   pGen->gunNeutron_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<inb; i++){
  //   pGen->gunNeutronb_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunNeutronb_->SetParticleMomentum( gPfNb[i] );
  //   pGen->gunNeutronb_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikz; i++){
  //   pGen->gunK0_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0_->SetParticleMomentum( gPfKZ[i] );
  //   pGen->gunK0_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikzs; i++){
  //   pGen->gunK0S_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0S_->SetParticleMomentum( gPfKZS[i] );
  //   pGen->gunK0S_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikzl; i++){
  //   pGen->gunK0L_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0L_->SetParticleMomentum( gPfKZL[i] );
  //   pGen->gunK0L_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<igamma; i++){
  //   pGen->gunGamma_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunGamma_->SetParticleMomentum( gPfG[i] );
  //   pGen->gunGamma_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<iep; i++){
  //   pGen->gunEP_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunEP_->SetParticleMomentum( gPfEP[i] );
  //   pGen->gunEP_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<iem; i++){
  //   pGen->gunEM_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunEM_->SetParticleMomentum( gPfEM[i] );
  //   pGen->gunEM_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<imup; i++){
  //   pGen->gunMuP_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunMuP_->SetParticleMomentum( gPfMuP[i] );
  //   pGen->gunMuP_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<imum; i++){
  //   pGen->gunMuM_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunMuM_->SetParticleMomentum( gPfMuM[i] );
  //   pGen->gunMuM_->GeneratePrimaryVertex( anEvent );
  // }

  if( pGen->anaMan_ ){
    pGen->SetPrimaryInformation( LPos, MF1, MF2, 
                                 LPf1, LPf2, LPf1CM, LPf2CM, pb*LBeamDir, 0. );
  } 
  
  //parent1->deleteTree();
  parent2->deleteTree();

  //XiKStarMesonFromKMinusProton1
  //K- p -> Xi*- K+
}


////////////////////
//K- p -> Xi*- K*+
void PGReaction5::XiKStarMesonFromKMinusProton2( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  // 

  BeamParam *BP = pGen->BP_;

  G4double Mi1=G4KaonMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  //K*+
  G4double Mf1=0.89166*GeV;//K*+ mass
  G4double Wid = 0.05080*GeV;
  static EvtId KS1 = EvtPDL::getId(std::string("K*+"));  
  G4double MF1=BreitWigner( Mf1, Wid );

  //Hyperons
  G4double Mass=BP->RecMass;
  G4double Width=BP->RecWidth;

  static EvtId Xi1;
  if( 1.320*GeV<Mass && Mass<1.330*GeV ) Xi1 = EvtPDL::getId(std::string("Xi-")); 
  else if( 1.530*GeV<Mass && Mass<1.540*GeV ) Xi1 = EvtPDL::getId(std::string("Xi*-")); 
  else if( 1.815*GeV<Mass && Mass<1.825*GeV ) Xi1 = EvtPDL::getId(std::string("Xi(1820)-")); 
  // else if( 1.595*GeV<Mass && Mass<1.605*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1600)0")); 
  // else if( 1.665*GeV<Mass && Mass<1.675*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1670)0")); 
  // else if( 1.685*GeV<Mass && Mass<1.695*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1690)0")); 
  // else if( 1.795*GeV<Mass && Mass<1.805*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1800)0")); 
  // else if( 1.805*GeV<Mass && Mass<1.815*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1810)0")); 
  // else if( 1.815*GeV<Mass && Mass<1.825*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1820)0")); 
  // else if( 1.825*GeV<Mass && Mass<1.835*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1830)0")); 
  // else if( 1.885*GeV<Mass && Mass<1.895*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1890)0")); 
  else{
    Mass = 1.32171*GeV;
    Width = 0.0;
    Xi1 = EvtPDL::getId(std::string("Xi-")); 
  }

  G4double MF2=BreitWigner( Mass, Width );

  //Beam
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4double Ei=sqrt( Mi1*Mi1+Mi2*Mi2+2.*Mi2*sqrt(Mi1*Mi1+pb*pb) );
 
  G4double MKmin=(0.494+0.140)*GeV;//K + pi
  G4double MYmin=(1.322+0.140)*GeV;//Xi + pi

  if( Mass<1.322*GeV ){
    while(1){
      MF1=BreitWigner( Mf1, Wid );
      if(MF1>MKmin ) break;
    }
  }
  else{
    while(1){
      MF1=BreitWigner( Mf1, Wid );
      MF2=BreitWigner( Mass, Width );
      if(MF1>MKmin && MF2>MYmin) break;
    }
  }

  //Reaction
  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;
  AGUniform gen1( BP->cost1, BP->cost2 ); 
  //AGForward2 gen1( BP->cost1, BP->cost2 );

  G4bool status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
				    LPf1, LPf2, LPf1CM, LPf2CM, gen1 );
  G4int n=0;
  while(!status1){
    if(++n>MaxTry){
      G4Exception("XiKStarMesonFromKMinusProton1: under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    while( !(Ei>MF1+MF2) ){
      if( Mass<1.322*GeV ){
	while(1){
	  MF1=BreitWigner( Mf1, Wid );
	  if(MF1>MKmin ) break;
	}
      }
      else{
	while(1){
	  MF1=BreitWigner( Mf1, Wid );
	  MF2=BreitWigner( Mass, Width );
	  if(MF1>MKmin && MF2>MYmin) break;
	}
      }
    }
    status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, gen1 );
  }

  G4ThreeVector gPos(LPos);

  //Beam
  pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM1_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );

  ///////K*+ && Xi* decay by EvtGen
  EvtParticle* parent1(0);
  EvtParticle* parent2(0);
  
  G4LorentzVector LvKS;
  LvKS.setVect(LPf1/GeV);
  LvKS.setE(sqrt(MF1*MF1+LPf1.mag2())/GeV);
  
  G4LorentzVector LvXi;
  LvXi.setVect(LPf2/GeV);
  LvXi.setE(sqrt(MF2*MF2+LPf2.mag2())/GeV);
  
  EvtVector4R pInit1( LvKS.e(),
		      LvKS.vect().x(), LvKS.vect().y(), LvKS.vect().z() );
  EvtVector4R pInit2( LvXi.e(),
		      LvXi.vect().x(), LvXi.vect().y(), LvXi.vect().z() );

  parent1 = EvtParticleFactory::particleFactory(KS1, pInit1); 
  evtGenerator->generateDecay(parent1);
  parent2 = EvtParticleFactory::particleFactory(Xi1, pInit2); 
  evtGenerator->generateDecay(parent2);

  //std::cout<< "******" << std::endl;

  static EvtId kp=EvtPDL::getId(std::string("K+"));
  static EvtId km=EvtPDL::getId(std::string("K-"));
  static EvtId pip=EvtPDL::getId(std::string("pi+"));
  static EvtId pim=EvtPDL::getId(std::string("pi-"));
  static EvtId pr=EvtPDL::getId(std::string("p+"));
  static EvtId apr=EvtPDL::getId(std::string("anti-p-"));
  static EvtId pi0=EvtPDL::getId(std::string("pi0"));
  static EvtId ne=EvtPDL::getId(std::string("n0"));
  static EvtId ane=EvtPDL::getId(std::string("anti-n0"));
  static EvtId kz=EvtPDL::getId(std::string("K0"));
  static EvtId ks=EvtPDL::getId(std::string("K_S0"));
  static EvtId kl=EvtPDL::getId(std::string("K_L0"));
  static EvtId gm=EvtPDL::getId(std::string("gamma"));
  static EvtId ep=EvtPDL::getId(std::string("e+"));
  static EvtId em=EvtPDL::getId(std::string("e-"));
  static EvtId mup=EvtPDL::getId(std::string("mu+"));
  static EvtId mum=EvtPDL::getId(std::string("mu-"));

  G4int nkp=0, nkm=0;
  G4int npip=0, npim=0;
  G4int nproton=0, nprotonb=0;
  G4int npiz=0;
  G4int nneutron=0, nneutronb=0;
  G4int nkz=0, nkzs=0, nkzl=0;
  G4int ngamma=0;
  G4int nep=0, nem=0;
  G4int nmup=0, nmum=0;

  EvtParticle *p1;
  EvtParticle *p2;
  p1=parent1->getDaug(0);
  p2=parent2->getDaug(0);

  do{
    EvtId type=p1->getId();
    
    if (type==kp)  nkp++;
    if (type==km)  nkm++;
    if (type==pip) npip++;
    if (type==pim) npim++;
    if (type==pr)  nproton++;
    if (type==apr) nprotonb++;
    if (type==pi0) npiz++;
    if (type==ne)  nneutron++;
    if (type==ane) nneutronb++;
    if (type==kz)  nkz++;
    if (type==ks)  nkzs++;
    if (type==kl)  nkzl++;
    if (type==gm)  ngamma++;
    if (type==ep)  nep++;
    if (type==em)  nem++;
    if (type==mup) nmup++;
    if (type==mum) nmum++;
    
    p1=p1->nextIter();
  }while(p1!=0);

  do{
    EvtId type=p2->getId();
    
    if (type==kp)  nkp++;
    if (type==km)  nkm++;
    if (type==pip) npip++;
    if (type==pim) npim++;
    if (type==pr)  nproton++;
    if (type==apr) nprotonb++;
    if (type==pi0) npiz++;
    if (type==ne)  nneutron++;
    if (type==ane) nneutronb++;
    if (type==kz)  nkz++;
    if (type==ks)  nkzs++;
    if (type==kl)  nkzl++;
    if (type==gm)  ngamma++;
    if (type==ep)  nep++;
    if (type==em)  nem++;
    if (type==mup) nmup++;
    if (type==mum) nmum++;
    
    p2=p2->nextIter();
  }while(p2!=0);

  const G4int nKP=nkp;
  const G4int nKM=nkm;
  const G4int nPiP=npip;
  const G4int nPiM=npim;
  const G4int nP=nproton;
  const G4int nPb=nprotonb;
  const G4int nPiZ=npiz;
  const G4int nN=nneutron;
  const G4int nNb=nneutronb;
  const G4int nKZ=nkz;
  const G4int nKZS=nkzs;
  const G4int nKZL=nkzl;
  const G4int nG=ngamma;
  const G4int nEP=nep;
  const G4int nEM=nem;
  const G4int nMuP=nmup;
  const G4int nMuM=nmum;

  G4ThreeVector gPfKP[nKP], gPfKM[nKM];
  G4ThreeVector gPfPiP[nPiP], gPfPiM[nPiM];
  G4ThreeVector gPfP[nP], gPfPb[nPb];
  G4ThreeVector gPfPiZ[nPiZ];
  G4ThreeVector gPfN[nN], gPfNb[nNb];
  G4ThreeVector gPfKZ[nKZ];
  G4ThreeVector gPfKZS[nKZS], gPfKZL[nKZL];
  G4ThreeVector gPfG[nG];
  G4ThreeVector gPfEP[nEP], gPfEM[nEM];
  G4ThreeVector gPfMuP[nMuP], gPfMuM[nMuM];

  G4int ikp=0, ikm=0;
  G4int ipip=0, ipim=0;
  G4int ip=0, ipb=0;
  G4int ipiz=0;
  G4int in=0, inb=0;
  G4int ikz=0;
  G4int ikzs=0, ikzl=0;
  G4int igamma=0;
  G4int iep=0, iem=0;
  G4int imup=0, imum=0;

  EvtParticle *q1;
  EvtParticle *q2;
  q1=parent1->getDaug(0);
  q2=parent2->getDaug(0);

  do{
    EvtId type=q1->getId();
    EvtVector4R P4Lab=q1->getP4Lab();

    if( type == kp ){
      gPfKP[ikp].setX(P4Lab.get(1)*GeV);
      gPfKP[ikp].setY(P4Lab.get(2)*GeV);
      gPfKP[ikp].setZ(P4Lab.get(3)*GeV);
      ikp++; 
    }
    if( type == km ){
      gPfKM[ikm].setX(P4Lab.get(1)*GeV);
      gPfKM[ikm].setY(P4Lab.get(2)*GeV);
      gPfKM[ikm].setZ(P4Lab.get(3)*GeV);
      ikm++;
    }
    if( type == pip ){
      gPfPiP[ipip].setX(P4Lab.get(1)*GeV);
      gPfPiP[ipip].setY(P4Lab.get(2)*GeV);
      gPfPiP[ipip].setZ(P4Lab.get(3)*GeV);
      ipip++;
    }
    if( type == pim ){
      gPfPiM[ipim].setX(P4Lab.get(1)*GeV);
      gPfPiM[ipim].setY(P4Lab.get(2)*GeV);
      gPfPiM[ipim].setZ(P4Lab.get(3)*GeV);
      ipim++;
    }
    if( type == pr ){
      gPfP[ip].setX(P4Lab.get(1)*GeV);
      gPfP[ip].setY(P4Lab.get(2)*GeV);
      gPfP[ip].setZ(P4Lab.get(3)*GeV);
      ip++;
    }
    if( type == apr ){
      gPfPb[ipb].setX(P4Lab.get(1)*GeV);
      gPfPb[ipb].setY(P4Lab.get(2)*GeV);
      gPfPb[ipb].setZ(P4Lab.get(3)*GeV);
      ipb++;
    }
    if( type == pi0 ){
      gPfPiZ[ipiz].setX(P4Lab.get(1)*GeV);
      gPfPiZ[ipiz].setY(P4Lab.get(2)*GeV);
      gPfPiZ[ipiz].setZ(P4Lab.get(3)*GeV);
      ipiz++;
    }
    if( type == ne ){
      gPfN[in].setX(P4Lab.get(1)*GeV);
      gPfN[in].setY(P4Lab.get(2)*GeV);
      gPfN[in].setZ(P4Lab.get(3)*GeV);
      in++;
    }
    if( type == ane ){
      gPfNb[inb].setX(P4Lab.get(1)*GeV);
      gPfNb[inb].setY(P4Lab.get(2)*GeV);
      gPfNb[inb].setZ(P4Lab.get(3)*GeV);
      inb++;
    }
    if( type == kz ){
      gPfKZ[ikz].setX(P4Lab.get(1)*GeV);
      gPfKZ[ikz].setY(P4Lab.get(2)*GeV);
      gPfKZ[ikz].setZ(P4Lab.get(3)*GeV);
      ikz++;
    }
    if( type == ks ){
      gPfKZS[ikzs].setX(P4Lab.get(1)*GeV);
      gPfKZS[ikzs].setY(P4Lab.get(2)*GeV);
      gPfKZS[ikzs].setZ(P4Lab.get(3)*GeV);
      ikzs++;
    }
    if( type == kl ){
      gPfKZL[ikzl].setX(P4Lab.get(1)*GeV);
      gPfKZL[ikzl].setY(P4Lab.get(2)*GeV);
      gPfKZL[ikzl].setZ(P4Lab.get(3)*GeV);
      ikzl++;
    }
    if( type == gm ){
      gPfG[igamma].setX(P4Lab.get(1)*GeV);
      gPfG[igamma].setY(P4Lab.get(2)*GeV);
      gPfG[igamma].setZ(P4Lab.get(3)*GeV);
      igamma++;
    }
    if( type == ep ){
      gPfEP[iep].setX(P4Lab.get(1)*GeV);
      gPfEP[iep].setY(P4Lab.get(2)*GeV);
      gPfEP[iep].setZ(P4Lab.get(3)*GeV);
      iep++;
    }
    if( type == em ){
      gPfEM[iem].setX(P4Lab.get(1)*GeV);
      gPfEM[iem].setY(P4Lab.get(2)*GeV);
      gPfEM[iem].setZ(P4Lab.get(3)*GeV);
      iem++;
    }
    if( type == mup ){
      gPfMuP[imup].setX(P4Lab.get(1)*GeV);
      gPfMuP[imup].setY(P4Lab.get(2)*GeV);
      gPfMuP[imup].setZ(P4Lab.get(3)*GeV);
      imup++;
    }
    if( type == mum ){
      gPfMuM[imum].setX(P4Lab.get(1)*GeV);
      gPfMuM[imum].setY(P4Lab.get(2)*GeV);
      gPfMuM[imum].setZ(P4Lab.get(3)*GeV);
      imum++;
    }

    q1=q1->nextIter();
  }while(q1!=0);

  do{
    EvtId type=q2->getId();
    EvtVector4R P4Lab=q2->getP4Lab();

    if( type == kp ){
      gPfKP[ikp].setX(P4Lab.get(1)*GeV);
      gPfKP[ikp].setY(P4Lab.get(2)*GeV);
      gPfKP[ikp].setZ(P4Lab.get(3)*GeV);
      ikp++; 
    }
    if( type == km ){
      gPfKM[ikm].setX(P4Lab.get(1)*GeV);
      gPfKM[ikm].setY(P4Lab.get(2)*GeV);
      gPfKM[ikm].setZ(P4Lab.get(3)*GeV);
      ikm++;
    }
    if( type == pip ){
      gPfPiP[ipip].setX(P4Lab.get(1)*GeV);
      gPfPiP[ipip].setY(P4Lab.get(2)*GeV);
      gPfPiP[ipip].setZ(P4Lab.get(3)*GeV);
      ipip++;
    }
    if( type == pim ){
      gPfPiM[ipim].setX(P4Lab.get(1)*GeV);
      gPfPiM[ipim].setY(P4Lab.get(2)*GeV);
      gPfPiM[ipim].setZ(P4Lab.get(3)*GeV);
      ipim++;
    }
    if( type == pr ){
      gPfP[ip].setX(P4Lab.get(1)*GeV);
      gPfP[ip].setY(P4Lab.get(2)*GeV);
      gPfP[ip].setZ(P4Lab.get(3)*GeV);
      ip++;
    }
    if( type == apr ){
      gPfPb[ipb].setX(P4Lab.get(1)*GeV);
      gPfPb[ipb].setY(P4Lab.get(2)*GeV);
      gPfPb[ipb].setZ(P4Lab.get(3)*GeV);
      ipb++;
    }
    if( type == pi0 ){
      gPfPiZ[ipiz].setX(P4Lab.get(1)*GeV);
      gPfPiZ[ipiz].setY(P4Lab.get(2)*GeV);
      gPfPiZ[ipiz].setZ(P4Lab.get(3)*GeV);
      ipiz++;
    }
    if( type == ne ){
      gPfN[in].setX(P4Lab.get(1)*GeV);
      gPfN[in].setY(P4Lab.get(2)*GeV);
      gPfN[in].setZ(P4Lab.get(3)*GeV);
      in++;
    }
    if( type == ane ){
      gPfNb[inb].setX(P4Lab.get(1)*GeV);
      gPfNb[inb].setY(P4Lab.get(2)*GeV);
      gPfNb[inb].setZ(P4Lab.get(3)*GeV);
      inb++;
    }
    if( type == kz ){
      gPfKZ[ikz].setX(P4Lab.get(1)*GeV);
      gPfKZ[ikz].setY(P4Lab.get(2)*GeV);
      gPfKZ[ikz].setZ(P4Lab.get(3)*GeV);
      ikz++;
    }
    if( type == ks ){
      gPfKZS[ikzs].setX(P4Lab.get(1)*GeV);
      gPfKZS[ikzs].setY(P4Lab.get(2)*GeV);
      gPfKZS[ikzs].setZ(P4Lab.get(3)*GeV);
      ikzs++;
    }
    if( type == kl ){
      gPfKZL[ikzl].setX(P4Lab.get(1)*GeV);
      gPfKZL[ikzl].setY(P4Lab.get(2)*GeV);
      gPfKZL[ikzl].setZ(P4Lab.get(3)*GeV);
      ikzl++;
    }
    if( type == gm ){
      gPfG[igamma].setX(P4Lab.get(1)*GeV);
      gPfG[igamma].setY(P4Lab.get(2)*GeV);
      gPfG[igamma].setZ(P4Lab.get(3)*GeV);
      igamma++;
    }
    if( type == ep ){
      gPfEP[iep].setX(P4Lab.get(1)*GeV);
      gPfEP[iep].setY(P4Lab.get(2)*GeV);
      gPfEP[iep].setZ(P4Lab.get(3)*GeV);
      iep++;
    }
    if( type == em ){
      gPfEM[iem].setX(P4Lab.get(1)*GeV);
      gPfEM[iem].setY(P4Lab.get(2)*GeV);
      gPfEM[iem].setZ(P4Lab.get(3)*GeV);
      iem++;
    }
    if( type == mup ){
      gPfMuP[imup].setX(P4Lab.get(1)*GeV);
      gPfMuP[imup].setY(P4Lab.get(2)*GeV);
      gPfMuP[imup].setZ(P4Lab.get(3)*GeV);
      imup++;
    }
    if( type == mum ){
      gPfMuM[imum].setX(P4Lab.get(1)*GeV);
      gPfMuM[imum].setY(P4Lab.get(2)*GeV);
      gPfMuM[imum].setZ(P4Lab.get(3)*GeV);
      imum++;
    }

    q2=q2->nextIter();
  }while(q2!=0);

  // G4cout<< "**********************" << G4endl;
  // G4cout<< ikp << " " << ikm << G4endl;
  // for( G4int i=0; i<ikp; i++) G4cout<<  gPfKP[i] << G4endl;
  // for( G4int i=0; i<ikm; i++) G4cout<<  gPfKM[i] << G4endl;
  // G4cout<< ipip << " " << ipim << G4endl;
  // for( G4int i=0; i<ipip; i++) G4cout<<  gPfPiP[i] << G4endl;
  // for( G4int i=0; i<ipim; i++) G4cout<<  gPfPiM[i] << G4endl;
  // G4cout<< ip << " " << ipb <<G4endl;
  // for( G4int i=0; i<ip; i++)  G4cout<<  gPfP[i] << G4endl;
  // for( G4int i=0; i<ipb; i++) G4cout<<  gPfPb[i] << G4endl;

  //SetGuns
  for( G4int i=0; i<ikp; i++){
    pGen->gunKP1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKP1_->SetParticleMomentum( gPfKP[i] );
    pGen->gunKP1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikm; i++){
    pGen->gunKM1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKM1_->SetParticleMomentum( gPfKM[i] );
    pGen->gunKM1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipip; i++){
    pGen->gunPiP1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiP1_->SetParticleMomentum( gPfPiP[i] );
    pGen->gunPiP1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipim; i++){
    pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM1_->SetParticleMomentum( gPfPiM[i] );
    pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ip; i++){
    pGen->gunProton_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunProton_->SetParticleMomentum( gPfP[i] );
    pGen->gunProton_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipb; i++){
    pGen->gunProtonb_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunProtonb_->SetParticleMomentum( gPfPb[i] );
    pGen->gunProtonb_->GeneratePrimaryVertex( anEvent );
  }
  // for( G4int i=0; i<ipiz; i++){
  //   pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunPi0_->SetParticleMomentum( gPfPiZ[i] );
  //   pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<in; i++){
  //   pGen->gunNeutron_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunNeutron_->SetParticleMomentum( gPfN[i] );
  //   pGen->gunNeutron_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<inb; i++){
  //   pGen->gunNeutronb_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunNeutronb_->SetParticleMomentum( gPfNb[i] );
  //   pGen->gunNeutronb_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikz; i++){
  //   pGen->gunK0_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0_->SetParticleMomentum( gPfKZ[i] );
  //   pGen->gunK0_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikzs; i++){
  //   pGen->gunK0S_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0S_->SetParticleMomentum( gPfKZS[i] );
  //   pGen->gunK0S_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikzl; i++){
  //   pGen->gunK0L_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0L_->SetParticleMomentum( gPfKZL[i] );
  //   pGen->gunK0L_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<igamma; i++){
  //   pGen->gunGamma_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunGamma_->SetParticleMomentum( gPfG[i] );
  //   pGen->gunGamma_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<iep; i++){
  //   pGen->gunEP_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunEP_->SetParticleMomentum( gPfEP[i] );
  //   pGen->gunEP_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<iem; i++){
  //   pGen->gunEM_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunEM_->SetParticleMomentum( gPfEM[i] );
  //   pGen->gunEM_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<imup; i++){
  //   pGen->gunMuP_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunMuP_->SetParticleMomentum( gPfMuP[i] );
  //   pGen->gunMuP_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<imum; i++){
  //   pGen->gunMuM_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunMuM_->SetParticleMomentum( gPfMuM[i] );
  //   pGen->gunMuM_->GeneratePrimaryVertex( anEvent );
  // }

  if( pGen->anaMan_ ){
    pGen->SetPrimaryInformation( LPos, MF1, MF2, 
                                 LPf1, LPf2, LPf1CM, LPf2CM, pb*LBeamDir, 0. );
  } 
  
  parent1->deleteTree();
  parent2->deleteTree();

  //XiKStarMesonFromKMinusProton2
  //K- p -> Xi*- K*+
}


////////////////////
//K- p -> Xi*0 K*0
void PGReaction5::XiKStarMesonFromKMinusProton3( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  // 

  BeamParam *BP = pGen->BP_;

  G4double Mi1=G4KaonMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  //K*0 
  G4double Mf1=0.8960*GeV;//K*0 mass
  G4double Wid = 0.05030*GeV;
  static EvtId KS1 = EvtPDL::getId(std::string("K*0"));  
  G4double MF1=BreitWigner( Mf1, Wid );

  //Hyperons
  G4double Mass=BP->RecMass;
  G4double Width=BP->RecWidth;

  static EvtId Xi1;
  if( 1.320*GeV<Mass && Mass<1.330*GeV ) Xi1 = EvtPDL::getId(std::string("Xi-")); 
  else if( 1.530*GeV<Mass && Mass<1.540*GeV ) Xi1 = EvtPDL::getId(std::string("Xi*-")); 
  else if( 1.815*GeV<Mass && Mass<1.825*GeV ) Xi1 = EvtPDL::getId(std::string("Xi(1820)-")); 
  // else if( 1.595*GeV<Mass && Mass<1.605*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1600)0")); 
  // else if( 1.665*GeV<Mass && Mass<1.675*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1670)0")); 
  // else if( 1.685*GeV<Mass && Mass<1.695*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1690)0")); 
  // else if( 1.795*GeV<Mass && Mass<1.805*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1800)0")); 
  // else if( 1.805*GeV<Mass && Mass<1.815*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1810)0")); 
  // else if( 1.815*GeV<Mass && Mass<1.825*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1820)0")); 
  // else if( 1.825*GeV<Mass && Mass<1.835*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1830)0")); 
  // else if( 1.885*GeV<Mass && Mass<1.895*GeV ) Xi1 = EvtPDL::getId(std::string("Lambda(1890)0")); 
  else{
    Mass = 1.32171*GeV;
    Width = 0.0;
    Xi1 = EvtPDL::getId(std::string("Xi-")); 
  }

  G4double MF2=BreitWigner( Mass, Width );

  //Beam
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4double Ei=sqrt( Mi1*Mi1+Mi2*Mi2+2.*Mi2*sqrt(Mi1*Mi1+pb*pb) );
 
  G4double MKmin=(0.494+0.140)*GeV;//K + pi
  G4double MYmin=(1.322+0.140)*GeV;//Xi + pi

  if( Mass<1.322*GeV ){
    while(1){
      MF1=BreitWigner( Mf1, Wid );
      if(MF1>MKmin ) break;
    }
  }
  else{
    while(1){
      MF1=BreitWigner( Mf1, Wid );
      MF2=BreitWigner( Mass, Width );
      if(MF1>MKmin && MF2>MYmin) break;
    }
  }

  //Reaction
  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;
  AGUniform gen1( BP->cost1, BP->cost2 ); 
  //AGForward2 gen1( BP->cost1, BP->cost2 );

  G4bool status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
				    LPf1, LPf2, LPf1CM, LPf2CM, gen1 );
  G4int n=0;
  while(!status1){
    if(++n>MaxTry){
      G4Exception("XiKStarMesonFromKMinusProton1: under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    while( !(Ei>MF1+MF2) ){
      if( Mass<1.322*GeV ){
	while(1){
	  MF1=BreitWigner( Mf1, Wid );
	  if(MF1>MKmin ) break;
	}
      }
      else{
	while(1){
	  MF1=BreitWigner( Mf1, Wid );
	  MF2=BreitWigner( Mass, Width );
	  if(MF1>MKmin && MF2>MYmin) break;
	}
      }
    }
    status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, gen1 );
  }

  G4ThreeVector gPos(LPos);

  //Beam
  pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM1_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );

  ///////K*0 && Xi* decay by EvtGen
  EvtParticle* parent1(0);
  EvtParticle* parent2(0);
  
  G4LorentzVector LvKS;
  LvKS.setVect(LPf1/GeV);
  LvKS.setE(sqrt(MF1*MF1+LPf1.mag2())/GeV);
  
  G4LorentzVector LvXi;
  LvXi.setVect(LPf2/GeV);
  LvXi.setE(sqrt(MF2*MF2+LPf2.mag2())/GeV);
  
  EvtVector4R pInit1( LvKS.e(),
		      LvKS.vect().x(), LvKS.vect().y(), LvKS.vect().z() );
  EvtVector4R pInit2( LvXi.e(),
		      LvXi.vect().x(), LvXi.vect().y(), LvXi.vect().z() );

  parent1 = EvtParticleFactory::particleFactory(KS1, pInit1); 
  evtGenerator->generateDecay(parent1);
  parent2 = EvtParticleFactory::particleFactory(Xi1, pInit2); 
  evtGenerator->generateDecay(parent2);

  //std::cout<< "******" << std::endl;

  static EvtId kp=EvtPDL::getId(std::string("K+"));
  static EvtId km=EvtPDL::getId(std::string("K-"));
  static EvtId pip=EvtPDL::getId(std::string("pi+"));
  static EvtId pim=EvtPDL::getId(std::string("pi-"));
  static EvtId pr=EvtPDL::getId(std::string("p+"));
  static EvtId apr=EvtPDL::getId(std::string("anti-p-"));
  static EvtId pi0=EvtPDL::getId(std::string("pi0"));
  static EvtId ne=EvtPDL::getId(std::string("n0"));
  static EvtId ane=EvtPDL::getId(std::string("anti-n0"));
  static EvtId kz=EvtPDL::getId(std::string("K0"));
  static EvtId ks=EvtPDL::getId(std::string("K_S0"));
  static EvtId kl=EvtPDL::getId(std::string("K_L0"));
  static EvtId gm=EvtPDL::getId(std::string("gamma"));
  static EvtId ep=EvtPDL::getId(std::string("e+"));
  static EvtId em=EvtPDL::getId(std::string("e-"));
  static EvtId mup=EvtPDL::getId(std::string("mu+"));
  static EvtId mum=EvtPDL::getId(std::string("mu-"));

  G4int nkp=0, nkm=0;
  G4int npip=0, npim=0;
  G4int nproton=0, nprotonb=0;
  G4int npiz=0;
  G4int nneutron=0, nneutronb=0;
  G4int nkz=0, nkzs=0, nkzl=0;
  G4int ngamma=0;
  G4int nep=0, nem=0;
  G4int nmup=0, nmum=0;

  EvtParticle *p1;
  EvtParticle *p2;
  p1=parent1->getDaug(0);
  p2=parent2->getDaug(0);

  do{
    EvtId type=p1->getId();
    
    if (type==kp)  nkp++;
    if (type==km)  nkm++;
    if (type==pip) npip++;
    if (type==pim) npim++;
    if (type==pr)  nproton++;
    if (type==apr) nprotonb++;
    if (type==pi0) npiz++;
    if (type==ne)  nneutron++;
    if (type==ane) nneutronb++;
    if (type==kz)  nkz++;
    if (type==ks)  nkzs++;
    if (type==kl)  nkzl++;
    if (type==gm)  ngamma++;
    if (type==ep)  nep++;
    if (type==em)  nem++;
    if (type==mup) nmup++;
    if (type==mum) nmum++;
    
    p1=p1->nextIter();
  }while(p1!=0);

  do{
    EvtId type=p2->getId();
    
    if (type==kp)  nkp++;
    if (type==km)  nkm++;
    if (type==pip) npip++;
    if (type==pim) npim++;
    if (type==pr)  nproton++;
    if (type==apr) nprotonb++;
    if (type==pi0) npiz++;
    if (type==ne)  nneutron++;
    if (type==ane) nneutronb++;
    if (type==kz)  nkz++;
    if (type==ks)  nkzs++;
    if (type==kl)  nkzl++;
    if (type==gm)  ngamma++;
    if (type==ep)  nep++;
    if (type==em)  nem++;
    if (type==mup) nmup++;
    if (type==mum) nmum++;
    
    p2=p2->nextIter();
  }while(p2!=0);

  const G4int nKP=nkp;
  const G4int nKM=nkm;
  const G4int nPiP=npip;
  const G4int nPiM=npim;
  const G4int nP=nproton;
  const G4int nPb=nprotonb;
  const G4int nPiZ=npiz;
  const G4int nN=nneutron;
  const G4int nNb=nneutronb;
  const G4int nKZ=nkz;
  const G4int nKZS=nkzs;
  const G4int nKZL=nkzl;
  const G4int nG=ngamma;
  const G4int nEP=nep;
  const G4int nEM=nem;
  const G4int nMuP=nmup;
  const G4int nMuM=nmum;

  G4ThreeVector gPfKP[nKP], gPfKM[nKM];
  G4ThreeVector gPfPiP[nPiP], gPfPiM[nPiM];
  G4ThreeVector gPfP[nP], gPfPb[nPb];
  G4ThreeVector gPfPiZ[nPiZ];
  G4ThreeVector gPfN[nN], gPfNb[nNb];
  G4ThreeVector gPfKZ[nKZ];
  G4ThreeVector gPfKZS[nKZS], gPfKZL[nKZL];
  G4ThreeVector gPfG[nG];
  G4ThreeVector gPfEP[nEP], gPfEM[nEM];
  G4ThreeVector gPfMuP[nMuP], gPfMuM[nMuM];

  G4int ikp=0, ikm=0;
  G4int ipip=0, ipim=0;
  G4int ip=0, ipb=0;
  G4int ipiz=0;
  G4int in=0, inb=0;
  G4int ikz=0;
  G4int ikzs=0, ikzl=0;
  G4int igamma=0;
  G4int iep=0, iem=0;
  G4int imup=0, imum=0;

  EvtParticle *q1;
  EvtParticle *q2;
  q1=parent1->getDaug(0);
  q2=parent2->getDaug(0);

  do{
    EvtId type=q1->getId();
    EvtVector4R P4Lab=q1->getP4Lab();

    if( type == kp ){
      gPfKP[ikp].setX(P4Lab.get(1)*GeV);
      gPfKP[ikp].setY(P4Lab.get(2)*GeV);
      gPfKP[ikp].setZ(P4Lab.get(3)*GeV);
      ikp++; 
    }
    if( type == km ){
      gPfKM[ikm].setX(P4Lab.get(1)*GeV);
      gPfKM[ikm].setY(P4Lab.get(2)*GeV);
      gPfKM[ikm].setZ(P4Lab.get(3)*GeV);
      ikm++;
    }
    if( type == pip ){
      gPfPiP[ipip].setX(P4Lab.get(1)*GeV);
      gPfPiP[ipip].setY(P4Lab.get(2)*GeV);
      gPfPiP[ipip].setZ(P4Lab.get(3)*GeV);
      ipip++;
    }
    if( type == pim ){
      gPfPiM[ipim].setX(P4Lab.get(1)*GeV);
      gPfPiM[ipim].setY(P4Lab.get(2)*GeV);
      gPfPiM[ipim].setZ(P4Lab.get(3)*GeV);
      ipim++;
    }
    if( type == pr ){
      gPfP[ip].setX(P4Lab.get(1)*GeV);
      gPfP[ip].setY(P4Lab.get(2)*GeV);
      gPfP[ip].setZ(P4Lab.get(3)*GeV);
      ip++;
    }
    if( type == apr ){
      gPfPb[ipb].setX(P4Lab.get(1)*GeV);
      gPfPb[ipb].setY(P4Lab.get(2)*GeV);
      gPfPb[ipb].setZ(P4Lab.get(3)*GeV);
      ipb++;
    }
    if( type == pi0 ){
      gPfPiZ[ipiz].setX(P4Lab.get(1)*GeV);
      gPfPiZ[ipiz].setY(P4Lab.get(2)*GeV);
      gPfPiZ[ipiz].setZ(P4Lab.get(3)*GeV);
      ipiz++;
    }
    if( type == ne ){
      gPfN[in].setX(P4Lab.get(1)*GeV);
      gPfN[in].setY(P4Lab.get(2)*GeV);
      gPfN[in].setZ(P4Lab.get(3)*GeV);
      in++;
    }
    if( type == ane ){
      gPfNb[inb].setX(P4Lab.get(1)*GeV);
      gPfNb[inb].setY(P4Lab.get(2)*GeV);
      gPfNb[inb].setZ(P4Lab.get(3)*GeV);
      inb++;
    }
    if( type == kz ){
      gPfKZ[ikz].setX(P4Lab.get(1)*GeV);
      gPfKZ[ikz].setY(P4Lab.get(2)*GeV);
      gPfKZ[ikz].setZ(P4Lab.get(3)*GeV);
      ikz++;
    }
    if( type == ks ){
      gPfKZS[ikzs].setX(P4Lab.get(1)*GeV);
      gPfKZS[ikzs].setY(P4Lab.get(2)*GeV);
      gPfKZS[ikzs].setZ(P4Lab.get(3)*GeV);
      ikzs++;
    }
    if( type == kl ){
      gPfKZL[ikzl].setX(P4Lab.get(1)*GeV);
      gPfKZL[ikzl].setY(P4Lab.get(2)*GeV);
      gPfKZL[ikzl].setZ(P4Lab.get(3)*GeV);
      ikzl++;
    }
    if( type == gm ){
      gPfG[igamma].setX(P4Lab.get(1)*GeV);
      gPfG[igamma].setY(P4Lab.get(2)*GeV);
      gPfG[igamma].setZ(P4Lab.get(3)*GeV);
      igamma++;
    }
    if( type == ep ){
      gPfEP[iep].setX(P4Lab.get(1)*GeV);
      gPfEP[iep].setY(P4Lab.get(2)*GeV);
      gPfEP[iep].setZ(P4Lab.get(3)*GeV);
      iep++;
    }
    if( type == em ){
      gPfEM[iem].setX(P4Lab.get(1)*GeV);
      gPfEM[iem].setY(P4Lab.get(2)*GeV);
      gPfEM[iem].setZ(P4Lab.get(3)*GeV);
      iem++;
    }
    if( type == mup ){
      gPfMuP[imup].setX(P4Lab.get(1)*GeV);
      gPfMuP[imup].setY(P4Lab.get(2)*GeV);
      gPfMuP[imup].setZ(P4Lab.get(3)*GeV);
      imup++;
    }
    if( type == mum ){
      gPfMuM[imum].setX(P4Lab.get(1)*GeV);
      gPfMuM[imum].setY(P4Lab.get(2)*GeV);
      gPfMuM[imum].setZ(P4Lab.get(3)*GeV);
      imum++;
    }

    q1=q1->nextIter();
  }while(q1!=0);

  do{
    EvtId type=q2->getId();
    EvtVector4R P4Lab=q2->getP4Lab();

    if( type == kp ){
      gPfKP[ikp].setX(P4Lab.get(1)*GeV);
      gPfKP[ikp].setY(P4Lab.get(2)*GeV);
      gPfKP[ikp].setZ(P4Lab.get(3)*GeV);
      ikp++; 
    }
    if( type == km ){
      gPfKM[ikm].setX(P4Lab.get(1)*GeV);
      gPfKM[ikm].setY(P4Lab.get(2)*GeV);
      gPfKM[ikm].setZ(P4Lab.get(3)*GeV);
      ikm++;
    }
    if( type == pip ){
      gPfPiP[ipip].setX(P4Lab.get(1)*GeV);
      gPfPiP[ipip].setY(P4Lab.get(2)*GeV);
      gPfPiP[ipip].setZ(P4Lab.get(3)*GeV);
      ipip++;
    }
    if( type == pim ){
      gPfPiM[ipim].setX(P4Lab.get(1)*GeV);
      gPfPiM[ipim].setY(P4Lab.get(2)*GeV);
      gPfPiM[ipim].setZ(P4Lab.get(3)*GeV);
      ipim++;
    }
    if( type == pr ){
      gPfP[ip].setX(P4Lab.get(1)*GeV);
      gPfP[ip].setY(P4Lab.get(2)*GeV);
      gPfP[ip].setZ(P4Lab.get(3)*GeV);
      ip++;
    }
    if( type == apr ){
      gPfPb[ipb].setX(P4Lab.get(1)*GeV);
      gPfPb[ipb].setY(P4Lab.get(2)*GeV);
      gPfPb[ipb].setZ(P4Lab.get(3)*GeV);
      ipb++;
    }
    if( type == pi0 ){
      gPfPiZ[ipiz].setX(P4Lab.get(1)*GeV);
      gPfPiZ[ipiz].setY(P4Lab.get(2)*GeV);
      gPfPiZ[ipiz].setZ(P4Lab.get(3)*GeV);
      ipiz++;
    }
    if( type == ne ){
      gPfN[in].setX(P4Lab.get(1)*GeV);
      gPfN[in].setY(P4Lab.get(2)*GeV);
      gPfN[in].setZ(P4Lab.get(3)*GeV);
      in++;
    }
    if( type == ane ){
      gPfNb[inb].setX(P4Lab.get(1)*GeV);
      gPfNb[inb].setY(P4Lab.get(2)*GeV);
      gPfNb[inb].setZ(P4Lab.get(3)*GeV);
      inb++;
    }
    if( type == kz ){
      gPfKZ[ikz].setX(P4Lab.get(1)*GeV);
      gPfKZ[ikz].setY(P4Lab.get(2)*GeV);
      gPfKZ[ikz].setZ(P4Lab.get(3)*GeV);
      ikz++;
    }
    if( type == ks ){
      gPfKZS[ikzs].setX(P4Lab.get(1)*GeV);
      gPfKZS[ikzs].setY(P4Lab.get(2)*GeV);
      gPfKZS[ikzs].setZ(P4Lab.get(3)*GeV);
      ikzs++;
    }
    if( type == kl ){
      gPfKZL[ikzl].setX(P4Lab.get(1)*GeV);
      gPfKZL[ikzl].setY(P4Lab.get(2)*GeV);
      gPfKZL[ikzl].setZ(P4Lab.get(3)*GeV);
      ikzl++;
    }
    if( type == gm ){
      gPfG[igamma].setX(P4Lab.get(1)*GeV);
      gPfG[igamma].setY(P4Lab.get(2)*GeV);
      gPfG[igamma].setZ(P4Lab.get(3)*GeV);
      igamma++;
    }
    if( type == ep ){
      gPfEP[iep].setX(P4Lab.get(1)*GeV);
      gPfEP[iep].setY(P4Lab.get(2)*GeV);
      gPfEP[iep].setZ(P4Lab.get(3)*GeV);
      iep++;
    }
    if( type == em ){
      gPfEM[iem].setX(P4Lab.get(1)*GeV);
      gPfEM[iem].setY(P4Lab.get(2)*GeV);
      gPfEM[iem].setZ(P4Lab.get(3)*GeV);
      iem++;
    }
    if( type == mup ){
      gPfMuP[imup].setX(P4Lab.get(1)*GeV);
      gPfMuP[imup].setY(P4Lab.get(2)*GeV);
      gPfMuP[imup].setZ(P4Lab.get(3)*GeV);
      imup++;
    }
    if( type == mum ){
      gPfMuM[imum].setX(P4Lab.get(1)*GeV);
      gPfMuM[imum].setY(P4Lab.get(2)*GeV);
      gPfMuM[imum].setZ(P4Lab.get(3)*GeV);
      imum++;
    }

    q2=q2->nextIter();
  }while(q2!=0);

  // G4cout<< "**********************" << G4endl;
  // G4cout<< ikp << " " << ikm << G4endl;
  // for( G4int i=0; i<ikp; i++) G4cout<<  gPfKP[i] << G4endl;
  // for( G4int i=0; i<ikm; i++) G4cout<<  gPfKM[i] << G4endl;
  // G4cout<< ipip << " " << ipim << G4endl;
  // for( G4int i=0; i<ipip; i++) G4cout<<  gPfPiP[i] << G4endl;
  // for( G4int i=0; i<ipim; i++) G4cout<<  gPfPiM[i] << G4endl;
  // G4cout<< ip << " " << ipb <<G4endl;
  // for( G4int i=0; i<ip; i++)  G4cout<<  gPfP[i] << G4endl;
  // for( G4int i=0; i<ipb; i++) G4cout<<  gPfPb[i] << G4endl;

  //SetGuns
  for( G4int i=0; i<ikp; i++){
    pGen->gunKP1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKP1_->SetParticleMomentum( gPfKP[i] );
    pGen->gunKP1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikm; i++){
    pGen->gunKM1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKM1_->SetParticleMomentum( gPfKM[i] );
    pGen->gunKM1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipip; i++){
    pGen->gunPiP1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiP1_->SetParticleMomentum( gPfPiP[i] );
    pGen->gunPiP1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipim; i++){
    pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM1_->SetParticleMomentum( gPfPiM[i] );
    pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ip; i++){
    pGen->gunProton_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunProton_->SetParticleMomentum( gPfP[i] );
    pGen->gunProton_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipb; i++){
    pGen->gunProtonb_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunProtonb_->SetParticleMomentum( gPfPb[i] );
    pGen->gunProtonb_->GeneratePrimaryVertex( anEvent );
  }
  // for( G4int i=0; i<ipiz; i++){
  //   pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunPi0_->SetParticleMomentum( gPfPiZ[i] );
  //   pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<in; i++){
  //   pGen->gunNeutron_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunNeutron_->SetParticleMomentum( gPfN[i] );
  //   pGen->gunNeutron_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<inb; i++){
  //   pGen->gunNeutronb_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunNeutronb_->SetParticleMomentum( gPfNb[i] );
  //   pGen->gunNeutronb_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikz; i++){
  //   pGen->gunK0_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0_->SetParticleMomentum( gPfKZ[i] );
  //   pGen->gunK0_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikzs; i++){
  //   pGen->gunK0S_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0S_->SetParticleMomentum( gPfKZS[i] );
  //   pGen->gunK0S_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<ikzl; i++){
  //   pGen->gunK0L_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunK0L_->SetParticleMomentum( gPfKZL[i] );
  //   pGen->gunK0L_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<igamma; i++){
  //   pGen->gunGamma_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunGamma_->SetParticleMomentum( gPfG[i] );
  //   pGen->gunGamma_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<iep; i++){
  //   pGen->gunEP_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunEP_->SetParticleMomentum( gPfEP[i] );
  //   pGen->gunEP_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<iem; i++){
  //   pGen->gunEM_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunEM_->SetParticleMomentum( gPfEM[i] );
  //   pGen->gunEM_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<imup; i++){
  //   pGen->gunMuP_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunMuP_->SetParticleMomentum( gPfMuP[i] );
  //   pGen->gunMuP_->GeneratePrimaryVertex( anEvent );
  // }
  // for( G4int i=0; i<imum; i++){
  //   pGen->gunMuM_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunMuM_->SetParticleMomentum( gPfMuM[i] );
  //   pGen->gunMuM_->GeneratePrimaryVertex( anEvent );
  // }

  if( pGen->anaMan_ ){
    pGen->SetPrimaryInformation( LPos, MF1, MF2, 
                                 LPf1, LPf2, LPf1CM, LPf2CM, pb*LBeamDir, 0. );
  } 
  
  parent1->deleteTree();
  parent2->deleteTree();

  //XiKStarMesonFromKMinusProton3
  //K- p -> Xi*0 K*0
}

