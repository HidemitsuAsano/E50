/*
  PGReaction2.cc 

  2016/1  K.Shirotori
*/

////Charmed Baryon production
////Yc* conditions changed

#include "PGReaction2.hh"
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

const int MaxTry=100;

////////////////////////////////////////
////////////////////
//30: pi- p -> Lambdac+ D*-: FS All 
void PGReaction2::YcDStarMesonFromPiMinusProtonWithDecay1( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  //

  ConfMan *confMan = ConfMan::GetConfManager();

  BeamParam *BP = pGen->BP_;

  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  ////D*- meson
  G4double Mf1=2.010270*GeV;//D*- mass
  static EvtId DS = EvtPDL::getId(std::string("D*-"));

  ////Yc*  
  // G4double Mf2=2.286460*GeV;//Lambda_c(g.s.)+ mass
  // static EvtId Yc = EvtPDL::getId(std::string("Lambda_c+"));

  // G4double Mf2=2.59540*GeV;//Lambda_c(2593)+ mass
  // static EvtId Yc = EvtPDL::getId(std::string("Lambda_c(2593)+"));

  //G4double Mf2=2.6260*GeV;//Lambda_c(2625)+ mass
  //static EvtId Yc = EvtPDL::getId(std::string("Lambda_c(2625)+"));

  // G4double Mf2=2.88150*GeV;//Lambda_c(2880)+ mass
  // static EvtId Yc = EvtPDL::getId(std::string("Lambda_c(2880)+"));

  G4double Mf2=2.9300*GeV;//Lambda_c(2940)+ mass
  static EvtId Yc = EvtPDL::getId(std::string("Lambda_c(2940)+"));

  // G4double Mf2=2.45290*GeV;//Sigma_c(2455)+ mass
  // static EvtId Yc = EvtPDL::getId(std::string("Sigma_c+"));

  // G4double Mf2=2.51750*GeV;//Sigma_c(2520)+ mass
  // static EvtId Yc = EvtPDL::getId(std::string("Sigma_c*+"));

  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;

  AGUniform  genU( BP->cost1, BP->cost2 ); 
  AGForward1 genF( BP->cost1, BP->cost2 );

  G4bool status1 = false;
  if( confMan->GetCMAngle()==1 ){
    status1 = Scattering2Body( Mi1, Mi2, Mf1, Mf2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, genU );
  }
  else if( confMan->GetCMAngle()==2 ){
    status1 = Scattering2Body( Mi1, Mi2, Mf1, Mf2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, genF );
  }

  G4int n1=0;
  while(!status1){
    if(++n1>MaxTry){
      G4Exception("Reaction under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    if( confMan->GetCMAngle()==1 ){
      status1 = Scattering2Body( Mi1, Mi2, Mf1, Mf2, pb*LBeamDir,
				 LPf1, LPf2, LPf1CM, LPf2CM, genU );
    }
    else if( confMan->GetCMAngle()==2 ){
      status1 = Scattering2Body( Mi1, Mi2, Mf1, Mf2, pb*LBeamDir,
				 LPf1, LPf2, LPf1CM, LPf2CM, genF );
    }
  }

  G4ThreeVector gPos(LPos);

  //Beam
  pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM_->GeneratePrimaryVertex( anEvent );

  ///////D*- Yc*+ + decay by EvtGen
  EvtParticle* parent1(0);
  EvtParticle* parent2(0);

  G4LorentzVector LvDS;
  LvDS.setVect(LPf1/GeV);
  LvDS.setE(sqrt(Mf1*Mf1+LPf1.mag2())/GeV);
  
  G4LorentzVector LvYc;
  LvYc.setVect(LPf2/GeV);
  LvYc.setE(sqrt(Mf2*Mf2+LPf2.mag2())/GeV);

  EvtVector4R pInit1( LvDS.e(),
		      LvDS.vect().x(), LvDS.vect().y(), LvDS.vect().z() );
  EvtVector4R pInit2( LvYc.e(),
		      LvYc.vect().x(), LvYc.vect().y(), LvYc.vect().z() );

  parent1 = EvtParticleFactory::particleFactory(DS, pInit1); 
  evtGenerator->generateDecay(parent1);
  parent2 = EvtParticleFactory::particleFactory(Yc, pInit2); 
  evtGenerator->generateDecay(parent2);

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
    pGen->gunKP_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKP_->SetParticleMomentum( gPfKP[i] );
    pGen->gunKP_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikm; i++){
    pGen->gunKM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKM_->SetParticleMomentum( gPfKM[i] );
    pGen->gunKM_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipip; i++){
    pGen->gunPiP_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiP_->SetParticleMomentum( gPfPiP[i] );
    pGen->gunPiP_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipim; i++){
    pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM_->SetParticleMomentum( gPfPiM[i] );
    pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
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
  for( G4int i=0; i<ipiz; i++){
    pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPi0_->SetParticleMomentum( gPfPiZ[i] );
    pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<in; i++){
    pGen->gunNeutron_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunNeutron_->SetParticleMomentum( gPfN[i] );
    pGen->gunNeutron_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<inb; i++){
    pGen->gunNeutronb_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunNeutronb_->SetParticleMomentum( gPfNb[i] );
    pGen->gunNeutronb_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikz; i++){
    pGen->gunK0_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunK0_->SetParticleMomentum( gPfKZ[i] );
    pGen->gunK0_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikzs; i++){
    pGen->gunK0S_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunK0S_->SetParticleMomentum( gPfKZS[i] );
    pGen->gunK0S_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikzl; i++){
    pGen->gunK0L_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunK0L_->SetParticleMomentum( gPfKZL[i] );
    pGen->gunK0L_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<igamma; i++){
    pGen->gunGamma_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunGamma_->SetParticleMomentum( gPfG[i] );
    pGen->gunGamma_->GeneratePrimaryVertex( anEvent );
  }
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
    pGen->SetPrimaryInformation( LPos, Mf1, Mf2, 
                                 LPf1, LPf2, LPf1CM, LPf2CM, pb*LBeamDir, 0. );
  } 

  parent1->deleteTree();
  parent2->deleteTree();

  //YcDStarMesonFromPiMinusProtonWithDecay1
  //30: pi- p -> Lambdac+ D*-: FS All 
}


////////////////////
//31: pi- p -> Yc*+ D*-: FS All 
void PGReaction2::YcDStarMesonFromPiMinusProtonWithDecay2( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  //
 
  ConfMan *confMan = ConfMan::GetConfManager();

  BeamParam *BP = pGen->BP_;

  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  ////D*- meson
  G4double Mf1=2.010270*GeV;//D*- mass
  static EvtId DS = EvtPDL::getId(std::string("D*-"));

  G4int nev=anEvent->GetEventID();
  G4double Mf2; 
  static EvtId Yc;
  switch( nev%7 ){
  case 0:
    Mf2=2.286460*GeV;//Lambda_c(g.s.)+ mass
    Yc = EvtPDL::getId(std::string("Lambda_c+"));
    break;
  case 1:
    Mf2=2.59540*GeV;//Lambda_c(2593)+ mass
    Yc = EvtPDL::getId(std::string("Lambda_c(2593)+"));
    break;
  case 2:
    Mf2=2.6260*GeV;//Lambda_c(2625)+ mass
    Yc = EvtPDL::getId(std::string("Lambda_c(2625)+"));
    break;
  case 3:
    Mf2=2.88150*GeV;//Lambda_c(2880)+ mass
    Yc = EvtPDL::getId(std::string("Lambda_c(2880)+"));
    break;
  case 4:
    Mf2=2.9300*GeV;//Lambda_c(2940)+ mass
    Yc = EvtPDL::getId(std::string("Lambda_c(2940)+"));
    break;
  case 5:
    Mf2=2.45290*GeV;//Sigma_c(2455)+ mass
    Yc = EvtPDL::getId(std::string("Sigma_c+"));
    break;
  case 6:
    Mf2=2.51750*GeV;//Sigma_c(2520)+ mass
    Yc = EvtPDL::getId(std::string("Sigma_c*+"));
    break;
  default:
    Mf2=2.286460*GeV;//Lambda_c(g.s.)+ mass
    Yc = EvtPDL::getId(std::string("Lambda_c+"));
    break;
  }

  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;

  AGUniform  genU( BP->cost1, BP->cost2 ); 
  AGForward1 genF( BP->cost1, BP->cost2 );

  G4bool status1 = false;
  if( confMan->GetCMAngle()==1 ){
    status1 = Scattering2Body( Mi1, Mi2, Mf1, Mf2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, genU );
  }
  else if( confMan->GetCMAngle()==2 ){
    status1 = Scattering2Body( Mi1, Mi2, Mf1, Mf2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, genF );
  }

  G4int n1=0;
  while(!status1){
    if(++n1>MaxTry){
      G4Exception("Reaction under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    if( confMan->GetCMAngle()==1 ){
      status1 = Scattering2Body( Mi1, Mi2, Mf1, Mf2, pb*LBeamDir,
				 LPf1, LPf2, LPf1CM, LPf2CM, genU );
    }
    else if( confMan->GetCMAngle()==2 ){
      status1 = Scattering2Body( Mi1, Mi2, Mf1, Mf2, pb*LBeamDir,
				 LPf1, LPf2, LPf1CM, LPf2CM, genF );
    }
  }

  G4ThreeVector gPos(LPos);

  //Beam
  pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM_->GeneratePrimaryVertex( anEvent );

  ///////D*- Yc*+ + decay by EvtGen
  EvtParticle* parent1(0);
  EvtParticle* parent2(0);

  G4LorentzVector LvDS;
  LvDS.setVect(LPf1/GeV);
  LvDS.setE(sqrt(Mf1*Mf1+LPf1.mag2())/GeV);
  
  G4LorentzVector LvYc;
  LvYc.setVect(LPf2/GeV);
  LvYc.setE(sqrt(Mf2*Mf2+LPf2.mag2())/GeV);

  EvtVector4R pInit1( LvDS.e(),
		      LvDS.vect().x(), LvDS.vect().y(), LvDS.vect().z() );
  EvtVector4R pInit2( LvYc.e(),
		      LvYc.vect().x(), LvYc.vect().y(), LvYc.vect().z() );

  parent1 = EvtParticleFactory::particleFactory(DS, pInit1); 
  evtGenerator->generateDecay(parent1);
  parent2 = EvtParticleFactory::particleFactory(Yc, pInit2); 
  evtGenerator->generateDecay(parent2);

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
    pGen->gunKP_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKP_->SetParticleMomentum( gPfKP[i] );
    pGen->gunKP_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikm; i++){
    pGen->gunKM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKM_->SetParticleMomentum( gPfKM[i] );
    pGen->gunKM_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipip; i++){
    pGen->gunPiP_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiP_->SetParticleMomentum( gPfPiP[i] );
    pGen->gunPiP_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipim; i++){
    pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM_->SetParticleMomentum( gPfPiM[i] );
    pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
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
  for( G4int i=0; i<ipiz; i++){
    pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPi0_->SetParticleMomentum( gPfPiZ[i] );
    pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<in; i++){
    pGen->gunNeutron_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunNeutron_->SetParticleMomentum( gPfN[i] );
    pGen->gunNeutron_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<inb; i++){
    pGen->gunNeutronb_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunNeutronb_->SetParticleMomentum( gPfNb[i] );
    pGen->gunNeutronb_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikz; i++){
    pGen->gunK0_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunK0_->SetParticleMomentum( gPfKZ[i] );
    pGen->gunK0_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikzs; i++){
    pGen->gunK0S_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunK0S_->SetParticleMomentum( gPfKZS[i] );
    pGen->gunK0S_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikzl; i++){
    pGen->gunK0L_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunK0L_->SetParticleMomentum( gPfKZL[i] );
    pGen->gunK0L_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<igamma; i++){
    pGen->gunGamma_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunGamma_->SetParticleMomentum( gPfG[i] );
    pGen->gunGamma_->GeneratePrimaryVertex( anEvent );
  }
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
    pGen->SetPrimaryInformation( LPos, Mf1, Mf2, 
                                 LPf1, LPf2, LPf1CM, LPf2CM, pb*LBeamDir, 0. );
  } 

  parent1->deleteTree();
  parent2->deleteTree();

  //YcDStarMesonFromPiMinusProtonWithDecay2
  //31: pi- p -> Yc*+ D*-: FS All 
}


////////////////////
//32: pi- p -> Lambdac+ D**: FS All 
void PGReaction2::YcDStarMesonFromPiMinusProtonWithDecay3( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  //

  ConfMan *confMan = ConfMan::GetConfManager();

  BeamParam *BP = pGen->BP_;

  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  ////D** mesons
  G4int nev=anEvent->GetEventID();
  G4double Mf1, Mf2, Wid; 
  static EvtId DS, Yc;

  //Lambdac+ & Sigmac0
  switch( nev%6 ){
  case 0:
    Mf1 = 2.42300*GeV;//D_1- mass
    Mf2 = 2.286460*GeV;//Lambda_c(g.s.)+ mass
    Wid = 0.0200*GeV;
    DS = EvtPDL::getId(std::string("D_1-"));  
    Yc = EvtPDL::getId(std::string("Lambda_c+"));
    break;
  case 1:
    Mf1 = 2.4450*GeV;//D'_1- mass
    Mf2 = 2.286460*GeV;//Lambda_c(g.s.)+ mass
    Wid = 0.250*GeV;
    DS = EvtPDL::getId(std::string("D'_1-"));
    Yc = EvtPDL::getId(std::string("Lambda_c+"));
    break;
  case 2:
    Mf1 = 2.46010*GeV;//D_2*- mass
    Mf2 = 2.286460*GeV;//Lambda_c(g.s.)+ mass
    Wid = 0.037*GeV;
    DS = EvtPDL::getId(std::string("D_2*-"));
    Yc = EvtPDL::getId(std::string("Lambda_c+"));
    break;
  case 3:
    Mf1 = 2.42230*GeV;//D_10 mass
    Mf2 = 2.45376*GeV;//Sigma_c(2455)0 mass
    Wid = 0.0204*GeV;
    DS = EvtPDL::getId(std::string("anti-D_10"));
    Yc = EvtPDL::getId(std::string("Sigma_c0"));
    break;
  case 4:
    Mf1 = 2.4450*GeV;//D'_10 mass
    Mf2 = 2.45376*GeV;//Sigma_c(2455)0 mass
    Wid = 0.250*GeV;
    DS = EvtPDL::getId(std::string("anti-D'_10"));
    Yc = EvtPDL::getId(std::string("Sigma_c0"));
    break;
  case 5:
    Mf1 = 2.46110*GeV;//D_2*0 mass
    Mf2 = 2.45376*GeV;//Sigma_c(2455)0 mass
    Wid = 0.043*GeV;
    DS = EvtPDL::getId(std::string("anti-D_2*0"));
    Yc = EvtPDL::getId(std::string("Sigma_c0"));    
    break;
  default:
    Mf1 =2.01027*GeV; 
    Mf2 =2.286460*GeV;//Lambda_c(g.s.)+ mass
    Wid = 0.0*GeV;
    DS = EvtPDL::getId(std::string("D*-"));
    Yc = EvtPDL::getId(std::string("Lambda_c+"));
    break;
  }

  //Lambdac+
  // switch( nev%3 ){
  // case 0:
  //   Mf1 = 2.42300*GeV;//D_1- mass
  //   Mf2 = 2.286460*GeV;//Lambda_c(g.s.)+ mass
  //   Wid = 0.0200*GeV;
  //   DS = EvtPDL::getId(std::string("D_1-"));  
  //   Yc = EvtPDL::getId(std::string("Lambda_c+"));
  //   break;
  // case 1:
  //   Mf1 = 2.4450*GeV;//D'_1- mass
  //   Mf2 = 2.286460*GeV;//Lambda_c(g.s.)+ mass
  //   Wid = 0.250*GeV;
  //   DS = EvtPDL::getId(std::string("D'_1-"));
  //   Yc = EvtPDL::getId(std::string("Lambda_c+"));
  //   break;
  // case 2:
  //   Mf1 = 2.46010*GeV;//D_2*- mass
  //   Mf2 = 2.286460*GeV;//Lambda_c(g.s.)+ mass
  //   Wid = 0.037*GeV;
  //   DS = EvtPDL::getId(std::string("D_2*-"));
  //   Yc = EvtPDL::getId(std::string("Lambda_c+"));
  //   break;
  // // case 3:
  // //   Mf1 = 2.580*GeV;//D(2S)- mass
  // //   Mf2 = 2.286460*GeV;//Lambda_c(g.s.)+ mass
  // //   Wid = 0.0*GeV;
  // //   DS = EvtPDL::getId(std::string("D(2S)-"));
  // //   Yc = EvtPDL::getId(std::string("Lambda_c+"));
  // //   break;
  // // case 4:
  // //   Mf1 = 2.640*GeV;//D*(2S)- mass
  // //   Mf2 = 2.286460*GeV;//Lambda_c(g.s.)+ mass
  // //   Wid = 0.0*GeV;
  // //   DS = EvtPDL::getId(std::string("D*(2S)-"));
  // //   Yc = EvtPDL::getId(std::string("Lambda_c+"));
  // //   break;
  // default:
  //   Mf1 =2.01027*GeV; 
  //   Mf2 =2.286460*GeV;//Lambda_c(g.s.)+ mass
  //   Wid = 0.0*GeV;
  //   DS = EvtPDL::getId(std::string("D*-"));
  //   Yc = EvtPDL::getId(std::string("Lambda_c+"));
  //   break;
  // }

  //Sigmac0
  // switch( nev%3 ){
  // case 0:
  //   Mf1 = 2.42230*GeV;//D_10 mass
  //   Mf2 = 2.45376*GeV;//Sigma_c(2455)0 mass
  //   Wid = 0.0204*GeV;
  //   DS = EvtPDL::getId(std::string("anti-D_10"));
  //   Yc = EvtPDL::getId(std::string("Sigma_c0"));
  //   break;
  // case 1:
  //   Mf1 = 2.4450*GeV;//D'_10 mass
  //   Mf2 = 2.45376*GeV;//Sigma_c(2455)0 mass
  //   Wid = 0.250*GeV;
  //   DS = EvtPDL::getId(std::string("anti-D'_10"));
  //   Yc = EvtPDL::getId(std::string("Sigma_c0"));
  //   break;
  // case 2:
  //   Mf1 = 2.46110*GeV;//D_2*0 mass
  //   Mf2 = 2.45376*GeV;//Sigma_c(2455)0 mass
  //   Wid = 0.043*GeV;
  //   DS = EvtPDL::getId(std::string("anti-D_2*0"));
  //   Yc = EvtPDL::getId(std::string("Sigma_c0"));    
  //   break;
  // // case 3:
  // //   Mf1 = 2.580*GeV;//D(2S)0 mass
  // //   Mf2 = 2.45376*GeV;//Sigma_c(2455)0 mass
  // //   Wid = 0.0*GeV;
  // //   DS = EvtPDL::getId(std::string("anti-D(2S)0"));
  // //   Yc = EvtPDL::getId(std::string("Sigma_c0"));
  // //   break;
  // // case 4:
  // //   Mf1 = 2.640*GeV;//D*(2S)0 mass
  // //   Mf2 = 2.45376*GeV;//Sigma_c(2455)0 mass
  // //   Wid = 0.0*GeV;
  // //   DS = EvtPDL::getId(std::string("anti-D*(2S)0"));
  // //   Yc = EvtPDL::getId(std::string("Sigma_c0"));
  // //   break;
  // default:
  //   Mf1 =2.01027*GeV; 
  //   Mf2 =2.286460*GeV;//Lambda_c(g.s.)+ mass
  //   Wid = 0.0*GeV;
  //   DS = EvtPDL::getId(std::string("D*-"));
  //   Yc = EvtPDL::getId(std::string("Lambda_c+"));
  //   break;
  // }

  // ////Lambdac+
  // G4double Mf2=2.286460*GeV;//Lambda_c(g.s.)+ mass
  // static EvtId Yc = EvtPDL::getId(std::string("Lambda_c+"));

  G4double MF=BreitWigner( Mf1, Wid );

  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4double Ei=sqrt( Mi1*Mi1+Mi2*Mi2+2.*Mi2*sqrt(Mi1*Mi1+pb*pb) );
  G4double Mmin=(2.0103+0.140)*GeV;

  while(1){
    MF=BreitWigner( Mf1, Wid );
    if(MF>Mmin) break;
  }

  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;

  AGUniform  genU( BP->cost1, BP->cost2 ); 
  AGForward1 genF( BP->cost1, BP->cost2 );

  G4bool status1 = false;
  if( confMan->GetCMAngle()==1 ){
    status1 = Scattering2Body( Mi1, Mi2, MF, Mf2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, genU );
  }
  else if( confMan->GetCMAngle()==2 ){
    status1 = Scattering2Body( Mi1, Mi2, MF, Mf2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, genF );
  }

  G4int n1=0;
  while(!status1){
    if(++n1>MaxTry){
      G4Exception("Reaction under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;
    
    while( !(Ei>MF+Mf2) ){
      while(1){
	MF=BreitWigner( Mf1, Wid );
	if(MF>Mmin) break;
      }
    }
    
    if( confMan->GetCMAngle()==1 ){
      status1 = Scattering2Body( Mi1, Mi2, MF, Mf2, pb*LBeamDir,
				 LPf1, LPf2, LPf1CM, LPf2CM, genU );
    }
    else if( confMan->GetCMAngle()==2 ){
      status1 = Scattering2Body( Mi1, Mi2, MF, Mf2, pb*LBeamDir,
				 LPf1, LPf2, LPf1CM, LPf2CM, genF );
    }
  }

  G4ThreeVector gPos(LPos);

  //Beam
  pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM1_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );

  ///////D*- Lambda_c(g.s.)+ decay by EvtGen
  EvtParticle* parent1(0);
  EvtParticle* parent2(0);

  G4LorentzVector LvDS;
  LvDS.setVect(LPf1/GeV);
  LvDS.setE(sqrt(MF*MF+LPf1.mag2())/GeV);
  
  G4LorentzVector LvYc;
  LvYc.setVect(LPf2/GeV);
  LvYc.setE(sqrt(Mf2*Mf2+LPf2.mag2())/GeV);

  EvtVector4R pInit1( LvDS.e(),
		      LvDS.vect().x(), LvDS.vect().y(), LvDS.vect().z() );
  EvtVector4R pInit2( LvYc.e(),
		      LvYc.vect().x(), LvYc.vect().y(), LvYc.vect().z() );

  parent1 = EvtParticleFactory::particleFactory(DS, pInit1); 
  evtGenerator->generateDecay(parent1);
  parent2 = EvtParticleFactory::particleFactory(Yc, pInit2); 
  evtGenerator->generateDecay(parent2);

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
    pGen->gunKP_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKP_->SetParticleMomentum( gPfKP[i] );
    pGen->gunKP_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikm; i++){
    pGen->gunKM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKM_->SetParticleMomentum( gPfKM[i] );
    pGen->gunKM_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipip; i++){
    pGen->gunPiP_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiP_->SetParticleMomentum( gPfPiP[i] );
    pGen->gunPiP_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipim; i++){
    pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM_->SetParticleMomentum( gPfPiM[i] );
    pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
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
  for( G4int i=0; i<ipiz; i++){
    pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPi0_->SetParticleMomentum( gPfPiZ[i] );
    pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<in; i++){
    pGen->gunNeutron_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunNeutron_->SetParticleMomentum( gPfN[i] );
    pGen->gunNeutron_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<inb; i++){
    pGen->gunNeutronb_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunNeutronb_->SetParticleMomentum( gPfNb[i] );
    pGen->gunNeutronb_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikz; i++){
    pGen->gunK0_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunK0_->SetParticleMomentum( gPfKZ[i] );
    pGen->gunK0_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikzs; i++){
    pGen->gunK0S_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunK0S_->SetParticleMomentum( gPfKZS[i] );
    pGen->gunK0S_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikzl; i++){
    pGen->gunK0L_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunK0L_->SetParticleMomentum( gPfKZL[i] );
    pGen->gunK0L_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<igamma; i++){
    pGen->gunGamma_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunGamma_->SetParticleMomentum( gPfG[i] );
    pGen->gunGamma_->GeneratePrimaryVertex( anEvent );
  }
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
    pGen->SetPrimaryInformation( LPos, Mf1, Mf2, 
                                 LPf1, LPf2, LPf1CM, LPf2CM, pb*LBeamDir, 0. );
  } 

  parent1->deleteTree();
  parent2->deleteTree();

  //YcDStarMesonFromPiMinusProtonWithDecay3
  //32: pi- p -> Lambdac+ D**: FS All 
}


////////////////////
//33: pi- p -> Lambdac+ D*- pi0 : FS All 
void PGReaction2::YcDStarMesonFromPiMinusProtonWithDecay4( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  //
 
  BeamParam *BP = pGen->BP_;

  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  ////D*- meson
  G4double Mf1=2.010270*GeV;//D*- mass
  static EvtId DS = EvtPDL::getId(std::string("D*-"));

  // ////Yc*  
  // G4double Mf2=2.286460*GeV;//Lambda_c(g.s.)+ mass
  // static EvtId Yc = EvtPDL::getId(std::string("Lambda_c+"));

  G4int nev=anEvent->GetEventID();
  G4double Mf2; 
  static EvtId Yc;
  switch( nev%5 ){
  case 0:
    Mf2=2.286460*GeV;//Lambda_c(g.s.)+ mass
    Yc = EvtPDL::getId(std::string("Lambda_c+"));
    break;
  case 1:
    Mf2=2.59540*GeV;//Lambda_c(2593)+ mass
    Yc = EvtPDL::getId(std::string("Lambda_c(2593)+"));
    break;
  case 2:
    Mf2=2.6260*GeV;//Lambda_c(2625)+ mass
    Yc = EvtPDL::getId(std::string("Lambda_c(2625)+"));
    break;
  // case 3:
  //   Mf2=2.88150*GeV;//Lambda_c(2880)+ mass
  //   Yc = EvtPDL::getId(std::string("Lambda_c(2880)+"));
  //   break;
  // case 4:
  //   Mf2=2.9300*GeV;//Lambda_c(2940)+ mass
  //   Yc = EvtPDL::getId(std::string("Lambda_c(2940)+"));
  //   break;
  case 3:
    Mf2=2.45290*GeV;//Sigma_c(2455)+ mass
    Yc = EvtPDL::getId(std::string("Sigma_c+"));
    break;
  case 4:
    Mf2=2.51750*GeV;//Sigma_c(2520)+ mass
    Yc = EvtPDL::getId(std::string("Sigma_c*+"));
    break;
  default:
    Mf2=2.286460*GeV;//Lambda_c(g.s.)+ mass
    Yc = EvtPDL::getId(std::string("Lambda_c+"));
    break;
  }

  G4double Mf3=G4PionZero::Definition()->GetPDGMass();//pi0 mass

  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4double Ei=sqrt( Mi1*Mi1+Mi2*Mi2+2.*Mi2*sqrt(Mi1*Mi1+pb*pb) );

  G4ThreeVector LPf1, LPf2, LPf3;
  G4ThreeVector LPf1CM, LPf2CM;

  G4bool status=
    Decay3BodyPhaseSpace( Ei, Mf1, Mf2, Mf3, pb*LBeamDir, LPf1, LPf2, LPf3 );
 
  G4int n=0;
  while(!status){
    if(++n>MaxTry){
      G4Exception("Reaction under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;
    status = 
      Decay3BodyPhaseSpace( Ei, Mf1, Mf2, Mf3, pb*LBeamDir, LPf1, LPf2, LPf3 );
  }

  G4ThreeVector gPos(LPos), gPf3(LPf3);

  //Beam
  pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM_->GeneratePrimaryVertex( anEvent );

  //pi0
  pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPi0_->SetParticleMomentum( gPf3 );
  pGen->gunPi0_->GeneratePrimaryVertex( anEvent );

  ///////D*- Yc*+ + decay by EvtGen
  EvtParticle* parent1(0);
  EvtParticle* parent2(0);

  G4LorentzVector LvDS;
  LvDS.setVect(LPf1/GeV);
  LvDS.setE(sqrt(Mf1*Mf1+LPf1.mag2())/GeV);
  
  G4LorentzVector LvYc;
  LvYc.setVect(LPf2/GeV);
  LvYc.setE(sqrt(Mf2*Mf2+LPf2.mag2())/GeV);

  EvtVector4R pInit1( LvDS.e(),
		      LvDS.vect().x(), LvDS.vect().y(), LvDS.vect().z() );
  EvtVector4R pInit2( LvYc.e(),
		      LvYc.vect().x(), LvYc.vect().y(), LvYc.vect().z() );

  parent1 = EvtParticleFactory::particleFactory(DS, pInit1); 
  evtGenerator->generateDecay(parent1);
  parent2 = EvtParticleFactory::particleFactory(Yc, pInit2); 
  evtGenerator->generateDecay(parent2);

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
    pGen->gunKP_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKP_->SetParticleMomentum( gPfKP[i] );
    pGen->gunKP_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikm; i++){
    pGen->gunKM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunKM_->SetParticleMomentum( gPfKM[i] );
    pGen->gunKM_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipip; i++){
    pGen->gunPiP_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiP_->SetParticleMomentum( gPfPiP[i] );
    pGen->gunPiP_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ipim; i++){
    pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM_->SetParticleMomentum( gPfPiM[i] );
    pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
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
  for( G4int i=0; i<ipiz; i++){
    pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPi0_->SetParticleMomentum( gPfPiZ[i] );
    pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<in; i++){
    pGen->gunNeutron_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunNeutron_->SetParticleMomentum( gPfN[i] );
    pGen->gunNeutron_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<inb; i++){
    pGen->gunNeutronb_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunNeutronb_->SetParticleMomentum( gPfNb[i] );
    pGen->gunNeutronb_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikz; i++){
    pGen->gunK0_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunK0_->SetParticleMomentum( gPfKZ[i] );
    pGen->gunK0_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikzs; i++){
    pGen->gunK0S_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunK0S_->SetParticleMomentum( gPfKZS[i] );
    pGen->gunK0S_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<ikzl; i++){
    pGen->gunK0L_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunK0L_->SetParticleMomentum( gPfKZL[i] );
    pGen->gunK0L_->GeneratePrimaryVertex( anEvent );
  }
  for( G4int i=0; i<igamma; i++){
    pGen->gunGamma_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunGamma_->SetParticleMomentum( gPfG[i] );
    pGen->gunGamma_->GeneratePrimaryVertex( anEvent );
  }
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
    pGen->SetPrimaryInformation( LPos, Mf1, Mf2, 
                                 LPf1, LPf2, LPf1CM, LPf2CM, pb*LBeamDir, 0. );
  } 

  parent1->deleteTree();
  parent2->deleteTree();

  //YcDStarMesonFromPiMinusProtonWithDecay4
  //33: pi- p -> Lambdac*+ D*- pi0 : FS All 
}
