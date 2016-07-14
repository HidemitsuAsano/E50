/*
  PGReaction3.cc 

  2016/1  K.Shirotori
*/

////DDbar pari production
////c-cbar Meson production

#include "PGReaction3.hh"
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
//40: pi- p -> D Dbar N : FS All state & All state (D*- included chs)
void PGReaction3::DStarMesonPairFromPiMinusProton1( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // Uniform angular distribution in CM 
  // Not real distribution
  //

  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  G4bool pflag=false;

  G4int nev=anEvent->GetEventID();
  G4double Mf1, Mf2, Mf3; 
  static EvtId DS1;
  static EvtId DS2;
  switch( nev%4 ){
  case 0:
    Mf1=2.010270*GeV;//D*+ mass
    Mf2=2.010270*GeV;//D*- mass
    Mf3=G4Neutron::Definition()->GetPDGMass();//Neutron mass
    DS1 = EvtPDL::getId(std::string("D*+"));  
    DS2 = EvtPDL::getId(std::string("D*-"));
    break;
  case 1:
    Mf1=2.006970*GeV;//D*0 mass
    Mf2=2.010270*GeV;//D*- mass
    Mf3=G4Proton::Definition()->GetPDGMass();//Proton mass
    DS1 = EvtPDL::getId(std::string("D*0"));  
    DS2 = EvtPDL::getId(std::string("D*-"));
    pflag=true;
    break;
  case 3:
    Mf1=1.869620*GeV;//D+ mass
    Mf2=2.010270*GeV;//D*- mass
    Mf3=G4Neutron::Definition()->GetPDGMass();//Neutron mass
    DS1 = EvtPDL::getId(std::string("D+"));  
    DS2 = EvtPDL::getId(std::string("D*-"));
    break;
  case 4:
    Mf1=1.864840*GeV;//D0 mass
    Mf2=2.010270*GeV;//D*- mass
    Mf3=G4Proton::Definition()->GetPDGMass();//Proton mass
    DS1 = EvtPDL::getId(std::string("D0"));  
    DS2 = EvtPDL::getId(std::string("D*-"));
    pflag=true;
    break;
  default:
    Mf1=2.006970*GeV;//D*0 mass
    Mf2=2.006970*GeV;//D*0bar mass
    Mf3=G4Neutron::Definition()->GetPDGMass();//Neutron mass
    DS1 = EvtPDL::getId(std::string("D*0"));  
    DS2 = EvtPDL::getId(std::string("anti-D*0"));
    break;
  }

  BeamParam *BP = pGen->BP_;
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4double Ei=sqrt( Mi1*Mi1+Mi2*Mi2+2.*Mi2*sqrt(Mi1*Mi1+pb*pb) );
  
  //Decay 
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

  //Nucleon
  if( pflag ){
    pGen->gunProton_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunProton_->SetParticleMomentum( gPf3 );
    pGen->gunProton_->GeneratePrimaryVertex( anEvent );
  }
  else{
    pGen->gunNeutron_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunNeutron_->SetParticleMomentum( gPf3 );
    pGen->gunNeutron_->GeneratePrimaryVertex( anEvent );
  }

  ///////D D decay by EvtGen
  EvtParticle* parent1(0);
  EvtParticle* parent2(0);

  G4LorentzVector LvDS1;
  LvDS1.setVect(LPf1/GeV);
  LvDS1.setE(sqrt(Mf1*Mf1+LPf1.mag2())/GeV);

  G4LorentzVector LvDS2;
  LvDS2.setVect(LPf2/GeV);
  LvDS2.setE(sqrt(Mf2*Mf2+LPf2.mag2())/GeV);

  EvtVector4R pInit1( LvDS1.e(),
		      LvDS1.vect().x(), LvDS1.vect().y(), LvDS1.vect().z() );
  EvtVector4R pInit2( LvDS2.e(),
		      LvDS2.vect().x(), LvDS2.vect().y(), LvDS2.vect().z() );

  parent1 = EvtParticleFactory::particleFactory(DS1, pInit1); 
  evtGenerator->generateDecay(parent1);
  parent2 = EvtParticleFactory::particleFactory(DS2, pInit2); 
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

  //DStarMesonPairFromPiMinusProton1
  //40: pi- p -> D Dbar N : FS All state & All state
}


////////////////////
//41: pi- p -> M(ccbar) n : FS All state
void PGReaction3::DStarMesonPairFromPiMinusProton2( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) angular distribution in CM 
  // Not real distribution
  //

  ConfMan *confMan = ConfMan::GetConfManager();

  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  G4int nev=anEvent->GetEventID();
  G4double Mf1,Wid; 
  static EvtId M;
  switch( nev%12 ){
  case 0:
    Mf1 = 2.9803*GeV;//eta_c 
    Wid = 0.0276*GeV;
    M = EvtPDL::getId(std::string("eta_c"));
    break;
  case 1:
    Mf1=3.6370*GeV;//eta_c(2S)
    Wid = 0.0139*GeV;
    M = EvtPDL::getId(std::string("eta_c(2S)"));
    break;
  case 2:
    Mf1=3.096916*GeV;//J/psi
    Wid = 0.0*GeV;
    M = EvtPDL::getId(std::string("J/psi"));
    break;
  case 3:
    Mf1=3.68609*GeV;//psi(2S)
    Wid = 0.000317*GeV;
    M = EvtPDL::getId(std::string("psi(2S)"));
    break;
  case 4:
    Mf1=3.77292*GeV;//psi(3770)
    Wid = 0.0272*GeV;
    M = EvtPDL::getId(std::string("psi(3770)"));
    break;
  case 5:
    Mf1=4.0390*GeV;//psi(4040)
    Wid = 0.080*GeV;
    M = EvtPDL::getId(std::string("psi(4040)"));
    break;
  case 6:
    Mf1=4.1530*GeV;//psi(4160)
    Wid = 0.078*GeV;
    M = EvtPDL::getId(std::string("psi(4160)"));
    break;
  case 7:
    Mf1=4.4210*GeV;//psi(4415)
    Wid = 0.062*GeV;
    M = EvtPDL::getId(std::string("psi(4415)"));
    break;
  case 8:
    Mf1=3.41475*GeV;//chi_c0
    Wid = 0.0102*GeV;
    M = EvtPDL::getId(std::string("chi_c0"));
    break;
  case 9:
    Mf1=3.51066*GeV;//chi_c1
    Wid = 0.00089*GeV;
    M = EvtPDL::getId(std::string("chi_c1"));
    break;
  case 10:
    Mf1=3.55620*GeV;//chi_c2
    Wid = 0.00203*GeV;
    M = EvtPDL::getId(std::string("chi_c2"));
    break;
  case 11:
    Mf1=3.525930*GeV;//h_c
    Wid = 0.0*GeV;
    M = EvtPDL::getId(std::string("h_c"));
    break;
  default:
    Mf1=3.096916*GeV;//J/psi
    Wid = 0.0*GeV;
    M = EvtPDL::getId(std::string("J/psi"));
    break;
  }

  G4double Mf2=G4Neutron::Definition()->GetPDGMass();//Neutron mass

  G4double MF=BreitWigner( Mf1, Wid );
 
  BeamParam *BP = pGen->BP_;
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4double Ei=sqrt( Mi1*Mi1+Mi2*Mi2+2.*Mi2*sqrt(Mi1*Mi1+pb*pb) );
  G4double Mmin=2.5*GeV;

  if(MF<Mmin){
    while(1){
      MF=BreitWigner( Mf1, Wid );
      if(MF>Mmin) break;
    }
  }
  // G4cout<< "****1" << G4endl;
  // G4cout<< MF << G4endl;

  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;

  AGUniform  genU( BP->cost1, BP->cost2 );
  AGForward1 genF( BP->cost1, BP->cost2 );

  G4bool status = false;
  if( confMan->GetCMAngle()==1 ){
    status = Scattering2Body( Mi1, Mi2, MF, Mf2, pb*LBeamDir,
			      LPf1, LPf2, LPf1CM, LPf2CM, genU );
  }
  else if( confMan->GetCMAngle()==2 ){
    status = Scattering2Body( Mi1, Mi2, MF, Mf2, pb*LBeamDir,
			      LPf1, LPf2, LPf1CM, LPf2CM, genF );
  }

  G4int n=0;
  while(!status){
    if(++n>MaxTry){
      G4Exception("Reaction under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    while( !(Ei>MF+Mf2) ){
      MF=BreitWigner( Mf1, Wid );
      if(MF<Mmin){
	while(1){
	  MF=BreitWigner( Mf1, Wid );
	  if(MF>Mmin) break;
	}
      }
    }

    if( confMan->GetCMAngle()==1 ){
      status = Scattering2Body( Mi1, Mi2, MF, Mf2, pb*LBeamDir,
				LPf1, LPf2, LPf1CM, LPf2CM, genU );
    }
    else if( confMan->GetCMAngle()==2 ){
      status = Scattering2Body( Mi1, Mi2, MF, Mf2, pb*LBeamDir,
				LPf1, LPf2, LPf1CM, LPf2CM, genF );
    }
  }

  // G4cout<< "****2" << G4endl;
  // G4cout<< MF << " " << LPf1.mag() << G4endl;

  G4ThreeVector gPos(LPos);

  //Beam
  pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM_->GeneratePrimaryVertex( anEvent );

  ///////M(ccbar) decay by EvtGen
  EvtParticle* parent(0);

  G4LorentzVector LvM;
  LvM.setVect(LPf1/GeV);
  LvM.setE(sqrt(MF*MF+LPf1.mag2())/GeV);

  EvtVector4R pInit( LvM.e(),
		     LvM.vect().x(), LvM.vect().y(), LvM.vect().z() );

  // G4cout<< "****" << G4endl;
  // G4cout<< LvM.e() << " " << LvM.mag() << G4endl;

  parent = EvtParticleFactory::particleFactory(M, pInit); 
  evtGenerator->generateDecay(parent);

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

  EvtParticle *p;
  p=parent->getDaug(0);

  do{
    EvtId type=p->getId();
    
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
    
    p=p->nextIter();
  }while(p!=0);

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

  EvtParticle *q;
  q=parent->getDaug(0);

  do{
    EvtId type=q->getId();
    EvtVector4R P4Lab=q->getP4Lab();

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

    q=q->nextIter();
  }while(q!=0);

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

  parent->deleteTree();

  //DStarMesonPairFromPiMinusProton2
  //41: pi- p -> M(ccbar) n : FS All state
}

