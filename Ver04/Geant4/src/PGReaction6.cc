/*
  PGReaction6.cc 

  2016/1  K.Shirotori
*/

////Hyperon production
////Scattered particle (K*) conditions changed

#include "PGReaction6.hh"
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
//60: pi- p -> Y*0 K*0 for Lambda
void PGReaction6::HyperonRKStarZFromPiMinusProton1( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  // 

  ConfMan *confMan = ConfMan::GetConfManager();

  BeamParam *BP = pGen->BP_;

  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  //K*0
  G4double Mf1=0.8960*GeV;//K*0 mass
  G4double Wid = 0.05030*GeV;
  static EvtId KS1 = EvtPDL::getId(std::string("K*0"));  
  G4double MF1=BreitWigner( Mf1, Wid );

  //Hyperons
  G4double Mass=BP->RecMass;
  G4double Width=BP->RecWidth;

  static EvtId HY1;
  if( 1.110*GeV<Mass && Mass<1.120*GeV ) HY1 = EvtPDL::getId(std::string("Lambda0")); 
  else if( 1.400*GeV<Mass && Mass<1.410*GeV ) HY1 = EvtPDL::getId(std::string("Lambda(1405)0")); 
  else if( 1.515*GeV<Mass && Mass<1.525*GeV ) HY1 = EvtPDL::getId(std::string("Lambda(1520)0")); 
  else if( 1.595*GeV<Mass && Mass<1.605*GeV ) HY1 = EvtPDL::getId(std::string("Lambda(1600)0")); 
  else if( 1.665*GeV<Mass && Mass<1.675*GeV ) HY1 = EvtPDL::getId(std::string("Lambda(1670)0")); 
  else if( 1.685*GeV<Mass && Mass<1.695*GeV ) HY1 = EvtPDL::getId(std::string("Lambda(1690)0")); 
  else if( 1.795*GeV<Mass && Mass<1.805*GeV ) HY1 = EvtPDL::getId(std::string("Lambda(1800)0")); 
  else if( 1.805*GeV<Mass && Mass<1.815*GeV ) HY1 = EvtPDL::getId(std::string("Lambda(1810)0")); 
  else if( 1.815*GeV<Mass && Mass<1.825*GeV ) HY1 = EvtPDL::getId(std::string("Lambda(1820)0")); 
  else if( 1.825*GeV<Mass && Mass<1.835*GeV ) HY1 = EvtPDL::getId(std::string("Lambda(1830)0")); 
  else if( 1.885*GeV<Mass && Mass<1.895*GeV ) HY1 = EvtPDL::getId(std::string("Lambda(1890)0")); 
  else{
    Mass = 1.11568*GeV;
    Width = 0.0;
    HY1 = EvtPDL::getId(std::string("Lambda0")); 
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
  G4double MYmin=(1.189+0.140)*GeV;//Sigma + pi

  if( Mass<1.120*GeV ){
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

  AGUniform  genU( BP->cost1, BP->cost2 ); 
  AGForward2 genF( BP->cost1, BP->cost2 );

  G4bool status1 = false;
  if( confMan->GetCMAngle()==1 ){
    status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, genU );
  }
  else if( confMan->GetCMAngle()==2 ){
    status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, genF );
  }

  G4int n=0;
  while(!status1){
    if(++n>MaxTry){
      G4Exception("HyperonRKStarZFromPiMinusProton1: under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    while( !(Ei>MF1+MF2) ){
      if( Mass<1.120*GeV ){
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

    if( confMan->GetCMAngle()==1 ){
      status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
				 LPf1, LPf2, LPf1CM, LPf2CM, genU );
    }
    else if( confMan->GetCMAngle()==2 ){
      status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
				 LPf1, LPf2, LPf1CM, LPf2CM, genF );
    }
  }

  G4ThreeVector gPos(LPos);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  LBeamDir.rotateX(90.0*degree); LBeamDir.rotateZ( BP->TargetRotAngH );

  //Beam
  pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM_->GeneratePrimaryVertex( anEvent );

  ///////K*0 && L* decay by EvtGen
  EvtParticle* parent1(0);
  EvtParticle* parent2(0);
  
  G4LorentzVector LvKS;
  LvKS.setVect(LPf1/GeV);
  LvKS.setE(sqrt(MF1*MF1+LPf1.mag2())/GeV);
  
  G4LorentzVector LvHY;
  LvHY.setVect(LPf2/GeV);
  LvHY.setE(sqrt(MF2*MF2+LPf2.mag2())/GeV);
  
  EvtVector4R pInit1( LvKS.e(),
		      LvKS.vect().x(), LvKS.vect().y(), LvKS.vect().z() );
  EvtVector4R pInit2( LvHY.e(),
		      LvHY.vect().x(), LvHY.vect().y(), LvHY.vect().z() );

  parent1 = EvtParticleFactory::particleFactory(KS1, pInit1); 
  evtGenerator->generateDecay(parent1);
  parent2 = EvtParticleFactory::particleFactory(HY1, pInit2); 
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

  for( G4int i=0; i<ikp; i++){
    gPfKP[i].rotateX(90.0*degree); gPfKP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ikm; i++){
    gPfKM[i].rotateX(90.0*degree); gPfKM[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ipip; i++){
    gPfPiP[i].rotateX(90.0*degree); gPfPiP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ipim; i++){
    gPfPiM[i].rotateX(90.0*degree); gPfPiM[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ip; i++){
    gPfP[i].rotateX(90.0*degree); gPfP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ipb; i++){
    gPfPb[i].rotateX(90.0*degree); gPfPb[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ipiz; i++){
    gPfPiZ[i].rotateX(90.0*degree); gPfPiZ[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<in; i++){
    gPfN[i].rotateX(90.0*degree); gPfN[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<inb; i++){
    gPfNb[i].rotateX(90.0*degree); gPfNb[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ikz; i++){
    gPfKZ[i].rotateX(90.0*degree); gPfKZ[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ikzs; i++){
    gPfKZS[i].rotateX(90.0*degree); gPfKZS[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ikzl; i++){
    gPfKZL[i].rotateX(90.0*degree); gPfKZL[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<igamma; i++){
    gPfG[i].rotateX(90.0*degree); gPfG[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<iep; i++){
    gPfEP[i].rotateX(90.0*degree); gPfEP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<iem; i++){
    gPfEM[i].rotateX(90.0*degree); gPfEM[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<imup; i++){
    gPfMuP[i].rotateX(90.0*degree); gPfMuP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<imum; i++){
    gPfMuM[i].rotateX(90.0*degree); gPfMuM[i].rotateZ( BP->TargetRotAngH );
  }

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
    LBeamDir.rotateZ( -1.*(BP->TargetRotAngH) ); 
    LBeamDir.rotateX(-90.0*degree);
    pGen->SetPrimaryInformation( LPos, MF1, MF2, 
                                 LPf1, LPf2, LPf1CM, LPf2CM, pb*LBeamDir, 0. );
  } 
  
  parent1->deleteTree();
  parent2->deleteTree();

  //HyperonRKStarZFromPiMinusProton1
  //60: pi- p -> Y*0 K*0 for Lambda
}


////////////////////////////////////////
////////////////////
//61: pi- p -> Y*0 K*0 for Sigma0
void PGReaction6::HyperonRKStarZFromPiMinusProton2( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  // 

  ConfMan *confMan = ConfMan::GetConfManager();

  BeamParam *BP = pGen->BP_;

  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  //K*0
  G4double Mf1=0.8960*GeV;//K*0 mass
  G4double Wid = 0.05030*GeV;
  static EvtId KS1 = EvtPDL::getId(std::string("K*0"));  
  G4double MF1=BreitWigner( Mf1, Wid );

  //Hyperons
  G4double Mass=BP->RecMass;
  G4double Width=BP->RecWidth;

  static EvtId HY1;
  if( 1.190*GeV<Mass && Mass<1.195*GeV ) HY1 = EvtPDL::getId(std::string("Sigma0")); 
  else if( 1.380*GeV<Mass && Mass<1.390*GeV ) HY1 = EvtPDL::getId(std::string("Sigma*0")); 
  else if( 1.475*GeV<Mass && Mass<1.485*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1480)0")); 
  else if( 1.555*GeV<Mass && Mass<1.565*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1560)0")); 
  else if( 1.575*GeV<Mass && Mass<1.585*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1580)0")); 
  else if( 1.615*GeV<Mass && Mass<1.625*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1620)0")); 
  else if( 1.655*GeV<Mass && Mass<1.665*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1660)0")); 
  else if( 1.665*GeV<Mass && Mass<1.673*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1670)0")); 
  else if( 1.674*GeV<Mass && Mass<1.679*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1675)0")); 
  else if( 1.685*GeV<Mass && Mass<1.695*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1690)0")); 
  else if( 1.745*GeV<Mass && Mass<1.755*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1750)0")); 
  else if( 1.765*GeV<Mass && Mass<1.773*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1770)0")); 
  else if( 1.774*GeV<Mass && Mass<1.779*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1775)0")); 
  else if( 1.835*GeV<Mass && Mass<1.845*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1840)0")); 
  else if( 1.875*GeV<Mass && Mass<1.885*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1880)0")); 
  else if( 1.910*GeV<Mass && Mass<1.920*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1915)0")); 
  else if( 1.935*GeV<Mass && Mass<1.945*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1940)0")); 
  else if( 1.995*GeV<Mass && Mass<2.005*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(2000)0")); 
  else if( 2.025*GeV<Mass && Mass<2.035*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(2030)0")); 
  else{
    Mass = 1.19264*GeV;
    Width = 0.0;
    HY1 = EvtPDL::getId(std::string("Sigma0")); 
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
  G4double MYmin=(1.189+0.140)*GeV;//Sigma + pi
  
  if( Mass<1.195*GeV ){
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

  AGUniform  genU( BP->cost1, BP->cost2 ); 
  AGForward2 genF( BP->cost1, BP->cost2 );

  G4bool status1 = false;
  if( confMan->GetCMAngle()==1 ){
    status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, genU );
  }
  else if( confMan->GetCMAngle()==2 ){
    status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, genF );
  }

  G4int n=0;
  while(!status1){
    if(++n>MaxTry){
      G4Exception("HyperonRKStarZFromPiMinusProton1: under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    while( !(Ei>MF1+MF2) ){
      if( Mass<1.195*GeV ){
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

    if( confMan->GetCMAngle()==1 ){
      status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
				 LPf1, LPf2, LPf1CM, LPf2CM, genU );
    }
    else if( confMan->GetCMAngle()==2 ){
      status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
				 LPf1, LPf2, LPf1CM, LPf2CM, genF );
    }
  }

  G4ThreeVector gPos(LPos);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  LBeamDir.rotateX(90.0*degree); LBeamDir.rotateZ( BP->TargetRotAngH );

  //Beam
  pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM_->GeneratePrimaryVertex( anEvent );

  ///////K*0 && L* decay by EvtGen
  EvtParticle* parent1(0);
  EvtParticle* parent2(0);
  
  G4LorentzVector LvKS;
  LvKS.setVect(LPf1/GeV);
  LvKS.setE(sqrt(MF1*MF1+LPf1.mag2())/GeV);
  
  G4LorentzVector LvHY;
  LvHY.setVect(LPf2/GeV);
  LvHY.setE(sqrt(MF2*MF2+LPf2.mag2())/GeV);
  
  EvtVector4R pInit1( LvKS.e(),
		      LvKS.vect().x(), LvKS.vect().y(), LvKS.vect().z() );
  EvtVector4R pInit2( LvHY.e(),
		      LvHY.vect().x(), LvHY.vect().y(), LvHY.vect().z() );

  parent1 = EvtParticleFactory::particleFactory(KS1, pInit1); 
  evtGenerator->generateDecay(parent1);
  parent2 = EvtParticleFactory::particleFactory(HY1, pInit2); 
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

  for( G4int i=0; i<ikp; i++){
    gPfKP[i].rotateX(90.0*degree); gPfKP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ikm; i++){
    gPfKM[i].rotateX(90.0*degree); gPfKM[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ipip; i++){
    gPfPiP[i].rotateX(90.0*degree); gPfPiP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ipim; i++){
    gPfPiM[i].rotateX(90.0*degree); gPfPiM[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ip; i++){
    gPfP[i].rotateX(90.0*degree); gPfP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ipb; i++){
    gPfPb[i].rotateX(90.0*degree); gPfPb[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ipiz; i++){
    gPfPiZ[i].rotateX(90.0*degree); gPfPiZ[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<in; i++){
    gPfN[i].rotateX(90.0*degree); gPfN[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<inb; i++){
    gPfNb[i].rotateX(90.0*degree); gPfNb[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ikz; i++){
    gPfKZ[i].rotateX(90.0*degree); gPfKZ[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ikzs; i++){
    gPfKZS[i].rotateX(90.0*degree); gPfKZS[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ikzl; i++){
    gPfKZL[i].rotateX(90.0*degree); gPfKZL[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<igamma; i++){
    gPfG[i].rotateX(90.0*degree); gPfG[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<iep; i++){
    gPfEP[i].rotateX(90.0*degree); gPfEP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<iem; i++){
    gPfEM[i].rotateX(90.0*degree); gPfEM[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<imup; i++){
    gPfMuP[i].rotateX(90.0*degree); gPfMuP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<imum; i++){
    gPfMuM[i].rotateX(90.0*degree); gPfMuM[i].rotateZ( BP->TargetRotAngH );
  }

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
    LBeamDir.rotateZ( -1.*(BP->TargetRotAngH) ); 
    LBeamDir.rotateX(-90.0*degree);
    pGen->SetPrimaryInformation( LPos, MF1, MF2, 
                                 LPf1, LPf2, LPf1CM, LPf2CM, pb*LBeamDir, 0. );
  } 
  
  parent1->deleteTree();
  parent2->deleteTree();

  //HyperonRKStarZFromPiMinusProton2
  //61: pi- p -> Y*0 K*0 for Sigma0
}


////////////////////////////////////////
////////////////////
//62: pi- p -> Y*- K*+ for Sigma-
void PGReaction6::HyperonRKStarPFromPiMinusProton3( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  // 

  ConfMan *confMan = ConfMan::GetConfManager();

  BeamParam *BP = pGen->BP_;

  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  //K*+
  G4double Mf1=0.89166*GeV;//K*+ mass
  G4double Wid = 0.05080*GeV;
  static EvtId KS1 = EvtPDL::getId(std::string("K*+"));  
  G4double MF1=BreitWigner( Mf1, Wid );

  //Hyperons
  G4double Mass=BP->RecMass;
  G4double Width=BP->RecWidth;

  static EvtId HY1;
  if( 1.195*GeV<Mass && Mass<1.199*GeV ) HY1 = EvtPDL::getId(std::string("Sigma-")); 
  else if( 1.380*GeV<Mass && Mass<1.390*GeV ) HY1 = EvtPDL::getId(std::string("Sigma*-")); 
  else if( 1.475*GeV<Mass && Mass<1.485*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1480)-")); 
  else if( 1.555*GeV<Mass && Mass<1.565*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1560)-")); 
  else if( 1.575*GeV<Mass && Mass<1.585*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1580)-")); 
  else if( 1.615*GeV<Mass && Mass<1.625*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1620)-")); 
  else if( 1.655*GeV<Mass && Mass<1.665*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1660)-")); 
  else if( 1.665*GeV<Mass && Mass<1.673*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1670)-")); 
  else if( 1.674*GeV<Mass && Mass<1.679*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1675)-")); 
  else if( 1.685*GeV<Mass && Mass<1.695*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1690)-")); 
  else if( 1.745*GeV<Mass && Mass<1.755*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1750)-")); 
  else if( 1.765*GeV<Mass && Mass<1.773*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1770)-")); 
  else if( 1.774*GeV<Mass && Mass<1.779*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1775)-")); 
  else if( 1.835*GeV<Mass && Mass<1.845*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1840)-")); 
  else if( 1.875*GeV<Mass && Mass<1.885*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1880)-")); 
  else if( 1.910*GeV<Mass && Mass<1.920*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1915)-")); 
  else if( 1.935*GeV<Mass && Mass<1.945*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(1940)-")); 
  else if( 1.995*GeV<Mass && Mass<2.005*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(2000)-")); 
  else if( 2.025*GeV<Mass && Mass<2.035*GeV ) HY1 = EvtPDL::getId(std::string("Sigma(2030)-")); 
  else{
    Mass = 1.19745*GeV;
    Width = 0.0;
    HY1 = EvtPDL::getId(std::string("Sigma-")); 
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
  G4double MYmin=(1.189+0.140)*GeV;//Sigma + pi

 if( Mass<1.199*GeV ){
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

  AGUniform  genU( BP->cost1, BP->cost2 ); 
  AGForward2 genF( BP->cost1, BP->cost2 );

  G4bool status1=false;
  if( confMan->GetCMAngle()==1 ){
    status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, genU );
  }
  else if( confMan->GetCMAngle()==2 ){
    status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, genF );
  }

  G4int n=0;
  while(!status1){
    if(++n>MaxTry){
      G4Exception("HyperonRKStarZFromPiMinusProton1: under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    while( !(Ei>MF1+MF2) ){
      if( Mass<1.199*GeV ){
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

    if( confMan->GetCMAngle()==1 ){
      status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
				 LPf1, LPf2, LPf1CM, LPf2CM, genU );
    }
    else if( confMan->GetCMAngle()==2 ){
      status1 = Scattering2Body( Mi1, Mi2, MF1, MF2, pb*LBeamDir,
			       LPf1, LPf2, LPf1CM, LPf2CM, genF );
    }
  }

  G4ThreeVector gPos(LPos);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  LBeamDir.rotateX(90.0*degree); LBeamDir.rotateZ( BP->TargetRotAngH );

  //Beam
  pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM_->GeneratePrimaryVertex( anEvent );

  ///////K*+ && L* decay by EvtGen
  EvtParticle* parent1(0);
  EvtParticle* parent2(0);
  
  G4LorentzVector LvKS;
  LvKS.setVect(LPf1/GeV);
  LvKS.setE(sqrt(MF1*MF1+LPf1.mag2())/GeV);
  
  G4LorentzVector LvHY;
  LvHY.setVect(LPf2/GeV);
  LvHY.setE(sqrt(MF2*MF2+LPf2.mag2())/GeV);
  
  EvtVector4R pInit1( LvKS.e(),
		      LvKS.vect().x(), LvKS.vect().y(), LvKS.vect().z() );
  EvtVector4R pInit2( LvHY.e(),
		      LvHY.vect().x(), LvHY.vect().y(), LvHY.vect().z() );

  parent1 = EvtParticleFactory::particleFactory(KS1, pInit1); 
  evtGenerator->generateDecay(parent1);
  parent2 = EvtParticleFactory::particleFactory(HY1, pInit2); 
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

  for( G4int i=0; i<ikp; i++){
    gPfKP[i].rotateX(90.0*degree); gPfKP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ikm; i++){
    gPfKM[i].rotateX(90.0*degree); gPfKM[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ipip; i++){
    gPfPiP[i].rotateX(90.0*degree); gPfPiP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ipim; i++){
    gPfPiM[i].rotateX(90.0*degree); gPfPiM[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ip; i++){
    gPfP[i].rotateX(90.0*degree); gPfP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ipb; i++){
    gPfPb[i].rotateX(90.0*degree); gPfPb[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ipiz; i++){
    gPfPiZ[i].rotateX(90.0*degree); gPfPiZ[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<in; i++){
    gPfN[i].rotateX(90.0*degree); gPfN[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<inb; i++){
    gPfNb[i].rotateX(90.0*degree); gPfNb[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ikz; i++){
    gPfKZ[i].rotateX(90.0*degree); gPfKZ[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ikzs; i++){
    gPfKZS[i].rotateX(90.0*degree); gPfKZS[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<ikzl; i++){
    gPfKZL[i].rotateX(90.0*degree); gPfKZL[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<igamma; i++){
    gPfG[i].rotateX(90.0*degree); gPfG[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<iep; i++){
    gPfEP[i].rotateX(90.0*degree); gPfEP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<iem; i++){
    gPfEM[i].rotateX(90.0*degree); gPfEM[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<imup; i++){
    gPfMuP[i].rotateX(90.0*degree); gPfMuP[i].rotateZ( BP->TargetRotAngH );
  }
  for( G4int i=0; i<imum; i++){
    gPfMuM[i].rotateX(90.0*degree); gPfMuM[i].rotateZ( BP->TargetRotAngH );
  }

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
    LBeamDir.rotateZ( -1.*(BP->TargetRotAngH) ); 
    LBeamDir.rotateX(-90.0*degree);
    pGen->SetPrimaryInformation( LPos, MF1, MF2, 
                                 LPf1, LPf2, LPf1CM, LPf2CM, pb*LBeamDir, 0. );
  } 
  
  parent1->deleteTree();
  parent2->deleteTree();

  //HyperonRKStarPFromPiMinusProton3
  //62: pi- p -> Y*- K*+ for Sigma-
}
