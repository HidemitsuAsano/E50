/*
  PGReaction1.cc 

  2016/1  K.Shirotori
*/

////Charmed Baryon production
////Scattered particle (D*) conditions changed

#include "PGReaction1.hh"
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
//20: pi- p -> Yc+ D*- : FS K+ pi- pi- 
void PGReaction1::YcDStarMesonFromPiMinusProton1( G4Event *anEvent )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  //

  ConfMan *confMan = ConfMan::GetConfManager();

  BeamParam *BP = pGen->BP_;
  G4double Mass=BP->RecMass;
  G4double Width=BP->RecWidth;
 
  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  G4double Mm1=2.01022*GeV;//D*- mass
  G4double Mm2=BreitWigner( Mass, Width );

  G4double Mf1=G4DMesonZero::Definition()->GetPDGMass();
  G4double Mf2=G4PionMinus::Definition()->GetPDGMass();

  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4double Ei=sqrt( Mi1*Mi1+Mi2*Mi2+2.*Mi2*sqrt(Mi1*Mi1+pb*pb) );

  G4ThreeVector LPm1, LPm2;
  G4ThreeVector LPm1CM, LPm2CM;

  AGUniform  genU( BP->cost1, BP->cost2 );
  AGForward1 genF( BP->cost1, BP->cost2 );

  G4bool status1 = false;
  if( confMan->GetCMAngle()==1 ){ 
    status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
			       LPm1, LPm2, LPm1CM, LPm2CM, genU );
  }
  else if( confMan->GetCMAngle()==2 ){ 
    status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
			       LPm1, LPm2, LPm1CM, LPm2CM, genF );
  }

  G4int n1=0;
  while(!status1){
    if(++n1>MaxTry){
      G4Exception("Reaction under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    while( !(Ei>Mm1+Mm2) ) Mm2=BreitWigner( Mass, Width );

    if( confMan->GetCMAngle()==1 ){ 
      status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
				 LPm1, LPm2, LPm1CM, LPm2CM, genU );
    }
    else if( confMan->GetCMAngle()==2 ){ 
      status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
				 LPm1, LPm2, LPm1CM, LPm2CM, genF );
    }
  }
  
  //D*- decay 
  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;
  AGUniform gen2;
  G4bool status2 = Decay2Body( Mm1, Mf1, Mf2, LPm1, LPf1, LPf2, LPf1CM, LPf2CM, gen2 );
  if(!status2){
    G4Exception("Reaction under Threshold");
  }

  G4ThreeVector gPos(LPos), gPf1(LPf1), gPf2(LPf2);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  gPf2.rotateX(90.0*degree); gPf2.rotateZ( BP->TargetRotAngH );
  LBeamDir.rotateX(90.0*degree); LBeamDir.rotateZ( BP->TargetRotAngH );

  pGen->gunDMesonZ1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunDMesonZ1_->SetParticleMomentum( gPf1 );
  pGen->gunDMesonZ1_->GeneratePrimaryVertex( anEvent );
  
  pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM1_->SetParticleMomentum( gPf2 );
  pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );
  // pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  // pGen->gunPiM_->SetParticleMomentum( gPf2 );
  // pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
  
  //Beam
  pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM1_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );
  
  if( pGen->anaMan_ ){
    LBeamDir.rotateZ( -1.*(BP->TargetRotAngH) ); 
    LBeamDir.rotateX(-90.0*degree);
    gPf2.rotateZ( -1.*(BP->TargetRotAngH) ); 
    gPf2.rotateX(-90.0*degree);
    pGen->SetPrimaryInformation( LPos, Mm1, Mm2, 
                                 LPm1, LPm2, LPm1CM, LPm2CM, pb*LBeamDir, gPf2 );
  }
  //YcDStarMesonFromPiMinusProton1
  //20: pi- p -> Yc+ D*- : FS K+ pi- pi- 
}


////////////////////
//21: pi- p -> Yc+ D*- : FS All state 
void PGReaction1::YcDStarMesonFromPiMinusProton2( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  //

  ConfMan *confMan = ConfMan::GetConfManager();

  BeamParam *BP = pGen->BP_;
  G4double Mass=BP->RecMass;
  G4double Width=BP->RecWidth;
 
  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  G4double Mm1=2.0102700*GeV;//D*- mass
  //G4double Mm1=2.006970*GeV;//D*0 mass
  //G4double Mm1=1.8648400*GeV;//D0 mass
  //G4double Mm1=1.8696200*GeV;//D- mass

  G4double Mm2=BreitWigner( Mass, Width );

  G4double Mf1=G4DMesonZero::Definition()->GetPDGMass();
  G4double Mf2=G4PionMinus::Definition()->GetPDGMass();
  // G4double Mf1=G4DMesonMinus::Definition()->GetPDGMass();
  // G4double Mf2=G4PionZero::Definition()->GetPDGMass();

  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4double Ei=sqrt( Mi1*Mi1+Mi2*Mi2+2.*Mi2*sqrt(Mi1*Mi1+pb*pb) );

  G4ThreeVector LPm1, LPm2;
  G4ThreeVector LPm1CM, LPm2CM;

  AGUniform  genU( BP->cost1, BP->cost2 );
  AGForward1 genF( BP->cost1, BP->cost2 );

  G4bool status1 = false;
  if( confMan->GetCMAngle()==1 ){
    status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
			       LPm1, LPm2, LPm1CM, LPm2CM, genU );
  }
  else if( confMan->GetCMAngle()==2 ){
    status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
			       LPm1, LPm2, LPm1CM, LPm2CM, genF );
  }

  G4int n1=0;
  while(!status1){
    if(++n1>MaxTry){
      G4Exception("Reaction under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    while( !(Ei>Mm1+Mm2) ) Mm2=BreitWigner( Mass, Width );

    if( confMan->GetCMAngle()==1 ){
      status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
				 LPm1, LPm2, LPm1CM, LPm2CM, genU );
    }
    else if( confMan->GetCMAngle()==2 ){
      status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
				 LPm1, LPm2, LPm1CM, LPm2CM, genF );
    }
  }
  
  //D*- decay 
  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;
  AGUniform gen2;
  G4bool status2 = Decay2Body( Mm1, Mf1, Mf2, LPm1, LPf1, LPf2, LPf1CM, LPf2CM, gen2 );
  if(!status2){
    G4Exception("Reaction under Threshold");
  }

  G4ThreeVector gPos(LPos), gPf1(LPf1), gPf2(LPf2);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  //gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  gPf2.rotateX(90.0*degree); gPf2.rotateZ( BP->TargetRotAngH );
  LBeamDir.rotateX(90.0*degree); LBeamDir.rotateZ( BP->TargetRotAngH );
  
  //Beam
  pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM_->GeneratePrimaryVertex( anEvent );

  ///////D*- decay by EvtGen
  EvtParticle* parent(0);
  ////cbar mesons 
  static EvtId DS = EvtPDL::getId(std::string("D*-"));  
  //static EvtId DS = EvtPDL::getId(std::string("anti-D*0"));
  //static EvtId DS = EvtPDL::getId(std::string("anti-D0"));
  //static EvtId DS = EvtPDL::getId(std::string("D-"));

  ////c mesons 
  //static EvtId DS = EvtPDL::getId(std::string("D*+"));  
  //static EvtId DS = EvtPDL::getId(std::string("D*0"));
  //static EvtId DS = EvtPDL::getId(std::string("D0"));
  //static EvtId DS = EvtPDL::getId(std::string("D+"));

  G4LorentzVector LvDS0;
  LvDS0.setVect(LPm1/GeV);
  LvDS0.setE(sqrt(Mm1*Mm1+LPm1.mag2())/GeV);

  //G4cout<< LvDS0.mag() << " " 
  //<< LvDS0.vect().x() << " " << LvDS0.vect().y() << " " << LvDS0.vect().z() << G4endl;

  EvtVector4R pInit( LvDS0.e(),
   		     LvDS0.vect().x(), 
  		     LvDS0.vect().y(), 
  		     LvDS0.vect().z() );

  parent = EvtParticleFactory::particleFactory(DS, pInit); 
  //parent->setVectorSpinDensity();
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
    // G4double life1=q->getLifetime();
    // G4cout << "Life1=" <<life1 << G4endl;
    // q->setLifetime(0.);    
    // G4double life2=q->getLifetime();
    // G4cout << "Life2=" <<life2 << G4endl;
    // G4cout << "*************" << G4endl;
    // G4cout << q->getNDaug() << G4endl;
    // G4cout << "mass=" << q->mass() << G4endl;
    // G4cout << q->getName() << G4endl;
    // std::cout<< "Px= " << P4Lab.get(1) 
    // 	     << " Py= "<< P4Lab.get(2) 
    // 	     << " Pz= "<< P4Lab.get(3) << std::endl;

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
    pGen->SetPrimaryInformation( LPos, Mm1, Mm2, 
                                 LPm1, LPm2, LPm1CM, LPm2CM, pb*LBeamDir, 0. );
  }

  parent->deleteTree();

  //YcDStarMesonFromPiMinusProton2
  //21: pi- p -> Yc+ D*- : FS All state 
}


////////////////////
//22: pi- p -> Yc+ D** : FS All state
void PGReaction1::YcDStarMesonFromPiMinusProton3( G4Event *anEvent, EvtGen *evtGenerator )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  //

  ConfMan *confMan = ConfMan::GetConfManager();

  BeamParam *BP = pGen->BP_;
  G4double Mass=BP->RecMass;
  G4double Width=BP->RecWidth;
 
  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  ////D** mass
  //G4double Mf1=2.4230*GeV;//D_1- mass
  //G4double Mf1=2.4223*GeV;//D_10 mass
  //G4double Mf1=2.4450*GeV;//D'_1- mass
  //G4double Mf1=2.4450*GeV;//D'_10 mass
  //G4double Mf1=2.4601*GeV;//D_2*- mass
  //G4double Mf1=2.4611*GeV;//D_2*0 mass
  //G4double Mf1=2.5800*GeV;//D(2S)0 mass
  //G4double Mf1=2.5800*GeV;//D(2S)- mass
  //G4double Mf1=2.6400*GeV;//D*(2S)0 mass
  //G4double Mf1=2.6400*GeV;//D*(2S)- mass

  ////cbar mesons
  G4int nev=anEvent->GetEventID();
  G4double Mf1, Wid; 
  static EvtId DS;
  switch( nev%6 ){
  case 0:
    Mf1 = 2.42300*GeV;//D_1- mass
    Wid = 0.0200*GeV;
    DS = EvtPDL::getId(std::string("D_1-"));  
    break;
  case 1:
    Mf1 = 2.42230*GeV;//D_10 mass
    Wid = 0.0204*GeV;
    DS = EvtPDL::getId(std::string("anti-D_10"));
    break;
  case 2:
    Mf1 = 2.4450*GeV;//D'_1- mass
    Wid = 0.250*GeV;
    DS = EvtPDL::getId(std::string("D'_1-"));
    break;
  case 3:
    Mf1 = 2.4450*GeV;//D'_10 mass
    Wid = 0.250*GeV;
    DS = EvtPDL::getId(std::string("anti-D'_10"));
    break;
  case 4:
    Mf1 = 2.46010*GeV;//D_2*- mass
    Wid = 0.037*GeV;
    DS = EvtPDL::getId(std::string("D_2*-"));
    break;
  case 5:
    Mf1 = 2.46110*GeV;//D_2*0 mass
    Wid = 0.043*GeV;
    DS = EvtPDL::getId(std::string("anti-D_2*0"));
    break;
  // case 6:
  //   Mf1=2.580*GeV;//D(2S)0 mass
  //   DS = EvtPDL::getId(std::string("anti-D(2S)0"));
  //   break;
  // case 7:
  //   Mf1=2.580*GeV;//D(2S)- mass
  //   DS = EvtPDL::getId(std::string("D(2S)-"));
  //   break;
  // case 8:
  //   Mf1=2.640*GeV;//D*(2S)0 mass
  //   DS = EvtPDL::getId(std::string("anti-D*(2S)0"));
  //   break;
  // case 9:
  //   Mf1=2.640*GeV;//D*(2S)- mass
  //   DS = EvtPDL::getId(std::string("D*(2S)-"));
  //   break;
  default:
    Mf1 = 2.01027*GeV; 
    Wid = 0.0*GeV;
    DS = EvtPDL::getId(std::string("D*-"));
    break;
  }

  //static EvtId DS = EvtPDL::getId(std::string("D_1-"));  
  //static EvtId DS = EvtPDL::getId(std::string("anti-D_10"));
  //static EvtId DS = EvtPDL::getId(std::string("D'_1-"));
  //static EvtId DS = EvtPDL::getId(std::string("anti-D'_10"));
  //static EvtId DS = EvtPDL::getId(std::string("D_2*-"));
  //static EvtId DS = EvtPDL::getId(std::string("anti-D_2*0"));
  //static EvtId DS = EvtPDL::getId(std::string("anti-D(2S)0"));
  //static EvtId DS = EvtPDL::getId(std::string("D(2S)-"));
  //static EvtId DS = EvtPDL::getId(std::string("anti-D*(2S)0"));
  //static EvtId DS = EvtPDL::getId(std::string("D*(2S)-"));
  //static EvtId DS = EvtPDL::getId(std::string("D_1-"));  

  ////c mesons
  // G4int nev=anEvent->GetEventID();
  // G4double Mf1; 
  // static EvtId DS;
  // switch( nev%10 ){
  // case 0:
  //   Mf1=2.42300*GeV;//D_1+ mass
  //   DS = EvtPDL::getId(std::string("D_1+"));  
  //   break;
  // case 1:
  //   Mf1=2.42230*GeV;//D_10 mass
  //   DS = EvtPDL::getId(std::string("D_10"));
  //   break;
  // case 2:
  //   Mf1=2.4450*GeV;//D'_1+ mass
  //   DS = EvtPDL::getId(std::string("D'_1+"));
  //   break;
  // case 3:
  //   Mf1=2.4450*GeV;//D'_10 mass
  //   DS = EvtPDL::getId(std::string("D'_10"));
  //   break;
  // case 4:
  //   Mf1=2.46010*GeV;//D_2*+ mass
  //   DS = EvtPDL::getId(std::string("D_2*+"));
  //   break;
  // case 5:
  //   Mf1=2.46110*GeV;//D_2*0 mass
  //   DS = EvtPDL::getId(std::string("D_2*0"));
  //   break;
  // case 6:
  //   Mf1=2.580*GeV;//D(2S)0 mass
  //   DS = EvtPDL::getId(std::string("D(2S)0"));
  //   break;
  // case 7:
  //   Mf1=2.580*GeV;//D(2S)+ mass
  //   DS = EvtPDL::getId(std::string("D(2S)+"));
  //   break;
  // case 8:
  //   Mf1=2.640*GeV;//D*(2S)0 mass
  //   DS = EvtPDL::getId(std::string("D*(2S)0"));
  //   break;
  // case 9:
  //   Mf1=2.640*GeV;//D*(2S)+ mass
  //   DS = EvtPDL::getId(std::string("D*(2S)+"));
  //   break;
  // default:
  //   Mf1=2.01027*GeV; 
  //   DS = EvtPDL::getId(std::string("D*+"));
  //   break;
  // }

  //static EvtId DS = EvtPDL::getId(std::string("D_10"));
  //static EvtId DS = EvtPDL::getId(std::string("D'_1+"));
  //static EvtId DS = EvtPDL::getId(std::string("D'_10"));
  //static EvtId DS = EvtPDL::getId(std::string("D_2*+"));
  //static EvtId DS = EvtPDL::getId(std::string("D_2*0"));
  //static EvtId DS = EvtPDL::getId(std::string("D(2S)0"));
  //static EvtId DS = EvtPDL::getId(std::string("D(2S)+"));
  //static EvtId DS = EvtPDL::getId(std::string("D*(2S)0"));
  //static EvtId DS = EvtPDL::getId(std::string("D*(2S)+"));

  G4double Mf2=BreitWigner( Mass, Width );
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
  
  // G4cout<< "****1" << G4endl;
  // G4cout<< MF << G4endl;

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
      Mf2=BreitWigner( Mass, Width );

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

  // G4cout<< "****2" << G4endl;
  // G4cout<< MF << " " << LPf1.mag() << G4endl;

  G4ThreeVector gPos(LPos), gPf1(LPf1), gPf2(LPf2);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  gPf2.rotateX(90.0*degree); gPf2.rotateZ( BP->TargetRotAngH );
  LBeamDir.rotateX(90.0*degree); LBeamDir.rotateZ( BP->TargetRotAngH );
  
  //Beam
  pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM_->GeneratePrimaryVertex( anEvent );

  ///////D*- decay by EvtGen
  EvtParticle* parent(0);

  G4LorentzVector LvDS0;
  LvDS0.setVect(LPf1/GeV);
  LvDS0.setE(sqrt(MF*MF+LPf1.mag2())/GeV);

  // G4cout<< LvDS0.mag() 
  // 	<< " " << LvDS0.vect().x() << " " << LvDS0.vect().y() << " " << LvDS0.vect().z() << G4endl;

  EvtVector4R pInit( LvDS0.e(),
   		     LvDS0.vect().x(), 
  		     LvDS0.vect().y(), 
  		     LvDS0.vect().z() );

  // G4cout<< "****" << G4endl;
  // G4cout<< LvDS0.e() << " " << LvDS0.mag() << G4endl;

  parent = EvtParticleFactory::particleFactory(DS, pInit); 
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
    EvtVector4R P4Lab=p->getP4Lab();

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
    // G4double life1=q->getLifetime();
    // G4cout << "Life1=" <<life1 << G4endl;
    // q->setLifetime(0.);    
    // G4double life2=q->getLifetime();
    // G4cout << "Life2=" <<life2 << G4endl;
    // G4cout << "*************" << G4endl;
    // G4cout << q->getNDaug() << G4endl;
    // G4cout << "mass=" << q->mass() << G4endl;
    // G4cout << q->getName() << G4endl;
    // std::cout<< "Px= " << P4Lab.get(1) 
    // 	     << " Py= "<< P4Lab.get(2) 
    // 	     << " Pz= "<< P4Lab.get(3) << std::endl;

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
    pGen->SetPrimaryInformation( LPos, Mf1, Mf2, 
                                 LPf1, LPf2, LPf1CM, LPf2CM, pb*LBeamDir, 0. );
  }

  parent->deleteTree();

  //YcDStarMesonFromPiMinusProton3
  //22: pi- p -> Yc+ D** : FS All state
}


////////////////////
//23: pi- p -> Yc+ D*- : FS K+ pi- pi+ pi- pi- 
void PGReaction1::YcDStarMesonFromPiMinusProton4( G4Event *anEvent )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  //

  ConfMan *confMan = ConfMan::GetConfManager();

  BeamParam *BP = pGen->BP_;
  G4double Mass=BP->RecMass;
  G4double Width=BP->RecWidth;
 
  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  G4double Mm1=2.01022*GeV;//D*- mass
  G4double Mm2=BreitWigner( Mass, Width );

  G4double Mf1=G4DMesonZero::Definition()->GetPDGMass();
  G4double Mf2=G4PionMinus::Definition()->GetPDGMass();

  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4double Ei=sqrt( Mi1*Mi1+Mi2*Mi2+2.*Mi2*sqrt(Mi1*Mi1+pb*pb) );

  G4ThreeVector LPm1, LPm2;
  G4ThreeVector LPm1CM, LPm2CM;

  AGUniform  genU( BP->cost1, BP->cost2 );
  AGForward1 genF( BP->cost1, BP->cost2 );

  G4bool status1 = false;
  if( confMan->GetCMAngle()==1 ){
    status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
			       LPm1, LPm2, LPm1CM, LPm2CM, genU );
  }
  else if( confMan->GetCMAngle()==2 ){
    status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
			       LPm1, LPm2, LPm1CM, LPm2CM, genF );
  }

  G4int n1=0;
  while(!status1){
    if(++n1>MaxTry){
      G4Exception("Reaction under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    while( !(Ei>Mm1+Mm2) ) Mm2=BreitWigner( Mass, Width );

    if( confMan->GetCMAngle()==1 ){
      status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
				 LPm1, LPm2, LPm1CM, LPm2CM, genU );
    }
    else if( confMan->GetCMAngle()==2 ){
      status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
				 LPm1, LPm2, LPm1CM, LPm2CM, genF );
    }
  }
  
  //D*- decay 
  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;
  AGUniform gen2;
  G4bool status2 = Decay2Body( Mm1, Mf1, Mf2, LPm1, LPf1, LPf2, LPf1CM, LPf2CM, gen2 );
  if(!status2){
    G4Exception("Reaction under Threshold");
  }

  G4ThreeVector gPos(LPos), gPf1(LPf1), gPf2(LPf2);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  gPf2.rotateX(90.0*degree); gPf2.rotateZ( BP->TargetRotAngH );
  LBeamDir.rotateX(90.0*degree); LBeamDir.rotateZ( BP->TargetRotAngH );

  pGen->gunDMesonZ2_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunDMesonZ2_->SetParticleMomentum( gPf1 );
  pGen->gunDMesonZ2_->GeneratePrimaryVertex( anEvent );
  
  pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM1_->SetParticleMomentum( gPf2 );
  pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );

  // pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  // pGen->gunPiM_->SetParticleMomentum( gPf2 );
  // pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
  
  //Beam
  pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM1_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );

  if( pGen->anaMan_ ){
    LBeamDir.rotateZ( -1.*(BP->TargetRotAngH) ); 
    LBeamDir.rotateX(-90.0*degree);
    gPf2.rotateZ( -1.*(BP->TargetRotAngH) ); 
    gPf2.rotateX(-90.0*degree);
    pGen->SetPrimaryInformation( LPos, Mm1, Mm2, 
                                 LPm1, LPm2, LPm1CM, LPm2CM, pb*LBeamDir, gPf2 );
  }
  //YcDStarMesonFromPiMinusProton4
  //23: pi- p -> Yc+ D*- : FS K+ pi- pi+ pi- pi-
}


////////////////////
//24: pi- p -> Yc+ D*- : FS K+ pi- pi-, Yc+ decay pi+ pi-
void PGReaction1::YcDStarMesonFromPiMinusProton5( G4Event *anEvent )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  //

  ConfMan *confMan = ConfMan::GetConfManager();

  BeamParam *BP = pGen->BP_;
  G4double Mass=BP->RecMass;
  G4double Width=BP->RecWidth;
 
  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  G4double Mm1=2.01022*GeV;//D*- mass
  G4double Mm2=BreitWigner( Mass, Width );

  G4double Mf1=G4DMesonZero::Definition()->GetPDGMass();
  G4double Mf2=G4PionMinus::Definition()->GetPDGMass();

  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4double Ei=sqrt( Mi1*Mi1+Mi2*Mi2+2.*Mi2*sqrt(Mi1*Mi1+pb*pb) );

  G4ThreeVector LPm1, LPm2;
  G4ThreeVector LPm1CM, LPm2CM;

  AGUniform  genU( BP->cost1, BP->cost2 );
  AGForward1 genF( BP->cost1, BP->cost2 );

  G4bool status1 = false;
  if( confMan->GetCMAngle()==1 ){
    status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
			       LPm1, LPm2, LPm1CM, LPm2CM, genU );
  }
  else if( confMan->GetCMAngle()==2 ){
    status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
			       LPm1, LPm2, LPm1CM, LPm2CM, genF );
  }

  G4int n1=0;
  while(!status1){
    if(++n1>MaxTry){
      G4Exception("Reaction under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    while( !(Ei>Mm1+Mm2) ) Mm2=BreitWigner( Mass, Width );

    if( confMan->GetCMAngle()==1 ){
      status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
				 LPm1, LPm2, LPm1CM, LPm2CM, genU );
    }
    else if( confMan->GetCMAngle()==2 ){
      status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
				 LPm1, LPm2, LPm1CM, LPm2CM, genF );
    }
  }
  
  //D*- decay 
  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;
  AGUniform gen2;
  G4bool status2 = Decay2Body( Mm1, Mf1, Mf2, LPm1, LPf1, LPf2, LPf1CM, LPf2CM, gen2 );
  if(!status2){
    G4Exception("Reaction under Threshold");
  }

  G4ThreeVector gPos(LPos), gPf1(LPf1), gPf2(LPf2);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  gPf2.rotateX(90.0*degree); gPf2.rotateZ( BP->TargetRotAngH );
  LBeamDir.rotateX(90.0*degree); LBeamDir.rotateZ( BP->TargetRotAngH );

  pGen->gunDMesonZ1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunDMesonZ1_->SetParticleMomentum( gPf1 );
  pGen->gunDMesonZ1_->GeneratePrimaryVertex( anEvent );
  
  pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM1_->SetParticleMomentum( gPf2 );
  pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );
  // pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  // pGen->gunPiM_->SetParticleMomentum( gPf2 );
  // pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
  
  //Beam
  pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM1_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );
  
  //Yc*+ decay 
  G4double Mfd1_1=2.4540*GeV;
  G4double Mfd1_2=G4PionPlus::Definition()->GetPDGMass();
  //G4double Mfd1_2=G4PionZero::Definition()->GetPDGMass();
  G4ThreeVector LPfd1_1, LPfd1_2;
  G4ThreeVector LPfd1_1CM, LPfd1_2CM;

  G4double Mfd2_1=2.2865*GeV;
  G4double Mfd2_2=G4PionPlus::Definition()->GetPDGMass();
  //G4double Mfd2_2=G4PionZero::Definition()->GetPDGMass();
  G4ThreeVector LPfd2_1, LPfd2_2;
  G4ThreeVector LPfd2_1CM, LPfd2_2CM;
  
  if( Mm2>2.6 ){
    AGUniform gen3;
    Decay2Body( Mm2, Mfd1_1, Mfd1_2, LPm2, 
		LPfd1_1, LPfd1_2, LPfd1_1CM, LPfd1_2CM, gen3 );

    AGUniform gen4;
    Decay2Body( Mfd1_1, Mfd2_1, Mfd2_2, LPfd1_1, 
		LPfd2_1, LPfd2_2, LPfd2_1CM, LPfd2_2CM, gen4 );
  }

  G4ThreeVector gPfd1_1(LPfd1_1), gPfd1_2(LPfd1_2);
  G4ThreeVector gPfd2_1(LPfd2_1), gPfd2_2(LPfd2_2);
  gPfd1_1.rotateX(90.0*degree); gPfd1_1.rotateZ( BP->TargetRotAngH );
  gPfd1_2.rotateX(90.0*degree); gPfd1_2.rotateZ( BP->TargetRotAngH );
  gPfd2_1.rotateX(90.0*degree); gPfd2_1.rotateZ( BP->TargetRotAngH );
  gPfd2_2.rotateX(90.0*degree); gPfd2_2.rotateZ( BP->TargetRotAngH );
  
  int nev=anEvent->GetEventID();
  switch( nev%2 ){
  case 0:
    pGen->gunPiP1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiP1_->SetParticleMomentum( gPfd1_2 );
    pGen->gunPiP1_->GeneratePrimaryVertex( anEvent );
    
    pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM1_->SetParticleMomentum( gPfd2_2 );
    pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );
    break;
  case 1:
    pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM1_->SetParticleMomentum( gPfd1_2 );
    pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );
    
    pGen->gunPiP1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiP1_->SetParticleMomentum( gPfd2_2 );
    pGen->gunPiP1_->GeneratePrimaryVertex( anEvent );
    break;
    default:
      break;
  }

  // int nev=anEvent->GetEventID();
  // switch( nev%2 ){
  // case 0:
  //   pGen->gunPiP_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunPiP_->SetParticleMomentum( gPfd1_2 );
  //   pGen->gunPiP_->GeneratePrimaryVertex( anEvent );
    
  //   pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunPiM_->SetParticleMomentum( gPfd2_2 );
  //   pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
  //   break;
  // case 1:
  //   pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunPiM_->SetParticleMomentum( gPfd1_2 );
  //   pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
    
  //   pGen->gunPiP_->SetParticlePosition( gPos+BP->gTargetPos );
  //   pGen->gunPiP_->SetParticleMomentum( gPfd2_2 );
  //   pGen->gunPiP_->GeneratePrimaryVertex( anEvent );
  //   break;
  //   default:
  //     break;
  // }

  //   switch( nev%2 ){
  //   case 0:
  //     pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
  //     pGen->gunPi0_->SetParticleMomentum( gPfd1_2 );
  //     pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
      
  //     pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
  //     pGen->gunPi0_->SetParticleMomentum( gPfd2_2 );
  //     pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  //     break;
  //   case 1:
  //     pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
  //     pGen->gunPi0_->SetParticleMomentum( gPfd1_2 );
  //     pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
      
  //     pGen->gunPi0_->SetParticlePosition( gPos+BP->gTargetPos );
  //     pGen->gunPi0_->SetParticleMomentum( gPfd2_2 );
  //     pGen->gunPi0_->GeneratePrimaryVertex( anEvent );
  //     break;
  //   default:
  //     break;
  //   }
  // }

  if( pGen->anaMan_ ){
    LBeamDir.rotateZ( -1.*(BP->TargetRotAngH) ); 
    LBeamDir.rotateX(-90.0*degree);
    gPf2.rotateZ( -1.*(BP->TargetRotAngH) ); 
    gPf2.rotateX(-90.0*degree);
    // pGen->SetPrimaryInformation( LPos, Mfd1_1, Mfd1_2, 
    //                              LPfd1_1, LPfd1_2, LPfd1_1CM, LPfd1_2CM, pb*LBeamDir, gPf2 );
    pGen->SetPrimaryInformation( LPos, Mm1, Mm2, 
                                 LPm1, LPm2, LPm1CM, LPm2CM, pb*LBeamDir, gPf2 );
  }
  //YcDStarMesonFromPiMinusProton5
  //24: pi- p -> Yc+ D*- : FS K+ pi- pi-, Yc+ decay pi+ pi-
}


////////////////////
//25: pi- p -> Yc+ D*- : FS K+ pi- pi+ pi- pi-, Yc+ decay pi+ pi-
void PGReaction1::YcDStarMesonFromPiMinusProton6( G4Event *anEvent )
{
  // Caution !!
  // exp(bt) OR uniform angular distribution in CM 
  // Not real distribution
  //
  
  ConfMan *confMan = ConfMan::GetConfManager();

  BeamParam *BP = pGen->BP_;
  G4double Mass=BP->RecMass;
  G4double Width=BP->RecWidth;
 
  G4double Mi1=G4PionMinus::Definition()->GetPDGMass();
  G4double Mi2=G4Proton::Definition()->GetPDGMass();

  G4double Mm1=2.01022*GeV;//D*- mass
  G4double Mm2=BreitWigner( Mass, Width );

  G4double Mf1=G4DMesonZero::Definition()->GetPDGMass();
  G4double Mf2=G4PionMinus::Definition()->GetPDGMass();

  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

  G4double Ei=sqrt( Mi1*Mi1+Mi2*Mi2+2.*Mi2*sqrt(Mi1*Mi1+pb*pb) );

  G4ThreeVector LPm1, LPm2;
  G4ThreeVector LPm1CM, LPm2CM;

  AGUniform  genU( BP->cost1, BP->cost2 );
  AGForward1 genF( BP->cost1, BP->cost2 );

  G4bool status1 = false;
  if( confMan->GetCMAngle()==1 ){
    status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
			       LPm1, LPm2, LPm1CM, LPm2CM, genU );
  }
  else if( confMan->GetCMAngle()==2 ){
    status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
			       LPm1, LPm2, LPm1CM, LPm2CM, genF );
  }

  G4int n1=0;
  while(!status1){
    if(++n1>MaxTry){
      G4Exception("Reaction under Threshold");
    }
    pb = BP->p0; if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;

    while( !(Ei>Mm1+Mm2) ) Mm2=BreitWigner( Mass, Width );

    if( confMan->GetCMAngle()==1 ){
      status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
				 LPm1, LPm2, LPm1CM, LPm2CM, genU );
    }
    else if( confMan->GetCMAngle()==2 ){
      status1 = Scattering2Body( Mi1, Mi2, Mm1, Mm2, pb*LBeamDir,
				 LPm1, LPm2, LPm1CM, LPm2CM, genF );
    }
  }
  
  //D*- decay 
  G4ThreeVector LPf1, LPf2;
  G4ThreeVector LPf1CM, LPf2CM;
  AGUniform gen2;
  G4bool status2 = Decay2Body( Mm1, Mf1, Mf2, LPm1, LPf1, LPf2, LPf1CM, LPf2CM, gen2 );
  if(!status2){
    G4Exception("Reaction under Threshold");
  }

  G4ThreeVector gPos(LPos), gPf1(LPf1), gPf2(LPf2);
  gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
  gPf1.rotateX(90.0*degree); gPf1.rotateZ( BP->TargetRotAngH );
  gPf2.rotateX(90.0*degree); gPf2.rotateZ( BP->TargetRotAngH );
  LBeamDir.rotateX(90.0*degree); LBeamDir.rotateZ( BP->TargetRotAngH );

  pGen->gunDMesonZ2_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunDMesonZ2_->SetParticleMomentum( gPf1 );
  pGen->gunDMesonZ2_->GeneratePrimaryVertex( anEvent );
  
  pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM1_->SetParticleMomentum( gPf2 );
  pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );
  // pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
  // pGen->gunPiM_->SetParticleMomentum( gPf2 );
  // pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
  
  //Beam
  pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
  pGen->gunPiM1_->SetParticleMomentum( -1.*pb*LBeamDir );
  pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );

  //Yc*+ decay 
  G4double Mfd1_1=2.4540*GeV;
  G4double Mfd1_2=G4PionPlus::Definition()->GetPDGMass();
  G4ThreeVector LPfd1_1, LPfd1_2;
  G4ThreeVector LPfd1_1CM, LPfd1_2CM;

  G4double Mfd2_1=2.2865*GeV;
  G4double Mfd2_2=G4PionPlus::Definition()->GetPDGMass();
  G4ThreeVector LPfd2_1, LPfd2_2;
  G4ThreeVector LPfd2_1CM, LPfd2_2CM;
  
  if( Mm2>2.6 ){
    AGUniform gen3;
    Decay2Body( Mm2, Mfd1_1, Mfd1_2, LPm2, 
		LPfd1_1, LPfd1_2, LPfd1_1CM, LPfd1_2CM, gen3 );

    AGUniform gen4;
    Decay2Body( Mfd1_1, Mfd2_1, Mfd2_2, LPfd1_1, 
		LPfd2_1, LPfd2_2, LPfd2_1CM, LPfd2_2CM, gen4 );
  }
   
  G4ThreeVector gPfd1_1(LPfd1_1), gPfd1_2(LPfd1_2);
  G4ThreeVector gPfd2_1(LPfd2_1), gPfd2_2(LPfd2_2);
  gPfd1_1.rotateX(90.0*degree); gPfd1_1.rotateZ( BP->TargetRotAngH );
  gPfd1_2.rotateX(90.0*degree); gPfd1_2.rotateZ( BP->TargetRotAngH );
  gPfd2_1.rotateX(90.0*degree); gPfd2_1.rotateZ( BP->TargetRotAngH );
  gPfd2_2.rotateX(90.0*degree); gPfd2_2.rotateZ( BP->TargetRotAngH );
   
  int nev=anEvent->GetEventID();
  switch( nev%2 ){
  case 0:
    pGen->gunPiP1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiP1_->SetParticleMomentum( gPfd1_2 );
    pGen->gunPiP1_->GeneratePrimaryVertex( anEvent );
    
    pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM1_->SetParticleMomentum( gPfd2_2 );
    pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );
    break;
  case 1:
    pGen->gunPiM1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiM1_->SetParticleMomentum( gPfd1_2 );
    pGen->gunPiM1_->GeneratePrimaryVertex( anEvent );
      
    pGen->gunPiP1_->SetParticlePosition( gPos+BP->gTargetPos );
    pGen->gunPiP1_->SetParticleMomentum( gPfd2_2 );
    pGen->gunPiP1_->GeneratePrimaryVertex( anEvent );
    break;
  default:
    break;
  }

  if( pGen->anaMan_ ){
    LBeamDir.rotateZ( -1.*(BP->TargetRotAngH) ); 
    LBeamDir.rotateX(-90.0*degree);
    gPf2.rotateZ( -1.*(BP->TargetRotAngH) ); 
    gPf2.rotateX(-90.0*degree);
    // pGen->SetPrimaryInformation( LPos, Mfd1_1, Mfd1_2, 
    //                              LPfd1_1, LPfd1_2, LPfd1_1CM, LPfd1_2CM, pb*LBeamDir, gPf2 );
    pGen->SetPrimaryInformation( LPos, Mm1, Mm2, 
                                 LPm1, LPm2, LPm1CM, LPm2CM, pb*LBeamDir, gPf2 );
  }
  //YcDStarMesonFromPiMinusProton6
  //25: pi- p -> Yc+ D*- : FS K+ pi- pi+ pi- pi-, Yc+ decay pi+ pi-
}
