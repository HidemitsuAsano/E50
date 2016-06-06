/*
  PGReaction4.cc 

  2016/1  K.Shirotori
*/

////JAM Data

#include "PGReaction4.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "Analysis.hh"
#include "PrimaryInfo.hh"
#include "GeneratorHelper.hh"
#include "AngDisGenerator.hh"
#include "ReadJAMData.hh"
#include "JAMData.hh"
#include "ParticleID.hh"
#include "ConfMan.hh"

#include "G4ios.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTypes.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "G4UImanager.hh"

const G4int MaxTry=100;

const G4double Deg2Rad = acos(-1.)/180.;
const G4double Rad2Deg = 180./acos(-1.);

//D0 Mass cut -30 -> +30 MeV
const G4double MinDMMass = 1835.;
const G4double MaxDMMass = 1895.;

//Q-value cut -5 -> +10 MeV
const G4double MinQV =  0.;
const G4double MaxQV = 15.;

////////////////////////////////////////Background from JAM
////////////////////
void PGReaction4::
BackgroundFromPiMinusProtonJAM1( G4Event *anEvent, ReadJAMData *readJAMData )
{
  // Event generated from JAM data
  // 

  ConfMan *confMan = ConfMan::GetConfManager();

  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4PionMinus::Definition()->GetPDGMass();

  BeamParam *BP = pGen->BP_;
  G4ThreeVector LPos = G4ThreeVector( BP->x0, BP->y0, BP->z0 )
    + GaussPosition( BP->dx, BP->dy, BP->dztarg );
  G4ThreeVector LBeamDir = 
    GaussDirectionInUV( BP->u0, BP->v0, BP->du, BP->dv );
  G4double pb = BP->p0;
  if(BP->dp!=0.0) pb += (G4UniformRand()-0.5)*BP->dp;
  
  G4int nev=anEvent->GetEventID();
  const JAMDataContainer &cont=readJAMData->GetJAMData();

  // std::cout<< nev <<std::endl;
  // std::cout<< cont.size() <<std::endl;

  G4int offset=0;
  if( nev<cont.size() ){

    G4int nkp=0, nkm=0;
    G4int npip=0, npim=0;
    G4int nproton=0, nprotonb=0;
    G4int npiz=0;
    G4int nneutron=0, nneutronb=0;
    G4int nkz=0, nkzs=0, nkzl=0;
    G4int ngamma=0;
    G4int nep=0, nem=0;
    G4int nmup=0, nmum=0;
    
    JAMData *hit=cont[nev+offset];
    G4int nhpid=hit->GetPidSize();
    for( G4int j=0; j<nhpid; j++ ){
      G4int pid=hit->GetPid(j);
      if( pid == IdKaonP )   nkp++;
      if( pid == IdKaonM )   nkm++;
      if( pid == IdPionP )   npip++;
      if( pid == IdPionM )   npim++;
      if( pid == IdProton )  nproton++;
      if( pid == IdProtonb ) nprotonb++;
      if( pid == IdPionZ )   npiz++;
      if( pid == IdNeutron ) nneutron++;
      if( pid == IdNeutronb) nneutronb++;
      if( pid == IdKaonZ )   nkz++;
      if( pid == IdKaonZS )  nkzs++;
      if( pid == IdKaonZL )  nkzl++;
      if( pid == IdGamma )   ngamma++;
      if( pid == IdPosi )    nep++;
      if( pid == IdElec )    nem++;
      if( pid == IdMuonP )   nmup++;
      if( pid == IdMuonM )   nmum++;
    }

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
    
    for( G4int j=0; j<nhpid; j++ ){
      G4int pid=hit->GetPid(j);
      if( pid == IdKaonP ){
	gPfKP[ikp].setX(hit->GetMomX(j));
	gPfKP[ikp].setY(hit->GetMomY(j));
	gPfKP[ikp].setZ(hit->GetMomZ(j));
	ikp++;
      }
      if( pid == IdKaonM ){
	gPfKM[ikm].setX(hit->GetMomX(j));
	gPfKM[ikm].setY(hit->GetMomY(j));
	gPfKM[ikm].setZ(hit->GetMomZ(j));
	ikm++;
      }
      if( pid == IdPionP ){
	gPfPiP[ipip].setX(hit->GetMomX(j));
	gPfPiP[ipip].setY(hit->GetMomY(j));
	gPfPiP[ipip].setZ(hit->GetMomZ(j));
	ipip++;
      }
      if( pid == IdPionM ){
	gPfPiM[ipim].setX(hit->GetMomX(j));
	gPfPiM[ipim].setY(hit->GetMomY(j));
	gPfPiM[ipim].setZ(hit->GetMomZ(j));
	ipim++;
      }
      if( pid == IdProton ){
	gPfP[ip].setX(hit->GetMomX(j));
	gPfP[ip].setY(hit->GetMomY(j));
	gPfP[ip].setZ(hit->GetMomZ(j));
	ip++;
      }
      if( pid == IdProtonb ){
	gPfPb[ipb].setX(hit->GetMomX(j));
	gPfPb[ipb].setY(hit->GetMomY(j));
	gPfPb[ipb].setZ(hit->GetMomZ(j));
	ipb++;
      }
      if( pid == IdPionZ ){
	gPfPiZ[ipiz].setX(hit->GetMomX(j));
	gPfPiZ[ipiz].setY(hit->GetMomY(j));
	gPfPiZ[ipiz].setZ(hit->GetMomZ(j));
	ipiz++;
      }
      if( pid == IdNeutron ){
	gPfN[in].setX(hit->GetMomX(j));
	gPfN[in].setY(hit->GetMomY(j));
	gPfN[in].setZ(hit->GetMomZ(j));
	in++;
      }
      if( pid == IdNeutronb ){
	gPfNb[inb].setX(hit->GetMomX(j));
	gPfNb[inb].setY(hit->GetMomY(j));
	gPfNb[inb].setZ(hit->GetMomZ(j));
	inb++;
      }
      if( pid == IdKaonZ ){
	gPfKZ[ikz].setX(hit->GetMomX(j));
	gPfKZ[ikz].setY(hit->GetMomY(j));
	gPfKZ[ikz].setZ(hit->GetMomZ(j));
	ikz++;
      }
      if( pid == IdKaonZS ){
	gPfKZS[ikzs].setX(hit->GetMomX(j));
	gPfKZS[ikzs].setY(hit->GetMomY(j));
	gPfKZS[ikzs].setZ(hit->GetMomZ(j));
	ikzs++;
      }
      if( pid == IdKaonZL ){
	gPfKZL[ikzl].setX(hit->GetMomX(j));
	gPfKZL[ikzl].setY(hit->GetMomY(j));
	gPfKZL[ikzl].setZ(hit->GetMomZ(j));
	ikzl++;
      }
      if( pid == IdGamma ){
	gPfG[igamma].setX(hit->GetMomX(j));
	gPfG[igamma].setY(hit->GetMomY(j));
	gPfG[igamma].setZ(hit->GetMomZ(j));
	igamma++;
      }
      if( pid == IdPosi ){
	gPfEP[iep].setX(hit->GetMomX(j));
	gPfEP[iep].setY(hit->GetMomY(j));
	gPfEP[iep].setZ(hit->GetMomZ(j));
	iep++;
      }
      if( pid == IdElec ){
	gPfEM[iem].setX(hit->GetMomX(j));
	gPfEM[iem].setY(hit->GetMomY(j));
	gPfEM[iem].setZ(hit->GetMomZ(j));
	iem++;
      }
      if( pid == IdMuonP ){
	gPfMuP[imup].setX(hit->GetMomX(j));
	gPfMuP[imup].setY(hit->GetMomY(j));
	gPfMuP[imup].setZ(hit->GetMomZ(j));
	imup++;
      }
      if( pid == IdMuonM ){
	gPfMuM[imum].setX(hit->GetMomX(j));
	gPfMuM[imum].setY(hit->GetMomY(j));
	gPfMuM[imum].setZ(hit->GetMomZ(j));
	imum++;
      }
    }
    // G4cout<< ikp << " " << ikm << G4endl;
    // for( G4int i=0; i<ikp; i++) G4cout<<  gPfKP[i] << G4endl;
    // for( G4int i=0; i<ikm; i++) G4cout<<  gPfKM[i] << G4endl;
    // G4cout<< ipip << " " << ipim << G4endl;
    // for( G4int i=0; i<ipip; i++) G4cout<<  gPfPiP[i] << G4endl;
    // for( G4int i=0; i<ipim; i++) G4cout<<  gPfPiM[i] << G4endl;
    // G4cout<< ip << " " << ipb <<G4endl;
    // for( G4int i=0; i<ip; i++) G4cout<<  gPfP[i] << G4endl;
    // for( G4int i=0; i<ipb; i++) G4cout<<  gPfPb[i] << G4endl;
    
    G4ThreeVector gPos(LPos);
    gPos.rotateX(90.0*degree); gPos.rotateZ( BP->TargetRotAngH );
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
    LBeamDir.rotateX(90.0*degree); LBeamDir.rotateZ( BP->TargetRotAngH );
    
    //Data cut check
    G4bool flagD=true;
    G4bool flagDS=true;

    // G4bool flagD=false;
    // G4bool flagDS=false;

    // G4double KPMass =G4KaonPlus::Definition()->GetPDGMass();
    // G4double PiMMass=G4PionMinus::Definition()->GetPDGMass(); 
    // G4LorentzVector LvDM1, LvDSM1;
    // G4double QV=-1.;

    // //D0->K+ pi- mode: (K+, pi-, pis-)
    // if( confMan->DataCut()==1 ){      
    //   for( G4int ik=0; ik<ikp; ik++){
    // 	for( G4int ipi1=0; ipi1<ipim; ipi1++){
    // 	  for( G4int ipi2=0; ipi2<ipim; ipi2++){
    // 	    G4LorentzVector LvK( gPfKP[ik], sqrt(KPMass*KPMass+gPfKP[ik].mag2()) );
    // 	    G4LorentzVector LvPi1( gPfPiM[ipi1], sqrt(PiMMass*PiMMass+gPfPiM[ipi1].mag2()) );   
    // 	    G4LorentzVector LvPi2( gPfPiM[ipi2], sqrt(PiMMass*PiMMass+gPfPiM[ipi2].mag2()) );
	    
    // 	    LvDM1.setVect(gPfKP[ik]+gPfPiM[ipi1]);
    // 	    LvDM1.setE(LvK.e()+LvPi1.e());
	    
    // 	    if( MinDMMass<LvDM1.mag() && LvDM1.mag()<MaxDMMass ){
    // 	      flagD=true;
    // 	      LvDSM1.setVect(gPfKP[ik]+gPfPiM[ipi1]+gPfPiM[ipi2]);
    // 	      LvDSM1.setE(LvK.e()+LvPi1.e()+LvPi2.e());
    // 	      QV=LvDSM1.mag()-LvDM1.mag()-PiMMass;
    // 	      if( MinQV<QV && QV<MaxQV ){
    // 		flagDS=true;
    // 		// G4cout<< "*****D Mass " << LvDM1.mag()/GeV << G4endl;
    // 		// G4cout<< "*****D* Mass " << LvDSM1.mag()/GeV << G4endl;
    // 	      }
    // 	    }
    // 	  }
    // 	}
    //   }
    // }
    // //D0->K+ pi- mode: (pi+, pi-, pis-) for background
    // if( confMan->DataCut()==2 ){   
    //   for( G4int ik=0; ik<ipip; ik++){
    // 	for( G4int ipi1=0; ipi1<ipim; ipi1++){
    // 	  for( G4int ipi2=0; ipi2<ipim; ipi2++){
    // 	    G4LorentzVector LvK( gPfPiP[ik], sqrt(KPMass*KPMass+gPfPiP[ik].mag2()) );
    // 	    G4LorentzVector LvPi1( gPfPiM[ipi1], sqrt(PiMMass*PiMMass+gPfPiM[ipi1].mag2()) );   
    // 	    G4LorentzVector LvPi2( gPfPiM[ipi2], sqrt(PiMMass*PiMMass+gPfPiM[ipi2].mag2()) );
	    
    // 	    LvDM1.setVect(gPfKP[ik]+gPfPiM[ipi1]);
    // 	    LvDM1.setE(LvK.e()+LvPi1.e());

    // 	    if( MinDMMass<LvDM1.mag() && LvDM1.mag()<MaxDMMass ){
    // 	      flagD=true;
    // 	      LvDSM1.setVect(gPfKP[ik]+gPfPiM[ipi1]+gPfPiM[ipi2]);
    // 	      LvDSM1.setE(LvK.e()+LvPi1.e()+LvPi2.e());
    // 	      QV=LvDSM1.mag()-LvDM1.mag()-PiMMass;
    // 	      if( MinQV<QV && QV<MaxQV ){
    // 		flagDS=true;
    // 		// G4cout<< "*****D Mass " << LvDM1.mag()/GeV << G4endl;
    // 		// G4cout<< "*****D* Mass " << LvDSM1.mag()/GeV << G4endl;
    // 	      }
    // 	    }
    // 	  }
    // 	}
    //   }
    // }
    // //D0->K+ pi- mode: (p, pi-, pis-) for background    
    // if( confMan->DataCut()==3 ){        
    //   for( G4int ik=0; ik<ip; ik++){
    // 	for( G4int ipi1=0; ipi1<ipim; ipi1++){
    // 	  for( G4int ipi2=0; ipi2<ipim; ipi2++){
    // 	    G4LorentzVector LvK( gPfP[ik], sqrt(KPMass*KPMass+gPfP[ik].mag2()) );
    // 	    G4LorentzVector LvPi1( gPfPiM[ipi1], sqrt(PiMMass*PiMMass+gPfPiM[ipi1].mag2()) );   
    // 	    G4LorentzVector LvPi2( gPfPiM[ipi2], sqrt(PiMMass*PiMMass+gPfPiM[ipi2].mag2()) );
	    
    // 	    LvDM1.setVect(gPfKP[ik]+gPfPiM[ipi1]);
    // 	    LvDM1.setE(LvK.e()+LvPi1.e());
	    
    // 	    if( MinDMMass<LvDM1.mag() && LvDM1.mag()<MaxDMMass ){
    // 	      flagD=true;
    // 	      LvDSM1.setVect(gPfKP[ik]+gPfPiM[ipi1]+gPfPiM[ipi2]);
    // 	      LvDSM1.setE(LvK.e()+LvPi1.e()+LvPi2.e());
    // 	      QV=LvDSM1.mag()-LvDM1.mag()-PiMMass;
    // 	      if( MinQV<QV && QV<MaxQV ){
    // 		flagDS=true;
    // 		// G4cout<< "*****D Mass " << LvDM1.mag()/GeV << G4endl;
    // 		// G4cout<< "*****D* Mass " << LvDSM1.mag()/GeV << G4endl;
    // 	      }
    // 	    }
    // 	  }
    // 	}
    //   }
    // }
    // //D0->K+ pi- pi+ pi- mode
    // if( confMan->DataCut()==4 ){        
    //   for( G4int ik=0; ik<ikp; ik++){
    // 	for( G4int ipip1=0; ipip1<ipip; ipip1++){
    // 	  for( G4int ipim1=0; ipim1<ipim; ipim1++){
    // 	    for( G4int ipim2=0; ipim2<ipim; ipim2++){
    // 	      for( G4int ipim3=0; ipim3<ipim; ipim3++){
		
    // 		G4LorentzVector LvK( gPfKP[ik], sqrt(KPMass*KPMass+gPfKP[ik].mag2()) );
    // 		G4LorentzVector LvPiP1( gPfPiP[ipip1], sqrt(PiMMass*PiMMass+gPfPiP[ipip1].mag2()) );   
    // 		G4LorentzVector LvPiM1( gPfPiM[ipim1], sqrt(PiMMass*PiMMass+gPfPiM[ipim1].mag2()) );   
    // 		G4LorentzVector LvPiM2( gPfPiM[ipim2], sqrt(PiMMass*PiMMass+gPfPiM[ipim2].mag2()) );
    // 		G4LorentzVector LvPiM3( gPfPiM[ipim3], sqrt(PiMMass*PiMMass+gPfPiM[ipim3].mag2()) );
		
    // 		LvDM1.setVect(gPfKP[ik]+gPfPiP[ipip1]+gPfPiM[ipim1]+gPfPiM[ipim2]);
    // 		LvDM1.setE(LvK.e()+LvPiP1.e()+LvPiM1.e()+LvPiM2.e());
		
    // 		if( MinDMMass<LvDM1.mag() && LvDM1.mag()<MaxDMMass ){
    // 		  flagD=true;
    // 		  LvDSM1.setVect(gPfKP[ik]+gPfPiP[ipip1]+gPfPiM[ipim1]+gPfPiM[ipim2]+gPfPiM[ipim3]);
    // 		  LvDSM1.setE(LvK.e()+LvPiP1.e()+LvPiM1.e()+LvPiM2.e()+LvPiM3.e());
    // 		  QV=LvDSM1.mag()-LvDM1.mag()-PiMMass;
    // 		  if( MinQV<QV && QV<MaxQV ){
    // 		    flagDS=true;
    // 		    // G4cout<< "*****D Mass " << LvDM1.mag()/GeV << G4endl;
    // 		    // G4cout<< "*****D* Mass " << LvDSM1.mag()/GeV << G4endl;
    // 		  }
    // 		}
    // 	      }
    // 	    }
    // 	  }
    // 	}
    //   }
    // }
    // //D0->K+ pi- mode with 13 GeV/c BNL condtions
    // if( confMan->DataCut()==5 ){   
    //   G4bool momflag=false;
    //   G4bool thetakpiflag=false;
    //   G4bool ukpiflag=false;
    //   for( G4int ik=0; ik<ikp; ik++){
    // 	for( G4int ipi1=0; ipi1<ipim; ipi1++){
    // 	  for( G4int ipi2=ipi1; ipi2<ipim; ipi2++){
    // 	    G4LorentzVector LvK( gPfKP[ik], sqrt(KPMass*KPMass+gPfKP[ik].mag2()) );
    // 	    G4LorentzVector LvPi1( gPfPiM[ipi1], sqrt(PiMMass*PiMMass+gPfPiM[ipi1].mag2()) );   
    // 	    G4LorentzVector LvPi2( gPfPiM[ipi2], sqrt(PiMMass*PiMMass+gPfPiM[ipi2].mag2()) );
	    
    // 	    LvDM1.setVect( gPfKP[ik]+gPfPiM[ipi1] );
    // 	    LvDM1.setE( LvK.e()+LvPi1.e() );
	    
    // 	    LvDSM1.setVect( gPfKP[ik]+gPfPiM[ipi1]+gPfPiM[ipi2] );
    // 	    LvDSM1.setE( LvK.e()+LvPi1.e()+LvPi2.e() );
    
    // 	    G4double Q = LvDSM1.mag()-LvDM1.mag()-PiMMass;
    // 	    G4double DSMass = LvDSM1.mag()/GeV;	
    
    // 	    G4double pK  = gPfKP[ik].mag()/GeV;
    // 	    G4double pPi = gPfPiM[ipi1].mag()/GeV;
    // 	    G4double pKPi = pK+pPi;
	    
    // 	    if( (1.5 < pK && pK < 8.0) && 
    // 		(1.5 < pPi) && (6.0 < pKPi)){
    // 	      momflag=true;
    // 	    }
	    
    // 	    G4double costkpi1=gPfPiM[ipi1]*gPfKP[ik]/(gPfPiM[ipi1].mag()*gPfKP[ik].mag());
    // 	    G4double thetakpi1 = acos(costkpi1)*Rad2Deg;
	    
    // 	    if( thetakpi1 > 16.0 ){
    // 	      thetakpiflag=true;
    // 	    }
	    
    // 	    G4double uK = gPfKP[ik].x()/gPfKP[ik].mag();
    // 	    G4double uPi = gPfPiM[ipi1].x()/gPfPiM[ipi1].mag();
    // 	    G4double Up =  0.14;
    // 	    G4double Um = -0.14;
	    
    // 	    if(  !(Um<uK && uK<Up) && !(Um<uPi && uPi<Up) ){
    // 	      ukpiflag=true; 
    // 	    }
	    
    // 	    if( momflag && thetakpiflag && ukpiflag && 
    // 		3.2< Q && Q < 8.2){ 
    // 	      //2.000 < DSMass && DSMass <2.030 ){
    // 	      flagD=true; flagDS=true;
    // 	      // G4cout<< "*****D Mass " << LvDM1.mag()/GeV << G4endl;
    // 	      // G4cout<< "*****D* Mass " << LvDSM1.mag()/GeV << G4endl;
    // 	    }
    // 	  }
    // 	}
    //   }
    // }
    // else{
    //   flagD=true; flagDS=true;
    // }
    
    if( flagD && flagDS ){
      //G4cout<< "***" << G4endl;
      
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

      //Beam
      pGen->gunPiM_->SetParticlePosition( gPos+BP->gTargetPos );
      pGen->gunPiM_->SetParticleMomentum( -1.*pb*LBeamDir );
      pGen->gunPiM_->GeneratePrimaryVertex( anEvent );
      
      G4ThreeVector PisMom(0.,0.,0.);
      
      if( pGen->anaMan_ ){
	LBeamDir.rotateZ( -1.*(BP->TargetRotAngH) ); 
	LBeamDir.rotateX(-90.0*degree);
	PisMom.rotateZ( -1.*(BP->TargetRotAngH) ); 
	PisMom.rotateX(-90.0*degree);
	// pGen->SetPrimaryInformation( LPos, Mf1, Mf2, 
	// 				 PisMom, PisMom, PisMom, PisMom, 
	// 				 pb*LBeamDir, PisMom );
	pGen->SetPrimaryInformation( LPos, Mf1, Mf2, 
				     gPfKP[0], gPfPiM[0], gPfKP[0], gPfPiM[0], 
				     pb*LBeamDir, PisMom );
      }
    }
    //BackgroundFromPiMinusProtonJAM1
  }
}
