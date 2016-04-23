/*
  ReadJAMData.cc

  2013/4  K.Shirotori
*/

#include "ReadJAMData.hh"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "JAMData.hh"
#include "TemplateLib.hh"
#include "ParticleID.hh"

#include "ConfMan.hh"

#include "G4ParticleTypes.hh"
#include "G4LorentzVector.hh"

const G4double Deg2Rad = acos(-1.)/180.;
const G4double Rad2Deg = 180./acos(-1.);

//D0 Mass cut -30 -> +30 MeV
const G4double MinDMMass = 1835.;
const G4double MaxDMMass = 1895.;

//K*0 Mass cut
const G4double MinKSMass = 840.;
const G4double MaxKSMass = 950.;

//Q-value cut -5 -> +10 MeV
const G4double MinQV =  0.;
const G4double MaxQV = 15.;

ReadJAMData::ReadJAMData():
  Data(0), evnum_(0)
{}

ReadJAMData::~ReadJAMData()
{
  clearAll();
}

bool ReadJAMData::AddData( JAMDataContainer& cont,
			   G4int EventNum, G4int Pid, 
			   G4double MomX, G4double MomY, G4double MomZ )
{
  static const std::string funcname = "[ReadJAMData::AddData]";
  //std::cout << "**************************=" << EventNum << std::endl;

  JAMData *p=0;
  G4int nh=cont.size();
  for( G4int i=EventNum; i<nh; ++i ){
    JAMData *q=cont[i];
    if( q->EventNo()==EventNum ){
      p=q; break;
    }
  }
  if(!p){
    p = new JAMData( EventNum );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetPid( Pid );
    p->SetMomX( MomX );
    p->SetMomY( MomY );
    p->SetMomZ( MomZ );

    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

void ReadJAMData::clearAll()
{
  std::for_each(Data.begin(), Data.end(), DeleteObject());
  Data.clear();
    
  return;
}

bool ReadJAMData::Processing( std::istream &In, G4int EventNum )
{
  //clearAll();

  ConfMan *confMan = ConfMan::GetConfManager();

  //std::cout << "**************************" << std::endl;
  G4int event, numpart;
  G4int numbary, nummeson;
  
  //////////////////////////////////////////////////////////////
  In >> event >> numpart >> numbary >> nummeson;
  if( event==-9 ) return false;
  
  G4int lnum;
  In >> lnum;
  volatile G4int lnumold=lnum;

  std::vector <G4ThreeVector> pKP, pKM;
  std::vector <G4ThreeVector> pPiP, pPiM;
  std::vector <G4ThreeVector> pP, pPb;
  std::vector <G4ThreeVector> pPiZ;
  std::vector <G4ThreeVector> pN, pNb;
  std::vector <G4ThreeVector> pKZ;
  std::vector <G4ThreeVector> pKZS, pKZL;
  std::vector <G4ThreeVector> pG;
  std::vector <G4ThreeVector> pEP, pEM;
  std::vector <G4ThreeVector> pMuP, pMuM;

  while( lnum!=-1 ){
    //std::cout << lnumold << "**" << lnum << std::endl;

    //Kaon+
    if( lnumold==IdKaonP ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pKP.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;
	
    //Kaon-
    if( lnumold==IdKaonM ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pKM.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Pion+
    if( lnumold==IdPionP ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pPiP.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Pion-
    if( lnumold==IdPionM ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pPiM.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Proton
    if( lnumold==IdProton ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pP.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Proton bar
    if( lnumold==IdProtonb ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pPb.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Pion0
    if( lnumold==IdPionZ ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pPiZ.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Neutron
    if( lnumold==IdNeutron ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pN.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Neutron bar
    if( lnumold==IdNeutronb ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pNb.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Kaon0
    if( lnumold==IdKaonZ ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pKZ.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Kaon0S
    if( lnumold==IdKaonZS ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pKZS.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Kaon0L
    if( lnumold==IdKaonZL ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pKZL.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Gamma
    if( lnumold==IdGamma ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pG.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Positron
    if( lnumold==IdPosi ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pEP.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Electron
    if( lnumold==IdElec ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pEM.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Muon+
    if( lnumold==IdMuonP ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pMuP.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;

    //Muon-
    if( lnumold==IdMuonM ){
      G4double px, py, pz;
      In >> px >> py >> pz >> lnum;      

      G4ThreeVector Mom(px*GeV, py*GeV, pz*GeV);
      pMuM.push_back(Mom);
      // std::cout << lnumold 
      // 		<< " " << px << " " << py << " " << pz << std::endl;
    }
    lnumold = lnum;
  }

  //Data cut for save PC memory
  G4bool flagD=false;
  G4bool flagDS=false;
  
  G4double KPMass =G4KaonPlus::Definition()->GetPDGMass();
  G4double PiMMass=G4PionMinus::Definition()->GetPDGMass(); 
  G4LorentzVector LvDM1, LvDSM1;
  G4double QV=-1.;

  //D0->K+ pi- mode: (K+, pi-, pis-)
  if( confMan->DataCut()==1 ){    
    G4int ikp = pKP.size(), ipim= pPiM.size();
    for( G4int ik=0; ik<ikp; ik++){
      for( G4int ipi1=0; ipi1<ipim; ipi1++){
	for( G4int ipi2=0; ipi2<ipim; ipi2++){
	  G4LorentzVector LvK( pKP[ik], sqrt(KPMass*KPMass+pKP[ik].mag2()) );
	  G4LorentzVector LvPi1( pPiM[ipi1], sqrt(PiMMass*PiMMass+pPiM[ipi1].mag2()) );   
	  G4LorentzVector LvPi2( pPiM[ipi2], sqrt(PiMMass*PiMMass+pPiM[ipi2].mag2()) );
	  
	  LvDM1.setVect(pKP[ik]+pPiM[ipi1]);
	  LvDM1.setE(LvK.e()+LvPi1.e());
	  
	  if( MinDMMass<LvDM1.mag() && LvDM1.mag()<MaxDMMass ){ 
	    flagD=true; 
	    LvDSM1.setVect(pKP[ik]+pPiM[ipi1]+pPiM[ipi2]);
	    LvDSM1.setE(LvK.e()+LvPi1.e()+LvPi2.e());
	    QV=LvDSM1.mag()-LvDM1.mag()-PiMMass;
	    if( MinQV<QV && QV<MaxQV ){
	      flagDS=true;
	      // G4cout<< "*****D Mass " << LvDM1.mag()/GeV << G4endl;
	      // G4cout<< "*****D* Mass " << LvDSM1.mag()/GeV << G4endl;
	    }
	  }
	}
      }
    }
  }
  //D0->K+ pi- mode: (pi+, pi-, pis-)
  if( confMan->DataCut()==2 ){   
    G4int ipip = pPiP.size(), ipim= pPiM.size();
    for( G4int ik=0; ik<ipip; ik++){
      for( G4int ipi1=0; ipi1<ipim; ipi1++){
	for( G4int ipi2=0; ipi2<ipim; ipi2++){
	  G4LorentzVector LvK( pPiP[ik], sqrt(KPMass*KPMass+pPiP[ik].mag2()) );
	  G4LorentzVector LvPi1( pPiM[ipi1], sqrt(PiMMass*PiMMass+pPiM[ipi1].mag2()) );   
	  G4LorentzVector LvPi2( pPiM[ipi2], sqrt(PiMMass*PiMMass+pPiM[ipi2].mag2()) );
	  
	  LvDM1.setVect(pPiP[ik]+pPiM[ipi1]);
	  LvDM1.setE(LvK.e()+LvPi1.e());
	  
	  if( MinDMMass<LvDM1.mag() && LvDM1.mag()<MaxDMMass ){ 
	    flagD=true; 
	    LvDSM1.setVect(pPiP[ik]+pPiM[ipi1]+pPiM[ipi2]);
	    LvDSM1.setE(LvK.e()+LvPi1.e()+LvPi2.e());
	    QV=LvDSM1.mag()-LvDM1.mag()-PiMMass;
	    if( MinQV<QV && QV<MaxQV ){
	      flagDS=true;
	      // G4cout<< "*****D Mass " << LvDM1.mag()/GeV << G4endl;
	      // G4cout<< "*****D* Mass " << LvDSM1.mag()/GeV << G4endl;
	    }
	  }
	}
      }
    }
  }
  //D0->K+ pi- mode: (p, pi-, pis-)
  if( confMan->DataCut()==3 ){   
    G4int ip = pP.size(), ipim= pPiM.size();
    for( G4int ik=0; ik<ip; ik++){
      for( G4int ipi1=0; ipi1<ipim; ipi1++){
	for( G4int ipi2=0; ipi2<ipim; ipi2++){
	  G4LorentzVector LvK( pP[ik], sqrt(KPMass*KPMass+pP[ik].mag2()) );
	  G4LorentzVector LvPi1( pPiM[ipi1], sqrt(PiMMass*PiMMass+pPiM[ipi1].mag2()) );   
	  G4LorentzVector LvPi2( pPiM[ipi2], sqrt(PiMMass*PiMMass+pPiM[ipi2].mag2()) );
	  
	  LvDM1.setVect(pP[ik]+pPiM[ipi1]);
	  LvDM1.setE(LvK.e()+LvPi1.e());
	  
	  if( MinDMMass<LvDM1.mag() && LvDM1.mag()<MaxDMMass ){ 
	    flagD=true; 
	    LvDSM1.setVect(pP[ik]+pPiM[ipi1]+pPiM[ipi2]);
	    LvDSM1.setE(LvK.e()+LvPi1.e()+LvPi2.e());
	    QV=LvDSM1.mag()-LvDM1.mag()-PiMMass;
	    if( MinQV<QV && QV<MaxQV ){
	      flagDS=true;
	      // G4cout<< "*****D Mass " << LvDM1.mag()/GeV << G4endl;
	      // G4cout<< "*****D* Mass " << LvDSM1.mag()/GeV << G4endl;
	    }
	  }
	}
      }
    }
  }
  //D0->K+ pi- pi+ pi- mode
  if( confMan->DataCut()==4 ){   
    G4int ikp = pKP.size(), ipim= pPiM.size();
    G4int ipip = pPiP.size();
    for( G4int ik=0; ik<ikp; ik++){
      for( G4int ipip1=0; ipip1<ipip; ipip1++){
	for( G4int ipim1=0; ipim1<ipim; ipim1++){
	  for( G4int ipim2=0; ipim2<ipim; ipim2++){
	    for( G4int ipim3=0; ipim3<ipim; ipim3++){
	      
	      G4LorentzVector LvK( pKP[ik], sqrt(KPMass*KPMass+pKP[ik].mag2()) );
	      G4LorentzVector LvPiP1( pPiP[ipip1], sqrt(PiMMass*PiMMass+pPiP[ipip1].mag2()) );   
	      G4LorentzVector LvPiM1( pPiM[ipim1], sqrt(PiMMass*PiMMass+pPiM[ipim1].mag2()) );   
	      G4LorentzVector LvPiM2( pPiM[ipim2], sqrt(PiMMass*PiMMass+pPiM[ipim2].mag2()) );
	      G4LorentzVector LvPiM3( pPiM[ipim3], sqrt(PiMMass*PiMMass+pPiM[ipim3].mag2()) );
	      
	      LvDM1.setVect(pKP[ik]+pPiP[ipip1]+pPiM[ipim1]+pPiM[ipim2]);
	      LvDM1.setE(LvK.e()+LvPiP1.e()+LvPiM1.e()+LvPiM2.e());
	      
	      if( MinDMMass<LvDM1.mag() && LvDM1.mag()<MaxDMMass ){ 
		flagD=true; 
		LvDSM1.setVect(pKP[ik]+pPiP[ipip1]+pPiM[ipim1]+pPiM[ipim2]+pPiM[ipim3]);
		LvDSM1.setE(LvK.e()+LvPiP1.e()+LvPiM1.e()+LvPiM2.e()+LvPiM3.e());
		QV=LvDSM1.mag()-LvDM1.mag()-PiMMass;
		if( MinQV<QV && QV<MaxQV ){
		  flagDS=true;
		  // G4cout<< "*****D Mass " << LvDM1.mag()/GeV << G4endl;
		  // G4cout<< "*****D* Mass " << LvDSM1.mag()/GeV << G4endl;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  // For 13 GeV/c BNL data 
  if( confMan->DataCut()==5 ){   
    G4bool momflag=false;
    G4bool thetakpiflag=false;
    G4bool ukpiflag=false;
    G4int ikp = pKP.size(), ipim= pPiM.size();
    for( G4int ik=0; ik<ikp; ik++){
      for( G4int ipi1=0; ipi1<ipim; ipi1++){
	for( G4int ipi2=0; ipi2<ipim; ipi2++){
	  G4LorentzVector LvK( pKP[ik], sqrt(KPMass*KPMass+pKP[ik].mag2()) );
	  G4LorentzVector LvPi1( pPiM[ipi1], sqrt(PiMMass*PiMMass+pPiM[ipi1].mag2()) );   
	  G4LorentzVector LvPi2( pPiM[ipi2], sqrt(PiMMass*PiMMass+pPiM[ipi2].mag2()) );
	  
	  LvDM1.setVect(pKP[ik]+pPiM[ipi1]);
	  LvDM1.setE(LvK.e()+LvPi1.e());
	  
	  LvDSM1.setVect(pKP[ik]+pPiM[ipi1]+pPiM[ipi2]);
	  LvDSM1.setE(LvK.e()+LvPi1.e()+LvPi2.e());
	  
	  QV=LvDSM1.mag()-LvDM1.mag()-PiMMass;
	  
	  G4double pK  = pKP[ik].mag()/GeV;
	  G4double pPi = pPiM[ipi1].mag()/GeV;
	  G4double pKPi = pK+pPi;
	  
	  if( (1.5 < pK && pK < 8.0) && (1.5 < pPi) && (6.0 < pKPi)){
	    momflag=true;
	  }  
	  
	  G4double costkpi1=pPiM[ipi1]*pKP[ik]/(pPiM[ipi1].mag()*pKP[ik].mag());
	  G4double thetakpi1 = acos(costkpi1)*Rad2Deg;
	  
	  // if( thetakpi1 > 10.0 ){
	  thetakpiflag=true;
	  // }
	  
	  G4double uK = pKP[ik].x()/pKP[ik].mag();
	  G4double uPi = pPiM[ipi1].x()/pPiM[ipi1].mag();
	  G4double Up =  0.10;
	  G4double Um = -0.10;
	  
	  // if(  !(Um<uK && uK<Up) && !(Um<uPi && uPi<Up) ){
	  ukpiflag=true; 
	  // }
	  
	  if( momflag && thetakpiflag && ukpiflag && 2.9< QV && QV < 8.5 ){ 
	    flagD=true; 
	    flagDS=true;
	    // G4cout<< "*****D Mass " << LvDM1.mag()/GeV << G4endl;
	    // G4cout<< "*****D* Mass " << LvDSM1.mag()/GeV << G4endl;
	  }
	}
      }
    }
  }
  //K*0 -> K+ pi- mode: (K+, pi-)
  if( confMan->DataCut()==6 ){   
    G4int ikp = pKP.size(), ipim= pPiM.size();
    for( G4int ik=0; ik<ikp; ik++){
      for( G4int ipi1=0; ipi1<ipim; ipi1++){
	G4LorentzVector LvK( pKP[ik], sqrt(KPMass*KPMass+pKP[ik].mag2()) );
	G4LorentzVector LvPi1( pPiM[ipi1], sqrt(PiMMass*PiMMass+pPiM[ipi1].mag2()) );   
	
	LvDM1.setVect(pKP[ik]+pPiM[ipi1]);
	LvDM1.setE(LvK.e()+LvPi1.e());
	
	if( MinDMMass<LvDM1.mag() && LvDM1.mag()<MaxDMMass ){ 
	  flagD=true; 
	  flagDS=true;
	  // G4cout<< "*****D Mass " << LvDM1.mag()/GeV << G4endl;
	  // G4cout<< "*****D* Mass " << LvDSM1.mag()/GeV << G4endl;
	}  	  
      }
    }
  } 
  if( confMan->DataCut()==0 ){   
    flagD=true; flagDS=true;
  }

  if( !flagD && !flagDS ){
    //AddData( Data, EventNum, 0, 0., 0., 0. );
    return true;
  }

  if( flagD && flagDS ){
    for( G4int i=0; i<G4int(pKP.size()); i++)  AddData( Data, evnum_, IdKaonP,   pKP[i].x(),  pKP[i].y(),  pKP[i].z() );
    for( G4int i=0; i<G4int(pKM.size()); i++)  AddData( Data, evnum_, IdKaonM,   pKM[i].x(),  pKM[i].y(),  pKM[i].z() );
    for( G4int i=0; i<G4int(pPiP.size()); i++) AddData( Data, evnum_, IdPionP,   pPiP[i].x(), pPiP[i].y(), pPiP[i].z() );
    for( G4int i=0; i<G4int(pPiM.size()); i++) AddData( Data, evnum_, IdPionM,   pPiM[i].x(), pPiM[i].y(), pPiM[i].z() );
    for( G4int i=0; i<G4int(pP.size()); i++)   AddData( Data, evnum_, IdProton,  pP[i].x(),   pP[i].y(),   pP[i].z() );
    for( G4int i=0; i<G4int(pPb.size()); i++)  AddData( Data, evnum_, IdProtonb, pPb[i].x(),  pPb[i].y(),  pPb[i].z() );  
    for( G4int i=0; i<G4int(pPiZ.size()); i++) AddData( Data, evnum_, IdPionZ,   pPiZ[i].x(), pPiZ[i].y(), pPiZ[i].z() );
    for( G4int i=0; i<G4int(pN.size()); i++)   AddData( Data, evnum_, IdNeutron, pN[i].x(),   pN[i].y(),   pN[i].z() );
    for( G4int i=0; i<G4int(pNb.size()); i++)  AddData( Data, evnum_, IdNeutronb,pNb[i].x(),  pNb[i].y(),  pNb[i].z() );  
    for( G4int i=0; i<G4int(pKZ.size()); i++)  AddData( Data, evnum_, IdKaonZ,   pKZ[i].x(),  pKZ[i].y(),  pKZ[i].z() );
    for( G4int i=0; i<G4int(pKZS.size()); i++) AddData( Data, evnum_, IdKaonZS,  pKZS[i].x(), pKZS[i].y(), pKZS[i].z() );
    for( G4int i=0; i<G4int(pKZL.size()); i++) AddData( Data, evnum_, IdKaonZL,  pKZL[i].x(), pKZL[i].y(), pKZL[i].z() );
    for( G4int i=0; i<G4int(pG.size()); i++)   AddData( Data, evnum_, IdGamma,   pG[i].x(),   pG[i].y(),   pG[i].z() );
    for( G4int i=0; i<G4int(pEP.size()); i++)  AddData( Data, evnum_, IdPosi,    pEP[i].x(),  pEP[i].y(),  pEP[i].z() );
    for( G4int i=0; i<G4int(pEM.size()); i++)  AddData( Data, evnum_, IdElec,    pEM[i].x(),  pEM[i].y(),  pEM[i].z() );
    for( G4int i=0; i<G4int(pMuP.size()); i++) AddData( Data, evnum_, IdMuonP,   pMuP[i].x(), pMuP[i].y(), pMuP[i].z() );
    for( G4int i=0; i<G4int(pMuM.size()); i++) AddData( Data, evnum_, IdMuonM,   pMuM[i].x(), pMuM[i].y(), pMuM[i].z() );  
    evnum_++;
  }

  return true;
}

const JAMDataContainer& ReadJAMData::GetJAMData() const
{
  return Data;
}
