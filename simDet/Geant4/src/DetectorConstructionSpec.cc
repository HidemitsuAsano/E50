/*
  DetectorConstructionSpec.hh

  2016/4  K.Shirotori
  2016/5  update H. Asano
*/

#include "DetectorConstructionSpec.hh"

#include "MaterialList.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"

#include "G4TransportationManager.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"

#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Paraboloid.hh"
#include "G4Torus.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4BooleanSolid.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"

#include "G4Material.hh"

#include "s_DetectorSize.hh"
#include "s_SFT.hh"
#include "s_T0Wall.hh"

#include "r_DetectorSize.hh"
//#include "r_SFT_square.hh" //square ,1mm single cladding fiber. probably, not used anymore 
#include "r_SFT_Round.hh" //round  ,double cladding fiber, r=0.5 mm
                          //layer z position is defined in this code.

//#include "r_SFT_C.hh" //square ,0.5 mm single caldding fiber. probably, not used anymore
//#include "r_SFT_D.hh" //round , double cladding fiber r=0.25 mm (or 0.15 mm)

#include "SFTSD.hh"
#include "T0SD.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

DetectorConstructionSpec::DetectorConstructionSpec()
{
}

DetectorConstructionSpec::~DetectorConstructionSpec()
{
}

G4VPhysicalVolume *
DetectorConstructionSpec::ConstructPayload( void )
{
  ///////////////////// World
  G4Box *solidWorld =
    new G4Box( "WorldBox", 0.5*m, 0.5*m, 0.5*m );

  G4LogicalVolume *logWorld =  //Containes all info of volume except position
    new G4LogicalVolume( solidWorld, mList_->Air, "World", 0, 0, 0, false );
  logWorld->SetVisAttributes( G4VisAttributes::Invisible );

  G4VPhysicalVolume *physWorld =
    new G4PVPlacement( 0, G4ThreeVector(), "World", logWorld, 0, false, 0 );
  
  MakeTrackers( physWorld );
  MakeCounters( physWorld );
  
  //output geometry file which is used for the current GEANT4 simulation
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  geomMan.CreateParamFile("param.file.new");
  geomMan.PrintParameter();

  return physWorld;
}


void DetectorConstructionSpec::
MakeTrackers( G4VPhysicalVolume *pMother )
{
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  ConfMan *confMan = ConfMan::GetConfManager();
  G4RotationMatrix RotationMatrix;  
  
  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////

  G4int DetType=0;
  if( confMan->DetectorType()==0 ) DetType=0;//simple detector
  //else if( confMan->DetectorType()==1 ) DetType=1;//square fiber 1mm
  else if( confMan->DetectorType()==1 ) DetType=1;//round fiber 1mm
  else if( confMan->DetectorType()==2 ) DetType=2;//round fiber 1mm,  M-line configulation    
  else { G4cout << "invalid Detector type ! (DetType = 0, 1, 2) " << G4endl; return;}
  G4cout << "   "  << G4endl;
  G4cout << __FILE__ << "  " << __LINE__<< "  DetType " << confMan->DetectorType() << G4endl;
  G4cout << "   "  << G4endl;
  
  //Type 0 simple detector
  if( DetType==0 ){
    ///////////////////// SFT(x,u,v,x,u,v, u,v,x,u,v,x)
    G4int detid[rSFT_nLayer];
    detid[0]=geomMan.GetDetectorId("SFT-x-0");
    detid[1]=geomMan.GetDetectorId("SFT-u-0");
    detid[2]=geomMan.GetDetectorId("SFT-v-0");
    detid[3]=geomMan.GetDetectorId("SFT-x-1");
    detid[4]=geomMan.GetDetectorId("SFT-u-1");
    detid[5]=geomMan.GetDetectorId("SFT-v-1");
    detid[6]=geomMan.GetDetectorId("SFT-u-2");
    detid[7]=geomMan.GetDetectorId("SFT-v-2");
    detid[8]=geomMan.GetDetectorId("SFT-x-2");
    detid[9]=geomMan.GetDetectorId("SFT-u-3");
    detid[10]=geomMan.GetDetectorId("SFT-v-3");
    detid[11]=geomMan.GetDetectorId("SFT-x-3");
    
    G4ThreeVector gPos[rSFT_nLayer];
    for(int ilr=0;ilr<rSFT_nLayer;ilr++){
      gPos[ilr] = geomMan.GetGlobalPosition(detid[ilr]);
    }
    G4ThreeVector OfsLSFT( 0.0, 0.0, 0.0 );
    
    s_SFT *s_SF1 =
      new s_SFT( "SFT", pMother, RotationMatrix, 
		 gPos,
		 OfsLSFT, 
		 detid, 
		 mList_->Scin, 
		 mList_->Al, 
		 mList_->Air, 
		 mList_->Air );
    
    G4VisAttributes *FiberLayerAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.0, 1.0, 1.0 ) );//Cyan
//    G4VisAttributes *FiberFrameAttrib =
//      new G4VisAttributes( G4Colour::G4Colour( 0.3, 0.3, 0.3 ) );//Blac
    
    // s_SF1->SetVisAttributes( FiberLayerAttrib, 
    // 			     FiberFrameAttrib, 
    // 			     &G4VisAttributes::Invisible,
    // 			     &G4VisAttributes::Invisible );

    s_SF1->SetVisAttributes( FiberLayerAttrib, 
			     &G4VisAttributes::Invisible,
			     &G4VisAttributes::Invisible,
			     &G4VisAttributes::Invisible );
    
    //////// Sensitive Detectors ////////
    G4SDManager *SDMan = G4SDManager::GetSDMpointer();
    SFTSD *sftSD = new SFTSD( "/spec/sft" );
    SDMan->AddNewDetector( sftSD );  
    s_SF1->SetSensitiveDetector( sftSD );
  }
  
  /*
  //realistic SFT Square Fiber 1mm Single Cladding
  else if( DetType==1 ){
    ///////////////////// SFT(x,u,v,x,u,v, u,v,x,u,v,x)
    G4int detid[rSFT_nLayer];
    detid[0]=geomMan.GetDetectorId("SFT-x-0");
    detid[1]=geomMan.GetDetectorId("SFT-u-0");
    detid[2]=geomMan.GetDetectorId("SFT-v-0");
    detid[3]=geomMan.GetDetectorId("SFT-x-1");
    detid[4]=geomMan.GetDetectorId("SFT-u-1");
    detid[5]=geomMan.GetDetectorId("SFT-v-1");
    detid[6]=geomMan.GetDetectorId("SFT-u-2");
    detid[7]=geomMan.GetDetectorId("SFT-v-2");
    detid[8]=geomMan.GetDetectorId("SFT-x-2");
    detid[9]=geomMan.GetDetectorId("SFT-u-3");
    detid[10]=geomMan.GetDetectorId("SFT-v-3");
    detid[11]=geomMan.GetDetectorId("SFT-x-3");
    
    G4ThreeVector gPos[rSFT_nLayer];
    for(int ilr=0;ilr<rSFT_nLayer;ilr++){
      gPos[ilr] = geomMan.GetGlobalPosition(detid[ilr]);
    }
    G4ThreeVector OfsLSFT( 0.0, 0.0, 0.0 );

    r_SFT_Square *r_SFT_Square =
      new r_SFT_Square( "SFT", pMother, RotationMatrix, 
		   gPos, 
		   OfsLSFT, 
		   detid, 
		   mList_->Scin, 
		   mList_->PMMA, 
		   mList_->Al, 
		   mList_->Air, 
		   mList_->Air );
    
    G4VisAttributes *FiberLayerAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.0, 1.0, 1.0 ) );//Cyan
    G4VisAttributes *FiberCradAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 1.0, 0.0, 0.0 ) );//Red
    G4VisAttributes *FiberFrameAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.3, 0.3, 0.3 ) );//Blac
    
    // r_SF1->SetVisAttributes( FiberLayerAttrib, 
    // 			     FiberFrameAttrib, 
    // 			     &G4VisAttributes::Invisible,
    // 			     &G4VisAttributes::Invisible );

    r_SF1->SetVisAttributes( FiberLayerAttrib, 
           FiberCradAttrib,
			     &G4VisAttributes::Invisible,
			     &G4VisAttributes::Invisible,
			     &G4VisAttributes::Invisible );

    //////// Sensitive Detectors ////////
    G4SDManager *SDMan = G4SDManager::GetSDMpointer();
    SFTSD *sftSD = new SFTSD( "/spec/sft" );
    SDMan->AddNewDetector( sftSD );  
    r_SF1->SetSensitiveDetector( sftSD );
  }
  */ 

  //realistic SFT Round Fiber, round and double cladding structure
  else if( DetType==1 || DetType==2){
    ///////////////////// SFT(x,u,v,x,u,v, u,v,x,u,v,x)
    G4int detid[rSFT_nLayer];
    detid[0]=geomMan.GetDetectorId("SFT-x-0");
    detid[1]=geomMan.GetDetectorId("SFT-u-0");
    detid[2]=geomMan.GetDetectorId("SFT-v-0");
    detid[3]=geomMan.GetDetectorId("SFT-x-1");
    detid[4]=geomMan.GetDetectorId("SFT-u-1");
    detid[5]=geomMan.GetDetectorId("SFT-v-1");
    detid[6]=geomMan.GetDetectorId("SFT-u-2");
    detid[7]=geomMan.GetDetectorId("SFT-v-2");
    detid[8]=geomMan.GetDetectorId("SFT-x-2");
    detid[9]=geomMan.GetDetectorId("SFT-u-3");
    detid[10]=geomMan.GetDetectorId("SFT-v-3");
    detid[11]=geomMan.GetDetectorId("SFT-x-3");
    

    G4ThreeVector gPos[rSFT_nLayer];
    for(int ilr=0;ilr<rSFT_nLayer;ilr++){
      gPos[ilr] = geomMan.GetGlobalPosition(detid[ilr]);
    }
    G4ThreeVector OfsLSFT( 0.0, 0.0, 0.0 );

    //default 

    G4double LzLayer[rSFT_nLayer]={0.};
    
    if(DetType==1){
      LzLayer[0]= -165.0*mm;//x
      LzLayer[1]= -135.0*mm;//u
      LzLayer[2]= -105.0*mm;//v
      LzLayer[3]= -75.0*mm;//x
      LzLayer[4]= -45.0*mm;//u
      LzLayer[5]= -15.0*mm;//v
      LzLayer[6]=  15.0*mm;//u
      LzLayer[7]=  45.0*mm;//v
      LzLayer[8]=  75.0*mm;//x
      LzLayer[9]= 105.0*mm;//u
      LzLayer[10] = 135.0*mm;//v
      LzLayer[11] =  165.0*mm;//x
    }else if(DetType==2){//x,u,v layers are close to each other
      LzLayer[0]= -166.0*mm;//x
      LzLayer[1]= -164.0*mm;//u
      LzLayer[2]= -162.0*mm;//v
      LzLayer[3]= -47.0*mm;//x
      LzLayer[4]= -45.0*mm;//u
      LzLayer[5]= -43.0*mm;//v
      LzLayer[6]=  43.0*mm;//u
      LzLayer[7]=  45.0*mm;//v
      LzLayer[8]=  47.0*mm;//x
      LzLayer[9]= 162.0*mm;//u
      LzLayer[10] = 164.0*mm;//v
      LzLayer[11] =  166.0*mm;//x
    }

    r_SFT_Round *r_SFT =
      new r_SFT_Round( "SFT", pMother, RotationMatrix, 
		   gPos, 
		   OfsLSFT, 
		   detid,
       LzLayer,
		   mList_->Scin, 
		   mList_->PMMA, 
		   mList_->FP, 
		   mList_->Al, 
		   mList_->Air, 
		   mList_->Air );
    
    G4VisAttributes *FiberLayerAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.0, 1.0, 1.0 ) );//Cyan
    G4VisAttributes *FiberInnerCradAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 1.0, 0.0, 0.0 ) );//Red
    G4VisAttributes *FiberOuterCradAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 1.0, 1.0, 0.0 ) );//yellow
//    G4VisAttributes *FiberFrameAttrib =
//      new G4VisAttributes( G4Colour::G4Colour( 0.3, 0.3, 0.3 ) );//Blac
    
    // r_SF1->SetVisAttributes( FiberLayerAttrib, 
    // 			     FiberFrameAttrib, 
    // 			     &G4VisAttributes::Invisible,
    // 			     &G4VisAttributes::Invisible );

    r_SFT->SetVisAttributes( FiberLayerAttrib, 
           FiberInnerCradAttrib,
           FiberOuterCradAttrib,
			     &G4VisAttributes::Invisible,
			     &G4VisAttributes::Invisible,
			     &G4VisAttributes::Invisible );

    //////// Sensitive Detectors ////////
    G4SDManager *SDMan = G4SDManager::GetSDMpointer();
    SFTSD *sftSD = new SFTSD( "/spec/sft" );
    SDMan->AddNewDetector( sftSD );  
    r_SFT->SetSensitiveDetector( sftSD );
  }
   
  /*
  //Type C real SFT Square Fiber 0.5 mm Single Cladding
  else if( DetType==3 ){
    ///////////////////// SFT(x,u,v,x,u,v, u,v,x,u,v,x)
    G4int detid[12];
    detid[0]=geomMan.GetDetectorId("SFT-x-1");
    detid[1]=geomMan.GetDetectorId("SFT-u-1");
    detid[2]=geomMan.GetDetectorId("SFT-v-1");
    detid[3]=geomMan.GetDetectorId("SFT-x-2");
    detid[4]=geomMan.GetDetectorId("SFT-u-2");
    detid[5]=geomMan.GetDetectorId("SFT-v-2");
    detid[6]=geomMan.GetDetectorId("SFT-u-3");
    detid[7]=geomMan.GetDetectorId("SFT-v-3");
    detid[8]=geomMan.GetDetectorId("SFT-x-3");
    detid[9]=geomMan.GetDetectorId("SFT-u-4");
    detid[10]=geomMan.GetDetectorId("SFT-v-4");
    detid[11]=geomMan.GetDetectorId("SFT-x-4");
    
    G4ThreeVector gPos[12];
    for(int ilr=0;ilr<12;ilr++){
      gPos[ilr] = geomMan.GetGlobalPosition(detid[ilr]);
    }
    G4ThreeVector OfsLSFT( 0.0, 0.0, 0.0 );

    r_SFT_C *r_SF1 =
      new r_SFT_C( "SFT", pMother, RotationMatrix, 
		   gPos, 
		   OfsLSFT, 
		   detid,  
		   mList_->Scin, 
		   mList_->PMMA, 
		   mList_->Al, 
		   mList_->Air, 
		   mList_->Air );
    
    G4VisAttributes *FiberLayerAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.0, 1.0, 1.0 ) );//Cyan
    G4VisAttributes *FiberCradAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 1.0, 0.0, 0.0 ) );//Red
    G4VisAttributes *FiberFrameAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.3, 0.3, 0.3 ) );//Blac
    
    // r_SF1->SetVisAttributes( FiberLayerAttrib, 
    // 			     FiberFrameAttrib, 
    // 			     &G4VisAttributes::Invisible,
    // 			     &G4VisAttributes::Invisible );

    r_SF1->SetVisAttributes( FiberLayerAttrib, 
           FiberCradAttrib,
			     &G4VisAttributes::Invisible,
			     &G4VisAttributes::Invisible,
			     &G4VisAttributes::Invisible );

    //////// Sensitive Detectors ////////
    G4SDManager *SDMan = G4SDManager::GetSDMpointer();
    SFTSD *sftSD = new SFTSD( "/spec/sft" );
    SDMan->AddNewDetector( sftSD );  
    r_SF1->SetSensitiveDetector( sftSD );
  }
  */

}

void DetectorConstructionSpec::
MakeCounters( G4VPhysicalVolume *pMother )
{
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  //ConfMan *confMan = ConfMan::GetConfManager();
  G4RotationMatrix RotationMatrix;  

  ///////////////////////////////////////// 
  ////////////////////////////////T0: Time Zero counter
  ///////////////////////////////////////// 
    
  G4int idt0=geomMan.GetDetectorId("T0");
  G4ThreeVector gPosT0=geomMan.GetGlobalPosition(idt0);
  G4ThreeVector OfsLT0( 0.0, 0.0, 0.0 );
  
  s_T0Wall *s_T0 =
    new s_T0Wall( "T0", pMother, RotationMatrix, gPosT0, OfsLT0, 
		  mList_->Scin, 
		  mList_->Air );
  
  G4VisAttributes *T0ScinAttrib =
    new G4VisAttributes( G4Colour::G4Colour( 1.0, 0.0, 1.0 ) );//Magenta
  
  s_T0->SetVisAttributes( T0ScinAttrib, 
			  &G4VisAttributes::Invisible );
  
  //////// Sensitive Detectors ////////
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  T0SD *t0SD = new T0SD( "/spec/t0" );
  SDMan->AddNewDetector( t0SD );
  s_T0->SetSensitiveDetector( t0SD );
  
}

