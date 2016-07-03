/*
  DetectorConstructionSpec.hh

  2016/4  K.Shirotori
  2016/4  update H. Asano
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
#include "r_SFT_A.hh" //square ,1mm signle cladding fiber
#include "r_SFT_B.hh" //round  ,double cladding fiber, r=0.5 mm
                      //layer z position is adjuested in this code.
//#include "r_SFT_C.hh" //square ,0.5 mm single caldding fiber
#include "r_SFT_D.hh" //round , double cladding fiber r=0.25 mm (or 0.15 mm)

#include "SFTSD.hh"
#include "T0SD.hh"
#include "TargetSD.hh"

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
  
//  MakeSpecMagnetandTarget( physWorld );
  MakeTrackers( physWorld );
  MakeCounters( physWorld );

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  geomMan.CreateParamFile("param.file.new");
  return physWorld;
}


void DetectorConstructionSpec::
MakeTrackers( G4VPhysicalVolume *pMother )
{
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  ConfMan *confMan = ConfMan::GetConfManager();
  G4RotationMatrix RM;  
  

  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////

  G4int DetType=0;
  if( confMan->DetectorType()==0 ) DetType=0;//simple detector
  else if( confMan->DetectorType()==1 ) DetType=1;//square fiber 1mm
  else if( confMan->DetectorType()==2 ) DetType=2;//round fiber 1mm
  else if( confMan->DetectorType()==3 ) DetType=3;//round fiber 1mm M-line configulation    //square fiber 0.5mm (obsolete)
  else { G4cout << "invalid Detector type ! (DetType = 0, 1, 2) " << G4endl; return;}
  G4cout << "   "  << G4endl;
  G4cout << "DetType " << confMan->DetectorType() << G4endl;
  G4cout << "   "  << G4endl;
  //Type 0
  if( DetType==0 ){
    ///////////////////// SFT(x,u,v,x,u,v, u,v,x,u,v,x)
    G4int id1=geomMan.GetDetectorId("SFT-x-1");
    G4int id2=geomMan.GetDetectorId("SFT-u-1");
    G4int id3=geomMan.GetDetectorId("SFT-v-1");
    G4int id4=geomMan.GetDetectorId("SFT-x-2");
    G4int id5=geomMan.GetDetectorId("SFT-u-2");
    G4int id6=geomMan.GetDetectorId("SFT-v-2");
    G4int id7=geomMan.GetDetectorId("SFT-u-3");
    G4int id8=geomMan.GetDetectorId("SFT-v-3");
    G4int id9=geomMan.GetDetectorId("SFT-x-3");
    G4int id10=geomMan.GetDetectorId("SFT-u-4");
    G4int id11=geomMan.GetDetectorId("SFT-v-4");
    G4int id12=geomMan.GetDetectorId("SFT-x-4");
    
    G4ThreeVector gPos1=geomMan.GetGlobalPosition(id1);
    G4ThreeVector gPos2=geomMan.GetGlobalPosition(id2);
    G4ThreeVector gPos3=geomMan.GetGlobalPosition(id3);
    G4ThreeVector gPos4=geomMan.GetGlobalPosition(id4);
    G4ThreeVector gPos5=geomMan.GetGlobalPosition(id5);
    G4ThreeVector gPos6=geomMan.GetGlobalPosition(id6);
    G4ThreeVector gPos7=geomMan.GetGlobalPosition(id7);
    G4ThreeVector gPos8=geomMan.GetGlobalPosition(id8);
    G4ThreeVector gPos9=geomMan.GetGlobalPosition(id9);
    G4ThreeVector gPos10=geomMan.GetGlobalPosition(id10);
    G4ThreeVector gPos11=geomMan.GetGlobalPosition(id11);
    G4ThreeVector gPos12=geomMan.GetGlobalPosition(id12);
    G4ThreeVector OfsLSFT( 0.0, 0.0, 0.0 );
    
    s_SFT *s_SF1 =
      new s_SFT( "SFT", pMother, RM, 
		 gPos1,  gPos2,  gPos3, 
		 gPos4,  gPos5,  gPos6,
		 gPos7,  gPos8,  gPos9, 
		 gPos10, gPos11, gPos12,
		 OfsLSFT, 
		 id1,  id2,  id3, 
		 id4,  id5,  id6,
		 id7,  id8,  id9, 
		 id10, id11, id12,
		 mList_->Scin, 
		 mList_->Al, 
		 mList_->Air, 
		 mList_->Air );
    
    G4VisAttributes *FiberLayerAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.0, 1.0, 1.0 ) );//Cyan
    G4VisAttributes *FiberFrameAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.3, 0.3, 0.3 ) );//Blac
    
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

  //Type A real SFT Square Fiber 1mm Single Cladding
  else if( DetType==1 ){
    ///////////////////// SFT(x,u,v,x,u,v, u,v,x,u,v,x)
    G4int id1=geomMan.GetDetectorId("SFT-x-1");
    G4int id2=geomMan.GetDetectorId("SFT-u-1");
    G4int id3=geomMan.GetDetectorId("SFT-v-1");
    G4int id4=geomMan.GetDetectorId("SFT-x-2");
    G4int id5=geomMan.GetDetectorId("SFT-u-2");
    G4int id6=geomMan.GetDetectorId("SFT-v-2");
    G4int id7=geomMan.GetDetectorId("SFT-u-3");
    G4int id8=geomMan.GetDetectorId("SFT-v-3");
    G4int id9=geomMan.GetDetectorId("SFT-x-3");
    G4int id10=geomMan.GetDetectorId("SFT-u-4");
    G4int id11=geomMan.GetDetectorId("SFT-v-4");
    G4int id12=geomMan.GetDetectorId("SFT-x-4");
    
    G4ThreeVector gPos1=geomMan.GetGlobalPosition(id1);
    G4ThreeVector gPos2=geomMan.GetGlobalPosition(id2);
    G4ThreeVector gPos3=geomMan.GetGlobalPosition(id3);
    G4ThreeVector gPos4=geomMan.GetGlobalPosition(id4);
    G4ThreeVector gPos5=geomMan.GetGlobalPosition(id5);
    G4ThreeVector gPos6=geomMan.GetGlobalPosition(id6);
    G4ThreeVector gPos7=geomMan.GetGlobalPosition(id7);
    G4ThreeVector gPos8=geomMan.GetGlobalPosition(id8);
    G4ThreeVector gPos9=geomMan.GetGlobalPosition(id9);
    G4ThreeVector gPos10=geomMan.GetGlobalPosition(id10);
    G4ThreeVector gPos11=geomMan.GetGlobalPosition(id11);
    G4ThreeVector gPos12=geomMan.GetGlobalPosition(id12);
    G4ThreeVector OfsLSFT( 0.0, 0.0, 0.0 );

    r_SFT_A *r_SF1 =
      new r_SFT_A( "SFT", pMother, RM, 
		   gPos1,  gPos2,  gPos3, 
		   gPos4,  gPos5,  gPos6,
		   gPos7,  gPos8,  gPos9, 
		   gPos10, gPos11, gPos12,
		   OfsLSFT, 
		   id1,  id2,  id3, 
		   id4,  id5,  id6,
		   id7,  id8,  id9, 
		   id10, id11, id12,
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
  
  //Type B real SFT Round Fiber Double Cladding
  else if( DetType==2 || DetType==3){
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

    //default 

    G4double LzLayer[rSFT_nLayer]={0.};
    
    if(DetType==2){
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
    }else if(DetType==3){//x,u,v layers are close to each other
      LzLayer[0]= -165.0*mm;//x
      LzLayer[1]= -164.0*mm;//u
      LzLayer[2]= -163.0*mm;//v
      LzLayer[3]= -46.0*mm;//x
      LzLayer[4]= -45.0*mm;//u
      LzLayer[5]= -44.0*mm;//v
      LzLayer[6]=  44.0*mm;//u
      LzLayer[7]=  45.0*mm;//v
      LzLayer[8]=  46.0*mm;//x
      LzLayer[9]= 163.0*mm;//u
      LzLayer[10] = 164.0*mm;//v
      LzLayer[11] =  165.0*mm;//x
    }

    r_SFT_B *r_SF1 =
      new r_SFT_B( "SFT", pMother, RM, 
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
    G4VisAttributes *FiberFrameAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.3, 0.3, 0.3 ) );//Blac
    
    // r_SF1->SetVisAttributes( FiberLayerAttrib, 
    // 			     FiberFrameAttrib, 
    // 			     &G4VisAttributes::Invisible,
    // 			     &G4VisAttributes::Invisible );

    r_SF1->SetVisAttributes( FiberLayerAttrib, 
           FiberInnerCradAttrib,
           FiberOuterCradAttrib,
			     &G4VisAttributes::Invisible,
			     &G4VisAttributes::Invisible,
			     &G4VisAttributes::Invisible );

    //////// Sensitive Detectors ////////
    G4SDManager *SDMan = G4SDManager::GetSDMpointer();
    SFTSD *sftSD = new SFTSD( "/spec/sft" );
    SDMan->AddNewDetector( sftSD );  
    r_SF1->SetSensitiveDetector( sftSD );
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
      new r_SFT_C( "SFT", pMother, RM, 
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
  ConfMan *confMan = ConfMan::GetConfManager();
  G4RotationMatrix RM;  

  ///////////////////////////////////////// 
  ////////////////////////////////T0
  ///////////////////////////////////////// 
    
  G4int idt0=geomMan.GetDetectorId("T0U");
  G4ThreeVector gPosT0=geomMan.GetGlobalPosition(idt0);
  G4ThreeVector OfsLT0( 0.0, 0.0, 0.0 );
  
  s_T0Wall *s_T0 =
    new s_T0Wall( "T0", pMother, RM, gPosT0, OfsLT0, 
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


void DetectorConstructionSpec::
MakeSpecMagnetandTarget( G4VPhysicalVolume *pMother )
{
  /////////////////////FM Magnet
  G4String nmPole="FMPole";
  G4String nmYoke="FMYoke";
  G4String nmCoil="FMCoil";
  G4RotationMatrix *zeRot= new G4RotationMatrix();
  zeRot->rotateY(0.*degree); 

  G4Material *mMag, *mAir;
  mMag= mList_->Fe;
  mAir= mList_->HeGas;
  G4int Id_=1;

  //Box
  G4Box* solid= new G4Box("FMM", 1500./2.*cm, 1500./2.*cm, 1500./2.*cm); 
  G4LogicalVolume *DetLV_= new G4LogicalVolume(solid, mAir, "FMM LV");
  DetLV_-> SetVisAttributes(G4VisAttributes::GetInvisible());

  //Yoke
  G4double dx_YokeTB=539.5*cm; 
  G4double dy_YokeTB= 83.0*cm; 
  G4double dz_YokeTB=212.0*cm; 

  G4Box* solidYokeTB= new G4Box(nmYoke, dx_YokeTB/2., dy_YokeTB/2., dz_YokeTB/2.); 
  G4LogicalVolume* LVYokeTB= new G4LogicalVolume(solidYokeTB, mMag, nmYoke);
  G4VisAttributes *va= new G4VisAttributes( G4Color( 0.0, 0.5, 1.0 ) ); //Blue
  LVYokeTB-> SetVisAttributes(va);
  new G4PVPlacement
    (zeRot, G4ThreeVector(-3.35*cm, 203.5*cm,0.*cm), LVYokeTB, nmYoke, DetLV_, false, Id_);
  new G4PVPlacement
    (zeRot, G4ThreeVector(-3.35*cm,-203.5*cm,0.*cm), LVYokeTB, nmYoke, DetLV_, false, Id_);
  
  G4Box* solidYokeRA= new G4Box(nmYoke, 84.4/2.*cm, 324./2.*cm, 212.0/2.*cm); 
  G4Box* solidYokeRB= new G4Box(nmYoke, 50.*cm, 324./2.*cm, 12.33*cm);
  G4RotationMatrix *rot0 = new G4RotationMatrix();
  rot0-> rotateY(-21.*degree);
  G4DisplacedSolid *solidYokeRC = new G4DisplacedSolid( nmYoke, solidYokeRB, rot0,
							G4ThreeVector(84.4/2.*cm,0.*cm,106.*cm));
  
  G4Box* solidYokeRD= new G4Box(nmYoke, 50.*cm, 324./2.*cm, 14.73*cm);
  G4RotationMatrix *rot1 = new G4RotationMatrix();
  rot1-> rotateY(47.74*degree);
  G4DisplacedSolid *solidYokeRE = new G4DisplacedSolid( nmYoke, solidYokeRD, rot1,
							G4ThreeVector(84.4/2.*cm,0.*cm,-106.*cm));
  G4SubtractionSolid* solidYokeR0=
    new G4SubtractionSolid(nmYoke,solidYokeRA,solidYokeRC);
  G4SubtractionSolid* solidYokeR1=
    new G4SubtractionSolid(nmYoke,solidYokeR0,solidYokeRE);

  G4LogicalVolume* LVYokeR= new G4LogicalVolume(solidYokeR1, mMag, nmYoke);
  LVYokeR-> SetVisAttributes(va);
  new G4PVPlacement(zeRot, G4ThreeVector(-237.2*cm,0.*cm,0.*cm), LVYokeR, nmYoke, DetLV_, false, Id_);
  
  G4Box* solidYokeLA= new G4Box(nmYoke, 91.1/2.*cm, 324./2.*cm, 212.0/2.*cm); 
  G4Box* solidYokeLB= new G4Box(nmYoke, 50.*cm, 324./2.*cm, 14.05*cm);
  G4RotationMatrix *rot2 = new G4RotationMatrix();
  rot2-> rotateY(39.70*degree);
  G4DisplacedSolid *solidYokeLC = new G4DisplacedSolid(nmYoke, solidYokeLB, rot2,
						       G4ThreeVector(-84.4/2.*cm,0.*cm,106.*cm));

  G4Box* solidYokeLD= new G4Box(nmYoke, 50.*cm, 324./2.*cm, 34.23*cm);
  G4RotationMatrix *rot3 = new G4RotationMatrix();
  rot3-> rotateY(-57.03*degree);
  G4DisplacedSolid *solidYokeLE = new G4DisplacedSolid(nmYoke, solidYokeLD, rot3,
						       G4ThreeVector(-84.4/2.*cm,0.*cm,-106.*cm));
  G4SubtractionSolid* solidYokeL0=
    new G4SubtractionSolid(nmYoke,solidYokeLA,solidYokeLC);
  G4SubtractionSolid* solidYokeL1=
    new G4SubtractionSolid(nmYoke,solidYokeL0,solidYokeLE);
  G4LogicalVolume* LVYokeL= new G4LogicalVolume(solidYokeL1, mMag, nmYoke);
  LVYokeL-> SetVisAttributes(va);
  new G4PVPlacement(zeRot, G4ThreeVector(240.55*cm,0.*cm,0.*cm), LVYokeL, nmYoke, DetLV_, false, Id_);
  
  //Pole
  G4Tubs *solidPole0 = new G4Tubs(nmPole,0.*cm,106.*cm,40./2.*cm,0,2.*3.14);
  G4LogicalVolume* LVPole0= new G4LogicalVolume(solidPole0, mMag, nmPole);
  LVPole0-> SetVisAttributes(va);
  G4RotationMatrix *rot4 = new G4RotationMatrix();
  rot4-> rotateX(-90.0*degree);
  new G4PVPlacement(rot4, G4ThreeVector(0.*cm,142.*cm,0.*cm), LVPole0, nmPole, DetLV_, false, Id_);
  new G4PVPlacement(rot4, G4ThreeVector(0.*cm,-142.*cm,0.*cm), LVPole0, nmPole, DetLV_, false, Id_);

  //Pole 1m version
  G4Tubs *solidPole1 = new G4Tubs(nmPole,0.*cm,88.*cm,72.0/2.*cm,0,2.*3.14);
  G4LogicalVolume* LVPole1= new G4LogicalVolume(solidPole1, mMag, nmPole);
  //va= new G4VisAttributes(G4Color(1.,0.5,0.));
  LVPole1-> SetVisAttributes(va);
  new G4PVPlacement(rot4, G4ThreeVector(0.*cm, 86.0*cm,0.*cm), LVPole1, nmPole, DetLV_, false, Id_);
  new G4PVPlacement(rot4, G4ThreeVector(0.*cm,-86.0*cm,0.*cm), LVPole1, nmPole, DetLV_, false, Id_);
  //   G4Tubs *solidPole1 = new G4Tubs(nmPole,0.*cm,88.*cm,30.6/2.*cm,0,2.*3.14);
  //   G4LogicalVolume* LVPole1= new G4LogicalVolume(solidPole1, mMag, nmPole);
  //   //va= new G4VisAttributes(G4Color(1.,0.5,0.));
  //   LVPole1-> SetVisAttributes(va);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,106.7*cm,0.*cm), LVPole1, nmPole, DetLV_, false, Id_);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,-106.7*cm,0.*cm), LVPole1, nmPole, DetLV_, false, Id_);

  G4Tubs *solidPole2 = new G4Tubs(nmPole,0.*cm,78.*cm,13.4/2.*cm,0,2.*3.14);
  G4LogicalVolume* LVPole2= new G4LogicalVolume(solidPole2, mMag, nmPole);
  LVPole2-> SetVisAttributes(va);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,84.7*cm,0.*cm), LVPole2, nmPole, DetLV_, false, Id_);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,-84.7*cm,0.*cm), LVPole2, nmPole, DetLV_, false, Id_);

  G4Paraboloid *solidPole3 = 
    new G4Paraboloid(nmPole,26./2.*cm,52.*cm,78.*cm);
  G4LogicalVolume* LVPole3= new G4LogicalVolume(solidPole3, mMag, nmPole);
  LVPole3-> SetVisAttributes(va);
  G4RotationMatrix *rot5 = new G4RotationMatrix();
  rot5-> rotateX(90.0*degree);
  //   new G4PVPlacement(rot5, G4ThreeVector(0.*cm,65.*cm,0.*cm), LVPole3, nmPole, DetLV_, false, Id_);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,-65.*cm,0.*cm), LVPole3, nmPole, DetLV_, false, Id_);

  G4Tubs *solidPole4 = new G4Tubs(nmPole,0.*cm,42.*cm,10.0/2.*cm,0,2.*3.14);
  G4LogicalVolume* LVPole4= new G4LogicalVolume(solidPole4, mMag, nmPole);
  LVPole4-> SetVisAttributes(va);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,47.0*cm,0.*cm), LVPole4, nmPole, DetLV_, false, Id_);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,-47.0*cm,0.*cm), LVPole4, nmPole, DetLV_, false, Id_);

  G4Paraboloid *solidPole5 = 
    new G4Paraboloid(nmPole,22./2.*cm,20.*cm,42.*cm);
  G4LogicalVolume* LVPole5= new G4LogicalVolume(solidPole5, mMag, nmPole);
  LVPole5-> SetVisAttributes(va);
  //  new G4PVPlacement(rot5, G4ThreeVector(0.*cm,31.*cm,0.*cm), LVPole5, nmPole, DetLV_, false, Id_);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,-31.*cm,0.*cm), LVPole5, nmPole, DetLV_, false, Id_);

  //Coil  
  G4Tubs *solidCoil = new G4Tubs(nmCoil,106.*cm,129.5*cm,30./2.*cm,0,2.*3.14);
  G4LogicalVolume* LVCoil= new G4LogicalVolume(solidCoil, mMag, nmCoil);
  va= new G4VisAttributes(G4Color(1.,0.,0.));
  LVCoil-> SetVisAttributes(va);
  new G4PVPlacement(rot5, G4ThreeVector(0.*cm,147.*cm,0.*cm), LVCoil, nmCoil, DetLV_, false, Id_);
  new G4PVPlacement(rot5, G4ThreeVector(0.*cm,-147.*cm,0.*cm), LVCoil, nmCoil, DetLV_, false, Id_);
  
  G4RotationMatrix rotMtx_; 
  rotMtx_.rotateX( -90.*degree);
  rotMtx_.rotateZ(  90.*degree);
  G4ThreeVector Pos( 0.0*mm, 0.0*mm, 0.0*mm );
  G4VPhysicalVolume *FMM = new G4PVPlacement( G4Transform3D(rotMtx_, Pos ), "FMM", DetLV_, pMother, true, Id_);

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  G4RotationMatrix RM;  
  RM.rotateX(  90.*degree );
  RM.rotateY( 180.*degree );

  //const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  ConfMan *confMan = ConfMan::GetConfManager();
  G4double tgtSizeX, tgtSizeY, tgtSizeZ; 
  G4int lnumTarg = geomMan.GetDetectorId("Target");
  G4ThreeVector targPos = geomMan.GetGlobalPosition( lnumTarg )*mm;
  G4double targAng2 = geomMan.GetRotAngle2( lnumTarg )*degree;
  G4double size=0.0;
  if( (size=confMan->GetTargetSizeX())>0.0 ) tgtSizeX=size;
  if( (size=confMan->GetTargetSizeY())>0.0 ) tgtSizeY=size;
  if( (size=confMan->GetTargetSizeZ())>0.0 ) tgtSizeZ=size;

  //LH2 cell
  G4Tubs *solidTGTBox =
    new G4Tubs( "TargetB", 0.0*mm, (tgtSizeX+50.0+1.0)/2.*mm, (tgtSizeZ+50.0+1.0)/2.*mm, 0.0*degree, 360.0*degree );
  G4Tubs *solidCFRP =
    new G4Tubs( "TargetC", (tgtSizeX+50.0-2.0)/2.*mm, (tgtSizeX+50.0)/2.*mm, (tgtSizeZ+50.0)/2.*mm, 
		0.0*degree, 360.0*degree );
  G4Tubs *solidTGT =
    new G4Tubs( "Target",  0.0*mm, tgtSizeX/2.*mm, tgtSizeZ/2.*mm, 0.0*degree, 360.0*degree );
  G4Tubs *solidMaylar =
    new G4Tubs( "TargetM", 0.0*mm, (tgtSizeX+50.0)/2.*mm, 0.3/2.*mm, 0.0*degree, 360.0*degree );
  
  G4LogicalVolume *logTGT, *logTGTC, *logTGTM, *logTGTB;
  if( confMan->ExistTarget() ){
    logTGT =
      new G4LogicalVolume( solidTGT, mList_->LiqH2, "Target", 0, 0, 0 );
    logTGTC =
      new G4LogicalVolume( solidCFRP, mList_->C, "TargetC", 0, 0, 0 );
    logTGTM =
      new G4LogicalVolume( solidMaylar, mList_->PET, "TargetM", 0, 0, 0 );
    logTGTB =
      new G4LogicalVolume( solidTGT, mList_->Air, "TargetB", 0, 0, 0 );
  }
  else{
    logTGT =
      new G4LogicalVolume( solidTGT, mList_->Air, "Target", 0, 0, 0 );
    logTGTC =
      new G4LogicalVolume( solidCFRP, mList_->Air, "TargetC", 0, 0, 0 );
    logTGTM =
      new G4LogicalVolume( solidMaylar, mList_->Air, "TargetM", 0, 0, 0 );
    logTGTB =
      new G4LogicalVolume( solidTGT, mList_->Air, "TargetB", 0, 0, 0 );
  }

  G4RotationMatrix RMTgt;  
  
  G4ThreeVector TgtPos( 0., 0.,  (-1.)*targPos.x()*mm );
  G4VPhysicalVolume *physTGTBox =
    new G4PVPlacement( G4Transform3D( RMTgt, TgtPos ),
		       logTGTB, "TargetBox", DetLV_, false, 0 );
  G4VPhysicalVolume *physTGTC = 
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
  		       "TargetC", logTGTC, physTGTBox, false, 0 );

  new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, (620.0+0.3)/2.*mm ),
  		     "TargetMu", logTGTM, physTGTC, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, (-1.)*(620.0+0.3)/2.*mm ),
  		     "TargetMd", logTGTM, physTGTC, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
		     "Target", logTGT, physTGTBox, false, 0 );

  G4VisAttributes *TgtAttrib =
    new G4VisAttributes( G4Colour::G4Colour( 0.0, 1.0, 1.0 ) );
  G4VisAttributes *TgtCAttrib =
    new G4VisAttributes( G4Colour::G4Colour( 0.0, 0.0, 0.0 ) );
  G4VisAttributes *TgtMAttrib =
    new G4VisAttributes( G4Colour::G4Colour( 1.0, 1.0, 1.0 ) );
  
  logTGT->SetVisAttributes(TgtAttrib);
  logTGTC->SetVisAttributes(TgtCAttrib);
  logTGTM->SetVisAttributes(TgtMAttrib);
  logTGTB->SetVisAttributes(&G4VisAttributes::Invisible);

  //////// Sensitive Detectors ////////
  //G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  TargetSD *tgtSD = new TargetSD( "/spec/target" );
  SDMan->AddNewDetector( tgtSD );
  logTGT->SetSensitiveDetector( tgtSD );



}
