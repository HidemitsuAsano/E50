/*
  DetectorConstructionSpec.hh

  2016/4  K.Shirotori
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
#include "r_SFT_A.hh"

#include "G4SDManager.hh"
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
  if( confMan->DetectorType()==0 ) DetType=0;
  if( confMan->DetectorType()==1 ) DetType=1;
  if( confMan->DetectorType()==2 ) DetType=2;

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

  //Type A real SFT Square Fiber Single Cladding
  if( DetType==1 ){
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
