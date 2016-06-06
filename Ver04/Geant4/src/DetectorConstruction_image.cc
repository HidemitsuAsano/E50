/*
  DetectorConstruction.cc

  2016/1  K.Shirotori
*/

#include "DetectorConstruction.hh"
#include "MaterialList.hh" 
#include "SpecField.hh"

#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"

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

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "s_DetectorSize.hh"
#include "s_T0Wall.hh"
#include "s_BFT.hh"
#include "s_SFT.hh"
#include "s_AFT.hh"
#include "s_ITWC1.hh"
#include "s_ITWC2L.hh"
#include "s_ITWC2R.hh"
#include "s_STWC1.hh"
#include "s_STWC2.hh"
#include "s_TofWall.hh"
#include "s_ITofWall.hh"
#include "s_ITofSWall.hh"
#include "s_PAD.hh"
#include "s_RICH.hh"
#include "s_PID1.hh"
#include "s_PID2.hh"
#include "s_MuonFilter.hh"

#include "r_T0Wall.hh"

#include "G4SDManager.hh"
#include "TargetSD.hh"
#include "T0SD.hh"
#include "BTrackerSD.hh"
#include "SITrackerSD.hh"
#include "AOTrackerSD.hh"
#include "IITrackerSD.hh"
#include "IOTrackerSD.hh"
#include "SOTrackerSD.hh"
#include "TofSD.hh"
#include "ITofSD.hh"
#include "PADSD.hh"
#include "VDSD.hh"
#include "RICHSD.hh"
#include "PID1SD.hh"
#include "PID2SD.hh"
#include "MuonFilterSD.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "ConfMan.hh"

DetectorConstruction::DetectorConstruction()
  : mList_(0), field_(0)
{
}

DetectorConstruction::~DetectorConstruction()
{
  delete field_;
  delete mList_;
}
MaterialList *DetectorConstruction::DefineMaterials()
{
  // not safe for exception 
  if(mList_) delete mList_;
  return new MaterialList();
}

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  mList_ = DefineMaterials();

  ConfMan *confMan = ConfMan::GetConfManager();

  if(confMan->ExistField()){
    field_=MakeMagFieldFromMap( confMan->FieldMapName(),
				confMan->MapScale() );
  }
  else{
    field_=MakeUniformMagField( 0.0*tesla );
  }

  G4VPhysicalVolume *world = ConstructPayload();

  return world;
}

G4MagneticField * DetectorConstruction::
MakeUniformMagField( G4double fieldValue )
{
  G4MagneticField *field=0;
  G4FieldManager* fieldMan
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  if( fieldValue!=0. ){ 
    field = new G4UniformMagField( G4ThreeVector( 0.,0.,fieldValue ));
    fieldMan->CreateChordFinder( field );
  }
  fieldMan->SetDetectorField( field );

  return field;
}

G4MagneticField * DetectorConstruction::
MakeMagFieldFromMap( const std::string & filename,
		     double NormFac )
{
  G4MagneticField *field = new SpecField( filename, NormFac );
  G4FieldManager *fieldManager =
    G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldManager->SetDetectorField( field );
  fieldManager->CreateChordFinder( field );
  fieldManager->GetChordFinder()->SetDeltaChord( 1.0E-3*mm );

  return field;
}

G4bool DetectorConstruction::IsVolumeStopper( G4VPhysicalVolume *physVol ) const
{
  // Check Magnet
  G4String name = physVol->GetName();
  if( name=="FMYoke" || name=="FMPole" || name=="FMCoil") return true;
  else
    return false;
}

//Main
void DetectorConstruction::
MakeSpecMagnet( G4VPhysicalVolume *pMother, G4Material *matGap )
{
  /////////////////////FM Magnet
  G4String nmYoke;
  G4String nmPole;
  G4String nmCoil;
  nmPole="FMPole";
  nmYoke="FMYoke";
  nmCoil="FMCoil";
  G4RotationMatrix *zeRot= new G4RotationMatrix();
  zeRot->rotateY(0.*degree); 

  G4Material *mMag_, *mAir_;
  mMag_= mList_->Fe;
  mAir_= mList_->HeGas;
  G4int Id_=1;

  //Box
  G4Box* solid= new G4Box("FMM", 1500./2.*cm, 1500./2.*cm, 1500./2.*cm); 
  G4LogicalVolume *DetLV_= new G4LogicalVolume(solid, mAir_, "FMM LV");
  DetLV_-> SetVisAttributes(G4VisAttributes::GetInvisible());

  //Yoke
  G4double dx_YokeTB=539.5*cm; 
  G4double dy_YokeTB= 83.0*cm; 
  G4double dz_YokeTB=212.0*cm; 

  G4Box* solidYokeTB= new G4Box(nmYoke, dx_YokeTB/2., dy_YokeTB/2., dz_YokeTB/2.); 
  G4LogicalVolume* LVYokeTB= new G4LogicalVolume(solidYokeTB, mMag_, nmYoke);
  G4VisAttributes *va= new G4VisAttributes( G4Color( 0.0, 0.5, 1.0 ) ); //Blue
  LVYokeTB-> SetVisAttributes(va);
  new G4PVPlacement
    (zeRot, G4ThreeVector(-3.35*cm, 203.5*cm,0.*cm), LVYokeTB, nmYoke, DetLV_, false, Id_);
  // new G4PVPlacement
  //   (zeRot, G4ThreeVector(-3.35*cm,-203.5*cm,0.*cm), LVYokeTB, nmYoke, DetLV_, false, Id_);
  
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

  G4LogicalVolume* LVYokeR= new G4LogicalVolume(solidYokeR1, mMag_, nmYoke);
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
  G4LogicalVolume* LVYokeL= new G4LogicalVolume(solidYokeL1, mMag_, nmYoke);
  LVYokeL-> SetVisAttributes(va);
  //new G4PVPlacement(zeRot, G4ThreeVector(240.55*cm,0.*cm,0.*cm), LVYokeL, nmYoke, DetLV_, false, Id_);
  
  //Pole
  G4Tubs *solidPole0 = new G4Tubs(nmPole,0.*cm,106.*cm,40./2.*cm,0,2.*3.14);
  G4LogicalVolume* LVPole0= new G4LogicalVolume(solidPole0, mMag_, nmPole);
  LVPole0-> SetVisAttributes(va);
  G4RotationMatrix *rot4 = new G4RotationMatrix();
  rot4-> rotateX(-90.0*degree);
  new G4PVPlacement(rot4, G4ThreeVector(0.*cm,142.*cm,0.*cm), LVPole0, nmPole, DetLV_, false, Id_);
  //new G4PVPlacement(rot4, G4ThreeVector(0.*cm,-142.*cm,0.*cm), LVPole0, nmPole, DetLV_, false, Id_);

  //Pole 1m version
  G4Tubs *solidPole1 = new G4Tubs(nmPole,0.*cm,88.*cm,72.0/2.*cm,0,2.*3.14);
  G4LogicalVolume* LVPole1= new G4LogicalVolume(solidPole1, mMag_, nmPole);
  //va= new G4VisAttributes(G4Color(1.,0.5,0.));
  LVPole1-> SetVisAttributes(va);
  new G4PVPlacement(rot4, G4ThreeVector(0.*cm, 86.0*cm,0.*cm), LVPole1, nmPole, DetLV_, false, Id_);
  //new G4PVPlacement(rot4, G4ThreeVector(0.*cm,-86.0*cm,0.*cm), LVPole1, nmPole, DetLV_, false, Id_);
  //   G4Tubs *solidPole1 = new G4Tubs(nmPole,0.*cm,88.*cm,30.6/2.*cm,0,2.*3.14);
  //   G4LogicalVolume* LVPole1= new G4LogicalVolume(solidPole1, mMag_, nmPole);
  //   //va= new G4VisAttributes(G4Color(1.,0.5,0.));
  //   LVPole1-> SetVisAttributes(va);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,106.7*cm,0.*cm), LVPole1, nmPole, DetLV_, false, Id_);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,-106.7*cm,0.*cm), LVPole1, nmPole, DetLV_, false, Id_);

  G4Tubs *solidPole2 = new G4Tubs(nmPole,0.*cm,78.*cm,13.4/2.*cm,0,2.*3.14);
  G4LogicalVolume* LVPole2= new G4LogicalVolume(solidPole2, mMag_, nmPole);
  LVPole2-> SetVisAttributes(va);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,84.7*cm,0.*cm), LVPole2, nmPole, DetLV_, false, Id_);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,-84.7*cm,0.*cm), LVPole2, nmPole, DetLV_, false, Id_);

  G4Paraboloid *solidPole3 = 
    new G4Paraboloid(nmPole,26./2.*cm,52.*cm,78.*cm);
  G4LogicalVolume* LVPole3= new G4LogicalVolume(solidPole3, mMag_, nmPole);
  LVPole3-> SetVisAttributes(va);
  G4RotationMatrix *rot5 = new G4RotationMatrix();
  rot5-> rotateX(90.0*degree);
  //   new G4PVPlacement(rot5, G4ThreeVector(0.*cm,65.*cm,0.*cm), LVPole3, nmPole, DetLV_, false, Id_);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,-65.*cm,0.*cm), LVPole3, nmPole, DetLV_, false, Id_);

  G4Tubs *solidPole4 = new G4Tubs(nmPole,0.*cm,42.*cm,10.0/2.*cm,0,2.*3.14);
  G4LogicalVolume* LVPole4= new G4LogicalVolume(solidPole4, mMag_, nmPole);
  LVPole4-> SetVisAttributes(va);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,47.0*cm,0.*cm), LVPole4, nmPole, DetLV_, false, Id_);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,-47.0*cm,0.*cm), LVPole4, nmPole, DetLV_, false, Id_);

  G4Paraboloid *solidPole5 = 
    new G4Paraboloid(nmPole,22./2.*cm,20.*cm,42.*cm);
  G4LogicalVolume* LVPole5= new G4LogicalVolume(solidPole5, mMag_, nmPole);
  LVPole5-> SetVisAttributes(va);
  //  new G4PVPlacement(rot5, G4ThreeVector(0.*cm,31.*cm,0.*cm), LVPole5, nmPole, DetLV_, false, Id_);
  //   new G4PVPlacement(rot4, G4ThreeVector(0.*cm,-31.*cm,0.*cm), LVPole5, nmPole, DetLV_, false, Id_);

  //Coil  
  G4Tubs *solidCoil = new G4Tubs(nmCoil,106.*cm,129.5*cm,30./2.*cm,0,2.*3.14);
  G4LogicalVolume* LVCoil= new G4LogicalVolume(solidCoil, mMag_, nmCoil);
  va= new G4VisAttributes(G4Color(1.,0.,0.));
  LVCoil-> SetVisAttributes(va);
  new G4PVPlacement(rot5, G4ThreeVector(0.*cm,147.*cm,0.*cm), LVCoil, nmCoil, DetLV_, false, Id_);
  //new G4PVPlacement(rot5, G4ThreeVector(0.*cm,-147.*cm,0.*cm), LVCoil, nmCoil, DetLV_, false, Id_);
  
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


  ////////////////////////////////////////
  /////////////////////////////////Target
  /////////////////////////////////////// 

  G4double tgtSizeX, tgtSizeY, tgtSizeZ; 
  G4int lnumTarg = geomMan.GetDetectorId("Target");
  G4ThreeVector targPos = geomMan.GetGlobalPosition( lnumTarg )*mm;
  G4double targAng2 = geomMan.GetRotAngle2( lnumTarg )*degree;
  ConfMan *confMan = ConfMan::GetConfManager();
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
  
  G4LogicalVolume *logTGT =
    new G4LogicalVolume( solidTGT, mList_->LiqH2, "Target", 0, 0, 0 );
  G4LogicalVolume *logTGTC =
    new G4LogicalVolume( solidCFRP, mList_->C, "TargetC", 0, 0, 0 );
  G4LogicalVolume *logTGTM =
    new G4LogicalVolume( solidMaylar, mList_->PET, "TargetM", 0, 0, 0 );
  G4LogicalVolume *logTGTB =
    new G4LogicalVolume( solidTGT, mList_->Air, "TargetB", 0, 0, 0 );

  G4RotationMatrix RMTgt;  
  
  G4ThreeVector TgtPos( 0., 0.,  (-1.)*targPos.x()*mm );
  G4VPhysicalVolume *physTGTBox =
    new G4PVPlacement( G4Transform3D( RMTgt, TgtPos ),
		       logTGTB, "TargetBox", DetLV_, false, 0 );
  // G4VPhysicalVolume *physTGTC = 
  //   new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
  // 		       "TargetC", logTGTC, physTGTBox, false, 0 );

  // new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, (620.0+0.3)/2.*mm ),
  // 		     "TargetMu", logTGTM, physTGTC, false, 0 );
  // new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, (-1.)*(620.0+0.3)/2.*mm ),
  // 		     "TargetMd", logTGTM, physTGTC, false, 0 );
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
  TargetSD *tgtSD = new TargetSD( "/spec/target" );
  SDMan->AddNewDetector( tgtSD );
  logTGT->SetSensitiveDetector( tgtSD );


  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////

  G4bool s_DetFlag=false;
  G4bool r_DetFlag=false;
  if( !confMan->DetectorType() ) s_DetFlag=true;
  if( confMan->DetectorType() )  r_DetFlag=true;

  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////
  //Simple detectors
  if( s_DetFlag ){

    ///////////////////////////////////////// 
    ////////////////////////////////T0
    ///////////////////////////////////////// 
    
    G4int idt0=geomMan.GetDetectorId("T0");
    G4ThreeVector gPosT0=geomMan.GetGlobalPosition(idt0);
    G4ThreeVector OfsLT0( 0.0, 0.0, 0.0 );
    
    s_T0Wall *s_T0 =
      new s_T0Wall( "T0", DetLV_, RM, gPosT0, OfsLT0, 
    		    mList_->Scin, 
		    mList_->Air );
    
    G4VisAttributes *T0ScinAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 1.0, 0.0, 1.0 ) );//Magenta

    s_T0->SetVisAttributes( T0ScinAttrib, 
			    &G4VisAttributes::Invisible );
  
    //////// Sensitive Detectors ////////
    T0SD *t0SD = new T0SD( "/spec/t0" );
    SDMan->AddNewDetector( t0SD );
    s_T0->SetSensitiveDetector( t0SD );


    //////////////////////////////////////////
    ////////////////////////////FiberTrackers
    ///////////////////////////////////////// 
    
    ///////////////////// BT1(x,u,v,x,u,v)
    G4int id101=geomMan.GetDetectorId("BT1-x-1");
    G4int id102=geomMan.GetDetectorId("BT1-u-1");
    G4int id103=geomMan.GetDetectorId("BT1-v-1");
    G4int id104=geomMan.GetDetectorId("BT1-x-2");
    G4int id105=geomMan.GetDetectorId("BT1-u-2");
    G4int id106=geomMan.GetDetectorId("BT1-v-2");
    
    G4ThreeVector gPos101=geomMan.GetGlobalPosition(id101);
    G4ThreeVector gPos102=geomMan.GetGlobalPosition(id102);
    G4ThreeVector gPos103=geomMan.GetGlobalPosition(id103);
    G4ThreeVector gPos104=geomMan.GetGlobalPosition(id104);
    G4ThreeVector gPos105=geomMan.GetGlobalPosition(id105);
    G4ThreeVector gPos106=geomMan.GetGlobalPosition(id106);
    G4ThreeVector OfsLBT1( 0.0, 0.0, 0.0 );
    
    s_BFT *s_BT1 =
      new s_BFT( "BT1", DetLV_, RM, 
    		 gPos101, gPos102, gPos103, 
    		 gPos104, gPos105, gPos106, 
    		 OfsLBT1, 
    		 id101, id102, id103, 
    		 id104, id105, id106, 
    		 mList_->Scin, 
		 mList_->Al, 
    		 mList_->Air, 
		 mList_->Air );
  
    ///////////////////// BT2(u,v,x,u,v,x)
    G4int id107=geomMan.GetDetectorId("BT2-u-1");
    G4int id108=geomMan.GetDetectorId("BT2-v-1");
    G4int id109=geomMan.GetDetectorId("BT2-x-1");
    G4int id110=geomMan.GetDetectorId("BT2-u-2");
    G4int id111=geomMan.GetDetectorId("BT2-v-2");
    G4int id112=geomMan.GetDetectorId("BT2-x-2");

    G4ThreeVector gPos107=geomMan.GetGlobalPosition(id107);
    G4ThreeVector gPos108=geomMan.GetGlobalPosition(id108);
    G4ThreeVector gPos109=geomMan.GetGlobalPosition(id109);
    G4ThreeVector gPos110=geomMan.GetGlobalPosition(id110);
    G4ThreeVector gPos111=geomMan.GetGlobalPosition(id111);
    G4ThreeVector gPos112=geomMan.GetGlobalPosition(id112);
    G4ThreeVector OfsLBT2( 0.0, 0.0, 0.0 );
    
    s_BFT *s_BT2 =
      new s_BFT( "BT2", DetLV_, RM, 
    		 gPos107, gPos108, gPos109, 
    		 gPos110, gPos111, gPos112, 
    		 OfsLBT2, 
    		 id107, id108, id109, 
    		 id110, id111, id112, 
    		 mList_->Scin, 
		 mList_->Al, 
    		 mList_->Air, 
		 mList_->Air );
  
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
      new s_SFT( "SF1", DetLV_, RM, 
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

    ///////////////////// AFTR(x,x)
    G4int id51=geomMan.GetDetectorId("AFTR-x-1");
    G4int id52=geomMan.GetDetectorId("AFTR-x-2");
    
    G4ThreeVector gPos51=geomMan.GetGlobalPosition(id51);
    G4ThreeVector gPos52=geomMan.GetGlobalPosition(id52);
    G4ThreeVector OfsLAFTR( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMAFTR;  
    RMAFTR.rotateX( 90.*degree );
    RMAFTR.rotateY( AFT_RotAngleR );
    
    s_AFT *s_AF1R =
      new s_AFT( "AF1R", DetLV_, RMAFTR, 
    		 gPos51,  
		 gPos52, 
    		 OfsLAFTR, 
    		 id51,  
		 id52, 
    		 mList_->Scin, 
		 mList_->C, 
    		 mList_->Air, 
		 mList_->Air );

    ///////////////////// AFTL(x,x)
    G4int id53=geomMan.GetDetectorId("AFTL-x-1");
    G4int id54=geomMan.GetDetectorId("AFTL-x-2");

    G4ThreeVector gPos53=geomMan.GetGlobalPosition(id53);
    G4ThreeVector gPos54=geomMan.GetGlobalPosition(id54);
    G4ThreeVector OfsLAFTL( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMAFTL;  
    RMAFTL.rotateX( 90.*degree );
    RMAFTL.rotateY( AFT_RotAngleL );
    
    s_AFT *s_AF1L =
      new s_AFT( "AF1L", DetLV_, RMAFTL, 
    		 gPos53,  
		 gPos54, 
    		 OfsLAFTL, 
    		 id53,  
		 id54, 
    		 mList_->Scin, 
		 mList_->C, 
    		 mList_->Air, 
		 mList_->Air );
    
    G4VisAttributes *FiberLayerAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.0, 1.0, 1.0 ) );//Cyan
    G4VisAttributes *FiberFrameAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.3, 0.3, 0.3 ) );//Black

    // s_BT1->SetVisAttributes( FiberLayerAttrib, 
    // 			     FiberFrameAttrib, 
    // 			     &G4VisAttributes::Invisible,
    // 			     &G4VisAttributes::Invisible );
    // s_BT2->SetVisAttributes( FiberLayerAttrib, 
    // 			     FiberFrameAttrib, 
    // 			     &G4VisAttributes::Invisible,
    // 			     &G4VisAttributes::Invisible );
    // s_SF1->SetVisAttributes( FiberLayerAttrib, 
    // 			     FiberFrameAttrib, 
    // 			     &G4VisAttributes::Invisible,
    // 			     &G4VisAttributes::Invisible );
    // s_AF1R->SetVisAttributes( FiberLayerAttrib, 
    // 			      FiberFrameAttrib, 
    // 			      &G4VisAttributes::Invisible,
    // 			      &G4VisAttributes::Invisible );
    // s_AF1L->SetVisAttributes( FiberLayerAttrib, 
    // 			      FiberFrameAttrib, 
    // 			      &G4VisAttributes::Invisible,
    // 			      &G4VisAttributes::Invisible );

    s_BT1->SetVisAttributes( FiberLayerAttrib, 
    			     &G4VisAttributes::Invisible,
    			     &G4VisAttributes::Invisible,
    			     &G4VisAttributes::Invisible );
    s_BT2->SetVisAttributes( FiberLayerAttrib, 
    			     &G4VisAttributes::Invisible,
    			     &G4VisAttributes::Invisible,
    			     &G4VisAttributes::Invisible );
    s_SF1->SetVisAttributes( FiberLayerAttrib, 
    			     &G4VisAttributes::Invisible,
    			     &G4VisAttributes::Invisible,
    			     &G4VisAttributes::Invisible );
    s_AF1R->SetVisAttributes( FiberLayerAttrib, 
    			      &G4VisAttributes::Invisible,
    			      &G4VisAttributes::Invisible,
    			      &G4VisAttributes::Invisible );
    s_AF1L->SetVisAttributes( FiberLayerAttrib, 
    			      &G4VisAttributes::Invisible,
    			      &G4VisAttributes::Invisible,
    			      &G4VisAttributes::Invisible );

    //////// Sensitive Detectors ////////
    BTrackerSD  *btrackerSD  = new BTrackerSD( "/spec/btr" );
    SITrackerSD *sitrackerSD = new SITrackerSD( "/spec/sitr" );
    AOTrackerSD *aotrackerSD = new AOTrackerSD( "/spec/aotr" );
    SDMan->AddNewDetector( btrackerSD );  
    SDMan->AddNewDetector( sitrackerSD );  
    SDMan->AddNewDetector( aotrackerSD );  

    s_BT1->SetSensitiveDetector( btrackerSD );
    s_BT2->SetSensitiveDetector( btrackerSD );
    s_SF1->SetSensitiveDetector( sitrackerSD );
    s_AF1R->SetSensitiveDetector( aotrackerSD );
    s_AF1L->SetSensitiveDetector( aotrackerSD );


    //////////////////////////////////////////
    ////////////////////////////Wire Chamber
    ///////////////////////////////////////// 

    ///////////////////// IT1(x,u,v,x,u,v,u,v,x,u,v,x) 
    G4int id15 =geomMan.GetDetectorId("IT1-x-1");
    G4int id16 =geomMan.GetDetectorId("IT1-u-1");
    G4int id17 =geomMan.GetDetectorId("IT1-v-1");
    G4int id18 =geomMan.GetDetectorId("IT1-x-2");
    G4int id19 =geomMan.GetDetectorId("IT1-u-2");
    G4int id20 =geomMan.GetDetectorId("IT1-v-2");
    G4int id21 =geomMan.GetDetectorId("IT1-u-3");
    G4int id22 =geomMan.GetDetectorId("IT1-v-3");
    G4int id23 =geomMan.GetDetectorId("IT1-x-3");
    G4int id24 =geomMan.GetDetectorId("IT1-u-4");
    G4int id25 =geomMan.GetDetectorId("IT1-v-4");
    G4int id26 =geomMan.GetDetectorId("IT1-x-4");
    G4ThreeVector gPos15 =geomMan.GetGlobalPosition(id15);
    G4ThreeVector gPos16 =geomMan.GetGlobalPosition(id16);
    G4ThreeVector gPos17 =geomMan.GetGlobalPosition(id17);
    G4ThreeVector gPos18 =geomMan.GetGlobalPosition(id18);
    G4ThreeVector gPos19 =geomMan.GetGlobalPosition(id19);
    G4ThreeVector gPos20 =geomMan.GetGlobalPosition(id20);
    G4ThreeVector gPos21 =geomMan.GetGlobalPosition(id21);
    G4ThreeVector gPos22 =geomMan.GetGlobalPosition(id22);
    G4ThreeVector gPos23 =geomMan.GetGlobalPosition(id23);
    G4ThreeVector gPos24 =geomMan.GetGlobalPosition(id24);
    G4ThreeVector gPos25 =geomMan.GetGlobalPosition(id25);
    G4ThreeVector gPos26 =geomMan.GetGlobalPosition(id26);
    G4ThreeVector OfsLIT1(  0.0, 0.0, 0.0 );
    
    s_ITWC1 *s_IT1 =
      new s_ITWC1( "IT1", DetLV_, RM, 
    		   gPos15, gPos16, gPos17, 
    		   gPos18, gPos19, gPos20, 
    		   gPos21, gPos22, gPos23, 
    		   gPos24, gPos25, gPos26, 
    		   OfsLIT1, 
    		   id15, id16, id17, 
    		   id18, id19, id20, 
    		   id21, id22, id23, 
    		   id24, id25, id26,
    		   mList_->W, 
		   mList_->Al,
    		   mList_->Ar80IsoButane20Gas, 
    		   mList_->Air );

    ///////////////////// IT2R(x,v,u,v,u,x,v,u,x) 
    G4int id31 =geomMan.GetDetectorId("IT2R-x-1");
    G4int id32 =geomMan.GetDetectorId("IT2R-v-1");
    G4int id33 =geomMan.GetDetectorId("IT2R-u-1");
    G4int id34 =geomMan.GetDetectorId("IT2R-v-2");
    G4int id35 =geomMan.GetDetectorId("IT2R-u-2");
    G4int id36 =geomMan.GetDetectorId("IT2R-x-2");
    G4int id37 =geomMan.GetDetectorId("IT2R-v-3");
    G4int id38 =geomMan.GetDetectorId("IT2R-u-3");
    G4int id39 =geomMan.GetDetectorId("IT2R-x-3");
    G4ThreeVector gPos31 =geomMan.GetGlobalPosition(id31);
    G4ThreeVector gPos32 =geomMan.GetGlobalPosition(id32);
    G4ThreeVector gPos33 =geomMan.GetGlobalPosition(id33);
    G4ThreeVector gPos34 =geomMan.GetGlobalPosition(id34);
    G4ThreeVector gPos35 =geomMan.GetGlobalPosition(id35);
    G4ThreeVector gPos36 =geomMan.GetGlobalPosition(id36);
    G4ThreeVector gPos37 =geomMan.GetGlobalPosition(id37);
    G4ThreeVector gPos38 =geomMan.GetGlobalPosition(id38);
    G4ThreeVector gPos39 =geomMan.GetGlobalPosition(id39);
    G4ThreeVector OfsLIT2R( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMIT2R;  
    RMIT2R.rotateX( 90.*degree );
    RMIT2R.rotateY( IT2_RotAngleR );

    s_ITWC2R *s_IT2R =
      new s_ITWC2R( "IT2R", DetLV_, RMIT2R, 
    		    gPos31, gPos32, gPos33, 
    		    gPos34, gPos35, gPos36, 
    		    gPos37, gPos38, gPos39,
    		    OfsLIT2R, 
    		    id31, id32, id33, 
    		    id34, id35, id36, 
    		    id37, id38, id39, 
    		    mList_->W, 
    		    mList_->Al,
    		    mList_->Ar80IsoButane20Gas, 
    		    mList_->Air );

    ///////////////////// IT2L(x,u,v,u,v,x,u,v,x) 
    G4int id41 =geomMan.GetDetectorId("IT2L-x-1");
    G4int id42 =geomMan.GetDetectorId("IT2L-u-1");
    G4int id43 =geomMan.GetDetectorId("IT2L-v-1");
    G4int id44 =geomMan.GetDetectorId("IT2L-u-2");
    G4int id45 =geomMan.GetDetectorId("IT2L-v-2");
    G4int id46 =geomMan.GetDetectorId("IT2L-x-2");
    G4int id47 =geomMan.GetDetectorId("IT2L-u-3");
    G4int id48 =geomMan.GetDetectorId("IT2L-v-3");
    G4int id49 =geomMan.GetDetectorId("IT2L-x-3");
    G4ThreeVector gPos41 =geomMan.GetGlobalPosition(id41);
    G4ThreeVector gPos42 =geomMan.GetGlobalPosition(id42);
    G4ThreeVector gPos43 =geomMan.GetGlobalPosition(id43);
    G4ThreeVector gPos44 =geomMan.GetGlobalPosition(id44);
    G4ThreeVector gPos45 =geomMan.GetGlobalPosition(id45);
    G4ThreeVector gPos46 =geomMan.GetGlobalPosition(id46);
    G4ThreeVector gPos47 =geomMan.GetGlobalPosition(id47);
    G4ThreeVector gPos48 =geomMan.GetGlobalPosition(id48);
    G4ThreeVector gPos49 =geomMan.GetGlobalPosition(id49);
    G4ThreeVector OfsLIT2L( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMIT2L;  
    RMIT2L.rotateX( 90.*degree );
    RMIT2L.rotateY( IT2_RotAngleL );

    s_ITWC2L *s_IT2L =
      new s_ITWC2L( "IT2L", DetLV_, RMIT2L, 
    		    gPos41, gPos42, gPos43, 
    		    gPos44, gPos45, gPos46, 
    		    gPos47, gPos48, gPos49,
    		    OfsLIT2L, 
    		    id41, id42, id43, 
    		    id44, id45, id46, 
    		    id47, id48, id49, 
    		    mList_->W, 
    		    mList_->Al,
    		    mList_->Ar80IsoButane20Gas, 
    		    mList_->Air );


    ///////////////////// ST1(x,u,v,x,u,v,u,v,x)
    G4int id61 =geomMan.GetDetectorId("ST1-x-1");
    G4int id62 =geomMan.GetDetectorId("ST1-u-1");
    G4int id63 =geomMan.GetDetectorId("ST1-v-1");
    G4int id64 =geomMan.GetDetectorId("ST1-x-2");
    G4int id65 =geomMan.GetDetectorId("ST1-u-2");
    G4int id66 =geomMan.GetDetectorId("ST1-v-2");
    G4int id67 =geomMan.GetDetectorId("ST1-u-3");
    G4int id68 =geomMan.GetDetectorId("ST1-v-3");
    G4int id69 =geomMan.GetDetectorId("ST1-x-3");
    G4ThreeVector gPos61 =geomMan.GetGlobalPosition(id61);
    G4ThreeVector gPos62 =geomMan.GetGlobalPosition(id62);
    G4ThreeVector gPos63 =geomMan.GetGlobalPosition(id63);
    G4ThreeVector gPos64 =geomMan.GetGlobalPosition(id64);
    G4ThreeVector gPos65 =geomMan.GetGlobalPosition(id65);
    G4ThreeVector gPos66 =geomMan.GetGlobalPosition(id66);
    G4ThreeVector gPos67 =geomMan.GetGlobalPosition(id67);
    G4ThreeVector gPos68 =geomMan.GetGlobalPosition(id68);
    G4ThreeVector gPos69 =geomMan.GetGlobalPosition(id69);
    G4ThreeVector OfsLST1(  0.0, 0.0, 0.0 );
    
    s_STWC1 *s_ST1 =
      new s_STWC1( "ST1", DetLV_, RM, 
    		   gPos61, gPos62, gPos63, 
    		   gPos64, gPos65, gPos66, 
    		   gPos67, gPos68, gPos69, 
    		   OfsLST1, 
    		   id61, id62, id63, 
    		   id64, id65, id66, 
    		   id67, id68, id69, 
    		   mList_->W, 
    		   mList_->Al,
    		   mList_->Ar80IsoButane20Gas, 
    		   mList_->Air );

    ///////////////////// ST2(x,u,v,u,v,x,u,v,x)
    G4int id71 =geomMan.GetDetectorId("ST2-x-1");
    G4int id72 =geomMan.GetDetectorId("ST2-u-1");
    G4int id73 =geomMan.GetDetectorId("ST2-v-1");
    G4int id74 =geomMan.GetDetectorId("ST2-u-2");
    G4int id75 =geomMan.GetDetectorId("ST2-v-2");
    G4int id76 =geomMan.GetDetectorId("ST2-x-2");
    G4int id77 =geomMan.GetDetectorId("ST2-u-3");
    G4int id78 =geomMan.GetDetectorId("ST2-v-3");
    G4int id79 =geomMan.GetDetectorId("ST2-x-3");
    G4ThreeVector gPos71 =geomMan.GetGlobalPosition(id71);
    G4ThreeVector gPos72 =geomMan.GetGlobalPosition(id72);
    G4ThreeVector gPos73 =geomMan.GetGlobalPosition(id73);
    G4ThreeVector gPos74 =geomMan.GetGlobalPosition(id74);
    G4ThreeVector gPos75 =geomMan.GetGlobalPosition(id75);
    G4ThreeVector gPos76 =geomMan.GetGlobalPosition(id76);
    G4ThreeVector gPos77 =geomMan.GetGlobalPosition(id77);
    G4ThreeVector gPos78 =geomMan.GetGlobalPosition(id78);
    G4ThreeVector gPos79 =geomMan.GetGlobalPosition(id79);
    G4ThreeVector OfsLST2(  0.0, 0.0, 0.0 );
    
    s_STWC2 *s_ST2 =
      new s_STWC2( "ST2", DetLV_, RM, 
    		   gPos71, gPos72, gPos73, 
    		   gPos74, gPos75, gPos76, 
    		   gPos77, gPos78, gPos79, 
    		   OfsLST2, 
    		   id71, id72, id73, 
    		   id74, id75, id76, 
    		   id77, id78, id79, 
    		   mList_->W, 
    		   mList_->Al,
    		   mList_->Ar80IsoButane20Gas, 
    		   mList_->Air );

    G4VisAttributes *WireAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 1.0, 1.0, 0.0 ) );//Yellow
    G4VisAttributes *WCFrameAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.7, 0.7, 0.7 ) );//Gray
    
    // s_IT1->SetVisAttributes( WireAttrib, 
    // 			     WCFrameAttrib,
    // 			     &G4VisAttributes::Invisible,
    // 			     &G4VisAttributes::Invisible );
    // s_IT2R->SetVisAttributes( WireAttrib, 
    // 			      WCFrameAttrib,
    // 			      &G4VisAttributes::Invisible,
    // 			      &G4VisAttributes::Invisible );
    // s_IT2L->SetVisAttributes( WireAttrib, 
    // 			      WCFrameAttrib,
    // 			      &G4VisAttributes::Invisible,
    // 			      &G4VisAttributes::Invisible );
    // s_ST1->SetVisAttributes( WireAttrib, 
    // 			     WCFrameAttrib,
    // 			     &G4VisAttributes::Invisible,
    // 			     &G4VisAttributes::Invisible );
    // s_ST2->SetVisAttributes( WireAttrib, 
    // 			     WCFrameAttrib,
    // 			     &G4VisAttributes::Invisible,
    // 			     &G4VisAttributes::Invisible );

    s_IT1->SetVisAttributes( WireAttrib,
    			     &G4VisAttributes::Invisible,
    			     &G4VisAttributes::Invisible,
    			     &G4VisAttributes::Invisible );
    s_IT2R->SetVisAttributes( WireAttrib,
    			      &G4VisAttributes::Invisible,
    			      &G4VisAttributes::Invisible,
    			      &G4VisAttributes::Invisible );
    s_IT2L->SetVisAttributes( WireAttrib,
    			      &G4VisAttributes::Invisible,
    			      &G4VisAttributes::Invisible,
    			      &G4VisAttributes::Invisible );
    s_ST1->SetVisAttributes( WireAttrib,
    			     &G4VisAttributes::Invisible,
    			     &G4VisAttributes::Invisible,
    			     &G4VisAttributes::Invisible );
    s_ST2->SetVisAttributes( WireAttrib,
    			     &G4VisAttributes::Invisible,
    			     &G4VisAttributes::Invisible,
    			     &G4VisAttributes::Invisible );
    
    //////// Sensitive Detectors ////////
    IITrackerSD *iitrackerSD = new IITrackerSD( "/spec/iitr" );
    IOTrackerSD *iotrackerSD = new IOTrackerSD( "/spec/iotr" );
    SOTrackerSD *sotrackerSD = new SOTrackerSD( "/spec/sotr" );
    SDMan->AddNewDetector( iitrackerSD );  
    SDMan->AddNewDetector( iotrackerSD );  
    SDMan->AddNewDetector( sotrackerSD );  

    s_IT1->SetSensitiveDetector( iitrackerSD );
    s_IT2R->SetSensitiveDetector( iotrackerSD );
    s_IT2L->SetSensitiveDetector( iotrackerSD );
    s_ST1->SetSensitiveDetector( sotrackerSD );
    s_ST2->SetSensitiveDetector( sotrackerSD );

    
    ///////////////////////////////////////// 
    ////////////////////////////////TOF
    ///////////////////////////////////////// 

    G4int idtof=geomMan.GetDetectorId("TOF");
    G4ThreeVector gPosTof=geomMan.GetGlobalPosition(idtof);
    G4ThreeVector OfsLTof( 0.0, 0.0, 0.0 );

    s_TofWall *s_Tof =
      new s_TofWall( "TOF", DetLV_, RM, gPosTof, OfsLTof, 
		     mList_->Scin, 
		     mList_->Air );
    
    G4VisAttributes *TofScinAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 1.0, 0.0, 1.0 ) );//Magenta

    s_Tof->SetVisAttributes( TofScinAttrib, &G4VisAttributes::Invisible );
  
    //////// Sensitive Detectors ////////
    TofSD *tofSD = new TofSD( "/spec/tof" );
    SDMan->AddNewDetector( tofSD );
    s_Tof->SetSensitiveDetector( tofSD );


    ///////////////////////////////////////////// 
    ////////////////////////////////Internal TOF
    //////////////////////////////////////////// 
    
    //ITOF
    G4int iditofr=geomMan.GetDetectorId("ITOFR");
    G4ThreeVector gPosITofR=geomMan.GetGlobalPosition(iditofr);
    G4ThreeVector OfsLITofR( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMITOFR;  
    RMITOFR.rotateX( 90.*degree );
    RMITOFR.rotateY( ITof_RotAngleR );

    s_ITofWall *s_ITofR =
      new s_ITofWall( "ITOFR", DetLV_, RMITOFR, gPosITofR, OfsLITofR, 
    		      mList_->Scin, mList_->Air);

    G4int iditofl=geomMan.GetDetectorId("ITOFL");
    G4ThreeVector gPosITofL=geomMan.GetGlobalPosition(iditofl);
    G4ThreeVector OfsLITofL( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMITOFL;  
    RMITOFL.rotateX( 90.*degree );
    RMITOFL.rotateY( ITof_RotAngleL );

    s_ITofWall *s_ITofL =
      new s_ITofWall( "ITOFL", DetLV_, RMITOFL, gPosITofL, OfsLITofL, 
    		      mList_->Scin, mList_->Air);
    
    //ITOF side wall
    G4int iditofsr=geomMan.GetDetectorId("ITOFSR");
    G4ThreeVector gPosITofSR=geomMan.GetGlobalPosition(iditofsr);
    G4ThreeVector OfsLITofSR( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMITOFSR;  
    RMITOFSR.rotateX( 90.*degree );
    RMITOFSR.rotateY( ITofS_RotAngleR );

    s_ITofSWall *s_ITofSR =
      new s_ITofSWall( "ITOFSR", DetLV_, RMITOFSR, gPosITofSR, OfsLITofSR, 
    		       mList_->Scin, mList_->Air);

    G4int iditofsl=geomMan.GetDetectorId("ITOFSL");
    G4ThreeVector gPosITofSL=geomMan.GetGlobalPosition(iditofsl);
    G4ThreeVector OfsLITofSL( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMITOFSL;  
    RMITOFSL.rotateX( 90.*degree );
    RMITOFSL.rotateY( ITofS_RotAngleL );

    s_ITofSWall *s_ITofSL =
      new s_ITofSWall( "ITOFSL", DetLV_, RMITOFSL, gPosITofSL, OfsLITofSL, 
    		       mList_->Scin, mList_->Air);

    G4VisAttributes *ITofScinAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 1.0, 0.0, 1.0 ) );//Magenta

    s_ITofR->SetVisAttributes( ITofScinAttrib, &G4VisAttributes::Invisible );
    s_ITofL->SetVisAttributes( ITofScinAttrib, &G4VisAttributes::Invisible );
    s_ITofSR->SetVisAttributes( ITofScinAttrib, &G4VisAttributes::Invisible );
    s_ITofSL->SetVisAttributes( ITofScinAttrib, &G4VisAttributes::Invisible );
  
    //////// Sensitive Detectors ////////
    ITofSD *itofSD = new ITofSD( "/spec/itof" );
    SDMan->AddNewDetector( itofSD );
    s_ITofR->SetSensitiveDetector( itofSD );
    s_ITofL->SetSensitiveDetector( itofSD );
    s_ITofSR->SetSensitiveDetector( itofSD );
    s_ITofSL->SetSensitiveDetector( itofSD );


    ///////////////////////////////////////////// 
    ////////////////////////////////PAD
    //////////////////////////////////////////// 
    
    //PAD
    G4int idpadu=geomMan.GetDetectorId("PADU");
    G4ThreeVector gPosPADU=geomMan.GetGlobalPosition(idpadu);
    G4ThreeVector OfsLPADU( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMPADU;  
    RMPADU.rotateX( 0.*degree );
    RMPADU.rotateY( 0.*degree );
    RMPADU.rotateZ( 180.*degree );

    // s_PAD *s_PADU =
    //   new s_PAD( "PADU", DetLV_, RMPADU, gPosPADU, OfsLPADU, 
    // 		 mList_->Scin, mList_->Air);

    G4int idpadd=geomMan.GetDetectorId("PADD");
    G4ThreeVector gPosPADD=geomMan.GetGlobalPosition(idpadd);
    G4ThreeVector OfsLPADD( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMPADD;  
    RMPADD.rotateX(  0.*degree );
    RMPADD.rotateY(  0.*degree );
    RMPADD.rotateZ( 180.*degree );
    
    s_PAD *s_PADD =
      new s_PAD( "PADD", DetLV_, RMPADD, gPosPADD, OfsLPADD, 
    		 mList_->Scin, mList_->Air);

    G4VisAttributes *PADScinAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 1.0, 0.0, 1.0 ) );//Magenta
    
    // s_PADU->SetVisAttributes( PADScinAttrib, &G4VisAttributes::Invisible );
    s_PADD->SetVisAttributes( PADScinAttrib, &G4VisAttributes::Invisible );
  
    //////// Sensitive Detectors ////////
    PADSD *padSD = new PADSD( "/spec/pad" );
    SDMan->AddNewDetector( padSD );
    // s_PADU->SetSensitiveDetector( padSD );
    s_PADD->SetSensitiveDetector( padSD );

    
    ///////////////////////////////////////// 
    ////////////////////////////////RICH
    ///////////////////////////////////////// 
    
    G4int idrich=geomMan.GetDetectorId("RICH");
    G4ThreeVector gPosRICH=geomMan.GetGlobalPosition(idrich);
    G4ThreeVector OfsLRICH( 0.0, 0.0, 0.0 );
    
    s_RICH *s_Rich =
      new s_RICH( "RICH", DetLV_, RM, gPosRICH, OfsLRICH, 
		  mList_->Aerogel, 
		  mList_->Air, //Make C4F10
		  mList_->Al, 
		  mList_->Air );

    G4VisAttributes *RICHAcAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 1.0, 0.5, 0.0 ) );//Orange
    G4VisAttributes *RICHBoxAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.5, 0.5, 1.0 ) );//BluePurple
    // G4VisAttributes *RICHGasAttrib =
    //   new G4VisAttributes( G4Colour::G4Colour( 0.0, 1.0, 0.0 ) );//Green
    
    s_Rich->SetVisAttributes( RICHAcAttrib, 
    			      RICHBoxAttrib, 
    			      &G4VisAttributes::Invisible,
    			      &G4VisAttributes::Invisible );

    // s_Rich->SetVisAttributes( RICHAcAttrib, 
    // 			      RICHBoxAttrib, 
    // 			      RICHGasAttrib, 
    // 			      0 );
  
    //////// Sensitive Detectors ////////
    RICHSD *richSD = new RICHSD( "/spec/rich" );
    SDMan->AddNewDetector( richSD );
    s_Rich->SetSensitiveDetector( richSD );


    /////////////////////////////////////////// 
    ////////////////////////////////PID1&&PID2
    /////////////////////////////////////////// 
    
    G4int idpid1r=geomMan.GetDetectorId("PID1R");
    G4ThreeVector gPosPID1R=geomMan.GetGlobalPosition(idpid1r);
    G4ThreeVector OfsLPID1R( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMPID1R;  
    RMPID1R.rotateX( 90.*degree );
    RMPID1R.rotateY( 180.*degree );
    
    s_PID1 *s_Pid1R =
      new s_PID1( "PID1R", DetLV_, RMPID1R, gPosPID1R, OfsLPID1R, 
		  mList_->Aerogel, 
		  mList_->Air,
		  mList_->Al, 
		  mList_->Air );

    G4int idpid1l=geomMan.GetDetectorId("PID1L");
    G4ThreeVector gPosPID1L=geomMan.GetGlobalPosition(idpid1l);
    G4ThreeVector OfsLPID1L( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMPID1L;  
    RMPID1L.rotateX( 90.*degree );
    RMPID1L.rotateY( 180.*degree );
    
    s_PID1 *s_Pid1L =
      new s_PID1( "PID1L", DetLV_, RMPID1L, gPosPID1L, OfsLPID1L, 
		  mList_->Aerogel, 
		  mList_->Air,
		  mList_->Al, 
		  mList_->Air );

    G4int idpid2r=geomMan.GetDetectorId("PID2R");
    G4ThreeVector gPosPID2R=geomMan.GetGlobalPosition(idpid2r);
    G4ThreeVector OfsLPID2R( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMPID2R;  
    RMPID2R.rotateX( 90.*degree );
    RMPID2R.rotateY( PID2_RotAngleR );
    
    s_PID2 *s_Pid2R =
      new s_PID2( "PID2R", DetLV_, RMPID2R, gPosPID2R, OfsLPID2R, 
    		  mList_->Aerogel, 
    		  mList_->Air,
    		  mList_->Al, 
    		  mList_->Air );

    G4int idpid2l=geomMan.GetDetectorId("PID2L");
    G4ThreeVector gPosPID2L=geomMan.GetGlobalPosition(idpid2l);
    G4ThreeVector OfsLPID2L( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMPID2L;  
    RMPID2L.rotateX( 90.*degree );
    RMPID2L.rotateY( PID2_RotAngleL );
    
    s_PID2 *s_Pid2L =
      new s_PID2( "PID2L", DetLV_, RMPID2L, gPosPID2L, OfsLPID2L, 
    		  mList_->Aerogel, 
    		  mList_->Air,
    		  mList_->Al, 
    		  mList_->Air );

    G4VisAttributes *PID1SegAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 1.0, 0.5, 0.0 ) );//Orange
    G4VisAttributes *PID1BoxAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.0, 1.0, 0.0 ) );//Green
    
    s_Pid1R->SetVisAttributes( PID1SegAttrib, 
			       PID1BoxAttrib, 
			       &G4VisAttributes::Invisible,
			       &G4VisAttributes::Invisible );
    s_Pid1L->SetVisAttributes( PID1SegAttrib, 
			       PID1BoxAttrib, 
			       &G4VisAttributes::Invisible,
			       &G4VisAttributes::Invisible );
    s_Pid2R->SetVisAttributes( PID1SegAttrib, 
    			       PID1BoxAttrib, 
    			       &G4VisAttributes::Invisible,
    			       &G4VisAttributes::Invisible );
    s_Pid2L->SetVisAttributes( PID1SegAttrib, 
    			       PID1BoxAttrib, 
    			       &G4VisAttributes::Invisible,
    			       &G4VisAttributes::Invisible );

    //////// Sensitive Detectors ////////
    PID1SD *pid1SD = new PID1SD( "/spec/pid1" );
    PID2SD *pid2SD = new PID2SD( "/spec/pid2" );
    SDMan->AddNewDetector( pid1SD );
    SDMan->AddNewDetector( pid2SD );
    s_Pid1R->SetSensitiveDetector( pid1SD );
    s_Pid1L->SetSensitiveDetector( pid1SD );
    s_Pid2R->SetSensitiveDetector( pid2SD );
    s_Pid2L->SetSensitiveDetector( pid2SD );

    ///////////////////////////////////////////// 
    ////////////////////////////////Muon Filter
    //////////////////////////////////////////// 

    G4int idmfu=geomMan.GetDetectorId("MFU");
    G4int idmfd=geomMan.GetDetectorId("MFD");
    G4ThreeVector gPosMFU=geomMan.GetGlobalPosition(idmfu);
    G4ThreeVector gPosMFD=geomMan.GetGlobalPosition(idmfd);
    G4ThreeVector OfsLMF( 0.0, 0.0, 0.0 );
    G4RotationMatrix RMMF;  
    RMMF.rotateX( 90.*degree );
    RMMF.rotateY( 180.*degree );
    
    s_MuonFilter *s_MF =
      new s_MuonFilter( "MuonFilter", DetLV_, RMMF, 
			gPosMFU, gPosMFD,
			OfsLPID1R, 
			mList_->Scin, 
			mList_->Fe,
			mList_->Air );

    G4VisAttributes *MFScinAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 1.0, 0.0, 1.0 ) );//Magenta
    G4VisAttributes *MFIronAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 1.0, 1.0, 1.0 ) );//White
    
    s_MF->SetVisAttributes( MFScinAttrib,
			    MFIronAttrib, 
			    &G4VisAttributes::Invisible );

    //////// Sensitive Detectors ////////
    MuonFilterSD *mfSD = new MuonFilterSD( "/spec/mf" );
    SDMan->AddNewDetector( mfSD );
    s_MF->SetSensitiveDetector( mfSD );

    ///////////////////////////////////////////////////
  }//s_DetFlag end

  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////
  //Real detectors
  if( r_DetFlag ){

    ///////////////////////////////////////// 
    ////////////////////////////////T0
    ///////////////////////////////////////// 
    G4int idt0=geomMan.GetDetectorId("T0");
    G4ThreeVector gPosT0=geomMan.GetGlobalPosition(idt0);
    G4ThreeVector OfsLT0( 0.0, 0.0, 0.0 );

    r_T0Wall *r_T0 =
      new r_T0Wall( "T0", DetLV_, RM, gPosT0, OfsLT0, 
  		    mList_->Scin, mList_->Scin,
  		    mList_->Air, mList_->Air );
    
    G4VisAttributes *T0ScinAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 1.0, 0.0, 1.0 ) );
    G4VisAttributes *T0MPPCAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.0, 0.0, 0.0 ) );
    G4VisAttributes *T0PCBAttrib =
      new G4VisAttributes( G4Colour::G4Colour( 0.0, 1.0, 0.0 ) );

    
    r_T0->SetVisAttributes( T0ScinAttrib, T0MPPCAttrib, T0PCBAttrib,
  			    &G4VisAttributes::Invisible );
  

    //////// Sensitive Detectors ////////
    T0SD *t0SD = new T0SD( "/spec/t0" );
    SDMan->AddNewDetector( t0SD );
    r_T0->SetSensitiveDetector( t0SD );
  }


  //////////////////////////////////////////
  ////////////////////////////////VD
  ///////////////////////////////////////// 

  G4VisAttributes *VDAttrib =
    new G4VisAttributes( G4Colour::G4Colour( 2.0, 2.0, 2.0 ) );//White
  
  //VD1 for Beam information
  G4Box *solidVD1 = 
    new G4Box( "VD1", 300./2., 1.0/2., 200./2. );
  G4LogicalVolume *logVD1 =
    new G4LogicalVolume( solidVD1, mList_->Air, "VD1", 0, 0, 0 );

  G4int idvd1=geomMan.GetDetectorId("SpecVp1");
  G4ThreeVector gPosVD1=geomMan.GetGlobalPosition(idvd1);
  G4ThreeVector OfsLVD1( 0.0, 0.0, 0.0 );

  G4VPhysicalVolume *physVD1 =
    new G4PVPlacement( G4Transform3D( RM,
                                      G4ThreeVector( gPosVD1.x()*mm,
                                                     gPosVD1.y()*mm,
                                                     gPosVD1.z()*mm ) ),
  		       logVD1, "VD1", DetLV_, false, 1 );

  //VD2 for IT2 downstream
  G4Box *solidVD2 = 
    new G4Box( "VD2", 2000./2., 1./2., 950./2. );
  G4LogicalVolume *logVD2 =
    new G4LogicalVolume( solidVD2, mList_->Air, "VD2", 0, 0, 0 );

  G4int idvd2=geomMan.GetDetectorId("SpecVp2");
  G4ThreeVector gPosVD2=geomMan.GetGlobalPosition(idvd2);
  G4ThreeVector OfsLVD2( 0.0, 0.0, 0.0 );

  G4VPhysicalVolume *physVD2 =
    new G4PVPlacement( G4Transform3D( RM,
                                      G4ThreeVector( gPosVD2.x()*mm,
                                                     gPosVD2.y()*mm,
                                                     gPosVD2.z()*mm ) ),
		       logVD2, "VD2", DetLV_, false, 2 );

  //VD3 for RICH information
  G4Box *solidVD3 = 
    new G4Box( "VD3", 5000./2., 1./2., 4000./2. );
  G4LogicalVolume *logVD3 =
    new G4LogicalVolume( solidVD3, mList_->Air, "VD3", 0, 0, 0 );

  G4int idvd3=geomMan.GetDetectorId("SpecVp3");
  G4ThreeVector gPosVD3=geomMan.GetGlobalPosition(idvd3);
  G4ThreeVector OfsLVD3( 0.0, 0.0, 0.0 );

  G4VPhysicalVolume *physVD3 =
    new G4PVPlacement( G4Transform3D( RM,
                                      G4ThreeVector( gPosVD3.x()*mm,
                                                     gPosVD3.y()*mm,
                                                     gPosVD3.z()*mm ) ),
		       logVD3, "VD3", DetLV_, false, 3 );

  // logVD1->SetVisAttributes(VDAttrib);
  // logVD2->SetVisAttributes(VDAttrib);
  // logVD3->SetVisAttributes(VDAttrib);

  logVD1->SetVisAttributes(&G4VisAttributes::Invisible);
  logVD2->SetVisAttributes(&G4VisAttributes::Invisible);
  logVD3->SetVisAttributes(&G4VisAttributes::Invisible);
  
  //////// Sensitive Detectors ////////
  VDSD *vdSD = new VDSD( "/Sks/VD" );
  SDMan->AddNewDetector( vdSD );
  logVD1->SetSensitiveDetector( vdSD );
  logVD2->SetSensitiveDetector( vdSD );
  logVD3->SetSensitiveDetector( vdSD );

  return;
}

