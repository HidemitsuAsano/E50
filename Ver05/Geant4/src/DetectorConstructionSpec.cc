/*
  DetectorConstructionSpec.hh

  FM version 2.0

  2016/1  K.Shirotori
*/

#include "DetectorConstructionSpec.hh"

#include "MaterialList.hh"
#include "SpecField.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "LiqTarget.hh"

#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"

#include "G4Box.hh" 
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"

const G4double DefTargSizX = 100.0*mm;
const G4double DefTargSizY = 100.0*mm;
const G4double DefTargSizZ =   1.0*mm;

DetectorConstructionSpec::DetectorConstructionSpec()
  : targSizX(DefTargSizX), targSizY(DefTargSizY), targSizZ(DefTargSizZ)
{
}

DetectorConstructionSpec::~DetectorConstructionSpec()
{
}


G4ThreeVector DetectorConstructionSpec::TargetPosition( void ) const
{
  return targPos;
}


G4double DetectorConstructionSpec::TargetAngle( void ) const
{
  return targAng2;
}


G4double DetectorConstructionSpec::TargetSizeX( void ) const
{ 
  return targSizX;
}

G4double DetectorConstructionSpec::TargetSizeY( void ) const
{ 
  return targSizY;
}

G4double DetectorConstructionSpec::TargetSizeZ( void ) const
{ 
  return targSizZ;
}

void DetectorConstructionSpec::SetRealMaterials( void )
{
  matWorld = mList_->HeGas;

  matSpecGap = mList_->HeGas;

  //BT1&2
  matBT1Base  = mList_->HeGas;
  matBT1Layer = mList_->Si;
  matBT2Base  = mList_->HeGas;
  matBT2Layer = mList_->Si;

  //T0
  matT0Scin = mList_->Scin;
  matT0LG   = mList_->Scin;
  matT0PMT  = mList_->Fe;  

  //ST1
  matST1Base  = mList_->HeGas;
  matST1Layer = mList_->Si;

  //ST2
  matST2Base  = mList_->HeGas;
  matST2Layer = mList_->Scin;

  //ST3&4
  matST3Gas = mList_->Ar80IsoButane20Gas;
  matST3Frame = mList_->Al;
  matST4Gas = mList_->Ar80IsoButane20Gas;
  matST4Frame = mList_->Al;

  //TOF
  matTofScin = mList_->Scin;
  matTofLG = mList_->Scin;
  matTofPMT = mList_->Fe;  
  matTofPMTBox = mList_->Fe;
  matTofFrame = mList_->Al;

  //VD (Virtual Detector)
  matVD = mList_->HeGas;

  ConfMan *confMan = ConfMan::GetConfManager();
  G4Material *mat=0;
  G4int idWorld=confMan->GetMatIdWorld();
  if((mat=mList_->chooseMaterial(idWorld))!=0) 
    matWorld=mat;
  matT0Box  = matWorld;
  matBT1Box = matWorld;
  matBT2Box = matWorld;
  matST1Box = matWorld;
  matST2Box = matWorld;
  matST3Box = matWorld;
  matST4Box = matWorld;
  matTofBox = matWorld;

  G4int idDcGas2=confMan->GetMatIdDcGas2();
  if((mat=mList_->chooseMaterial(idDcGas2))!=0) 
    matST3Gas=matST4Gas=mat;
  G4int idMagGap=confMan->GetMatIdMagGap();
  if((mat=mList_->chooseMaterial(idMagGap))!=0) 
    matSpecGap=mat;
  G4int idPMT=confMan->GetMatIdPMT();
  if((mat=mList_->chooseMaterial(idPMT))!=0) 
    matTofPMT=mat;
  G4int idPMTBox=confMan->GetMatIdPMTBox();
  if((mat=mList_->chooseMaterial(idPMTBox))!=0) 
    matTofPMTBox=mat;
  G4int idTgt=confMan->GetMatIdTarget();
  if((mat=mList_->chooseMaterial(idTgt))!=0) 
    matTarget=mat;
  else
    matTarget = matWorld;
  
  PrintRealMaterialName();
}


void DetectorConstructionSpec::PrintRealMaterialName( void ) const
{
  // G4cout << "========== Real Material Name ============" << G4endl;
  // G4cout << "World: [" << matWorld->GetName() << "] " 
  // 	 << "MagGap: [" << matSpecGap->GetName() << "] " << G4endl;
  // G4cout << "BT1Base: [" << matBT1Base->GetName() << "] " 
  // 	 << "BT1Layer: ["  << matBT1Layer->GetName() << "] " 
  // 	 << "Box(Region): [" << matBT1Box->GetName() << "] " << G4endl;
  // G4cout << "BT2Base: [" << matBT2Base->GetName() << "] " 
  // 	 << "BT2Layer: ["  << matBT2Layer->GetName() << "] " 
  // 	 << "Box(Region): [" << matBT2Box->GetName() << "] " << G4endl;
  // G4cout << "T0 Scin.: [" << matT0Scin->GetName() << "] " 
  // 	 << "LG: [" << matT0LG->GetName() << "] "
  // 	 << "PMT: [" << matT0PMT->GetName() << "] "
  // 	 << "Box(Region): [" << matT0Box->GetName() << "] " << G4endl;
  // G4cout << "ST1Base: [" << matST1Base->GetName() << "] " 
  // 	 << "ST1Layer: ["  << matST1Layer->GetName() << "] " 
  // 	 << "Box(Region): [" << matST1Box->GetName() << "] " << G4endl;
  // G4cout << "ST2Base: [" << matST2Base->GetName() << "] " 
  // 	 << "ST2Layer: ["  << matST2Layer->GetName() << "] " 
  // 	 << "Box(Region): [" << matST2Box->GetName() << "] " << G4endl;
  // G4cout << "ST3 Gas: [" << matST3Gas->GetName() << "] " 
  // 	 << "Frame: ["  << matST3Frame->GetName() << "] "
  // 	 << "Box(Region): [" << matST3Box->GetName() << "] " << G4endl;
  // G4cout << "ST4 Gas: [" << matST4Gas->GetName() << "] " 
  // 	 << "Frame: ["  << matST4Frame->GetName() << "] "
  // 	 << "Box(Region): [" << matST4Box->GetName() << "] " << G4endl;
  // G4cout << "Tof Scin.: [" << matTofScin->GetName() << "] " 
  // 	 << "LG: [" << matTofLG->GetName() << "] "
  // 	 << "PMT: [" << matTofPMT->GetName() << "] "
  // 	 << "PMTBox: [" << matTofPMTBox->GetName() << "] "
  // 	 << "Frame: [" << matTofFrame->GetName() << "] "
  // 	 << "Box(Region): [" << matTofBox->GetName() << "] " << G4endl;
  // G4cout << "VD: [" << matVD->GetName() << "] " << G4endl;
  // G4cout << "Target: [" << matTarget->GetName() << "] " << G4endl;
  // G4cout << "==========================================" << G4endl;

}


G4VPhysicalVolume *
DetectorConstructionSpec::ConstructPayload( void )
{
  SetRealMaterials();

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int lnumTarg = geomMan.GetDetectorId("Target");
  targPos = geomMan.GetGlobalPosition( lnumTarg )*mm;
  targAng2 = geomMan.GetRotAngle2( lnumTarg )*degree;

  ConfMan *confMan = ConfMan::GetConfManager();

  G4double size=0.0;
  if( (size=confMan->GetTargetSizeX())>0.0 ) targSizX=size;
  if( (size=confMan->GetTargetSizeY())>0.0 ) targSizY=size;
  if( (size=confMan->GetTargetSizeZ())>0.0 ) targSizZ=size;

  ///////////////////// World
  G4Box *solidWorld =
    new G4Box( "WorldBox", 1500.*cm/2., 1500.*cm/2., 1500.*cm/2. );
  G4LogicalVolume *logWorld =
    new G4LogicalVolume( solidWorld, matWorld, "World", 0, 0, 0, false );
  logWorld->SetVisAttributes( G4VisAttributes::Invisible );
  G4VPhysicalVolume *physWorld =
    new G4PVPlacement( 0, G4ThreeVector(), "World", logWorld, 0, false, 0 );

  MakeSpecMagnet( physWorld, matSpecGap );
  if(confMan->ExistLiqTarget()){
    MakeLiqTarget( physWorld );
  }
  else if(confMan->ExistTarget()){
    MakeTarget( physWorld );
  }
  MakeTrackers( physWorld );
  MakeTofCounters( physWorld );
  MakeVDetector( physWorld );

  return physWorld;
}


void DetectorConstructionSpec::
MakeTarget( G4VPhysicalVolume *pMother )
{
  // //G4UserLimits *volStepLimit = new G4UserLimits(0.01*mm);

  // G4VisAttributes *TgtAttrib =
  //  new G4VisAttributes( G4Colour::G4Colour( 0.0, 1.0, 1.0 ) );

  // //LH2 cell
  // G4Tubs *solidTGT =
  //   new G4Tubs( "Target",  0.0*mm, 100.0/2.*mm, 570./2.*mm, 0.0*degree, 360.0*degree );
  // G4LogicalVolume *logTGT =
  //   new G4LogicalVolume( solidTGT, matTarget, "Target", 0, 0, 0 );
  // G4RotationMatrix RM; RM.rotateX( 0.*degree );RM.rotateZ(targAng2);
  // // G4VPhysicalVolume *physTGT =
  // //   new G4PVPlacement( G4Transform3D( RM, targPos ), "Target",
  // // 		       logTGT, pMother, false, 0 );

  // logTGT->SetVisAttributes(TgtAttrib);
  //logTGT->SetUserLimits(volStepLimit);

  //////// Sensitive Detectors ////////
  // G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  // TargetSD *tgtSD = new TargetSD( "/spec/target" );
  // SDMan->AddNewDetector( tgtSD );
  // logTGT->SetSensitiveDetector( tgtSD );
}


void DetectorConstructionSpec::
MakeLiqTarget( G4VPhysicalVolume *pMother )
{
  // G4RotationMatrix RMTarg; RMTarg.rotateZ(targAng2);
  // G4ThreeVector OfsTarg( 0.0, 0.0, 0.0 );
  
  // LiqTarget *Targ =
  //   new LiqTarget( "LiqTarget", pMother, RMTarg, targPos, OfsTarg,
  // 		   matTarget, mList_->PET, mList_->Al, mList_->PET,
  // 		   mList_->Vacuum, matWorld );
  // G4VisAttributes *TgtAttrib =
  //  new G4VisAttributes( G4Colour::G4Colour( 0.0, 1.0, 1.0 ) );
  // G4VisAttributes *TgtCellAttrib =
  //  new G4VisAttributes( G4Colour::G4Colour( 0.0, 0.5, 0.5 ) );
  // G4VisAttributes *TgtVesselAttrib =
  //  new G4VisAttributes( G4Colour::G4Colour( 0.5, 0.5, 0.5 ) );

  // Targ->SetVisAttributes( TgtAttrib, TgtCellAttrib, TgtVesselAttrib,
  // 			  &G4VisAttributes::Invisible,
  // 			  TgtCellAttrib, &G4VisAttributes::Invisible );
  // //  Targ->SetVisAttributes( TgtAttrib, TgtCellAttrib, 
  // //			  &G4VisAttributes::Invisible,
  // //			  &G4VisAttributes::Invisible,
  // //			  TgtCellAttrib, &G4VisAttributes::Invisible );
  // //  Targ->SetVisAttributes( TgtAttrib, TgtCellAttrib, TgtVesselAttrib,
  // //			  0,TgtCellAttrib, &G4VisAttributes::Invisible );

  // targSizX=targSizY=Targ->PhiSize();
  // targSizZ=Targ->Length();
}


void DetectorConstructionSpec::
MakeTrackers( G4VPhysicalVolume *pMother )
{
  // const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  // G4VisAttributes *LayerAttrib =
  //   new G4VisAttributes( G4Colour::G4Colour( 0.0, 1.0, 1.0 ) );
  // G4VisAttributes *FrameAttrib =
  //   new G4VisAttributes( G4Colour::G4Colour( 1.0, 1.0, 1.0 ) );

  // ///////////////////// BT1(v,x,u,v,x,u)
  // G4int id101=geomMan.GetDetectorId("BT1-v-1");
  // G4int id102=geomMan.GetDetectorId("BT1-x-1");
  // G4int id103=geomMan.GetDetectorId("BT1-u-1");
  // G4int id104=geomMan.GetDetectorId("BT1-v-2");
  // G4int id105=geomMan.GetDetectorId("BT1-x-2");
  // G4int id106=geomMan.GetDetectorId("BT1-u-2");

  // G4ThreeVector gPos101=geomMan.GetGlobalPosition(id101);
  // G4ThreeVector gPos102=geomMan.GetGlobalPosition(id102);
  // G4ThreeVector gPos103=geomMan.GetGlobalPosition(id103);
  // G4ThreeVector gPos104=geomMan.GetGlobalPosition(id104);
  // G4ThreeVector gPos105=geomMan.GetGlobalPosition(id105);
  // G4ThreeVector gPos106=geomMan.GetGlobalPosition(id106);
  // G4double RotationAngleBT1=geomMan.GetRotAngle2(id103);
  // G4RotationMatrix RMBT1; RMBT1.rotateZ(RotationAngleBT1*degree);
  // G4ThreeVector OfsLBT1( 0.0, 0.0, 0.0 );

  // SSD6layers1 *BT1 =
  //   new SSD6layers1( "BT1", pMother, RMBT1, 
  // 		     gPos101, gPos102, gPos103, gPos104, gPos105, gPos106,
  // 		     OfsLBT1, 
  // 		     id101, id102, id103, id104, id105, id106, 
  // 		     matBT1Layer, matBT1Base, matBT1Box );
  
  // BT1->SetVisAttributes( LayerAttrib, FrameAttrib, &G4VisAttributes::Invisible );


  // ///////////////////// BT2(v,x,u,v,x,u)
  // G4int id107=geomMan.GetDetectorId("BT2-v-1");
  // G4int id108=geomMan.GetDetectorId("BT2-x-1");
  // G4int id109=geomMan.GetDetectorId("BT2-u-1");
  // G4int id110=geomMan.GetDetectorId("BT2-v-2");
  // G4int id111=geomMan.GetDetectorId("BT2-x-2");
  // G4int id112=geomMan.GetDetectorId("BT2-u-2");

  // G4ThreeVector gPos107=geomMan.GetGlobalPosition(id107);
  // G4ThreeVector gPos108=geomMan.GetGlobalPosition(id108);
  // G4ThreeVector gPos109=geomMan.GetGlobalPosition(id109);
  // G4ThreeVector gPos110=geomMan.GetGlobalPosition(id110);
  // G4ThreeVector gPos111=geomMan.GetGlobalPosition(id111);
  // G4ThreeVector gPos112=geomMan.GetGlobalPosition(id112);
  // G4double RotationAngleBT2=geomMan.GetRotAngle2(id109);
  // G4RotationMatrix RMBT2; RMBT2.rotateZ(RotationAngleBT2*degree);
  // G4ThreeVector OfsLBT2( 0.0, 0.0, 0.0 );

  // SSD6layers1 *BT2 =
  //   new SSD6layers1( "BT2", pMother, RMBT2, 
  // 		     gPos107, gPos108, gPos109, gPos110, gPos111, gPos112,
  // 		     OfsLBT2, 
  // 		     id107, id108, id109, id110, id111, id112, 
  // 		     matBT2Layer, matBT2Base, matBT2Box );
  
  // BT2->SetVisAttributes( LayerAttrib, FrameAttrib, &G4VisAttributes::Invisible);


  // ///////////////////// ST1(v,x,u)
  // G4int id1=geomMan.GetDetectorId("ST1-v-1");
  // G4int id2=geomMan.GetDetectorId("ST1-x-1");
  // G4int id3=geomMan.GetDetectorId("ST1-u-1");
  // G4int id4=geomMan.GetDetectorId("ST1-v-2");
  // G4int id5=geomMan.GetDetectorId("ST1-x-2");
  // G4int id6=geomMan.GetDetectorId("ST1-u-2");

  // G4ThreeVector gPos1=geomMan.GetGlobalPosition(id1);
  // G4ThreeVector gPos2=geomMan.GetGlobalPosition(id2);
  // G4ThreeVector gPos3=geomMan.GetGlobalPosition(id3);
  // G4ThreeVector gPos4=geomMan.GetGlobalPosition(id4);
  // G4ThreeVector gPos5=geomMan.GetGlobalPosition(id5);
  // G4ThreeVector gPos6=geomMan.GetGlobalPosition(id6);
  // G4double RotationAngleST1=geomMan.GetRotAngle2(id3);
  // G4RotationMatrix RMST1; RMST1.rotateZ(RotationAngleST1*degree);
  // G4ThreeVector OfsLST1( 0.0, 0.0, 0.0 );

  // SSD6layers2 *ST1 =
  //   new SSD6layers2( "ST1", pMother, RMST1, 
  // 		     gPos1, gPos2, gPos3,
  // 		     gPos4, gPos5, gPos6,
  // 		     OfsLST1, 
  // 		     id1, id2, id3,
  // 		     id4, id5, id6,
  // 		     matST1Layer, matST1Base, matST1Box );
  
  //  ST1->SetVisAttributes( LayerAttrib, FrameAttrib, &G4VisAttributes::Invisible );


  // ///////////////////// ST2(v,x,u,v,x,u,v,x,u)
  // G4int id7=geomMan.GetDetectorId("ST2-v-1");
  // G4int id8=geomMan.GetDetectorId("ST2-x-1");
  // G4int id9=geomMan.GetDetectorId("ST2-u-1");
  // G4int id10=geomMan.GetDetectorId("ST2-v-2");
  // G4int id11=geomMan.GetDetectorId("ST2-x-2");
  // G4int id12=geomMan.GetDetectorId("ST2-u-2");
  // G4int id13=geomMan.GetDetectorId("ST2-v-3");
  // G4int id14=geomMan.GetDetectorId("ST2-x-3");
  // G4int id15=geomMan.GetDetectorId("ST2-u-3");

  // G4ThreeVector gPos7=geomMan.GetGlobalPosition(id7);
  // G4ThreeVector gPos8=geomMan.GetGlobalPosition(id8);
  // G4ThreeVector gPos9=geomMan.GetGlobalPosition(id9);
  // G4ThreeVector gPos10=geomMan.GetGlobalPosition(id10);
  // G4ThreeVector gPos11=geomMan.GetGlobalPosition(id11);
  // G4ThreeVector gPos12=geomMan.GetGlobalPosition(id12);
  // G4ThreeVector gPos13=geomMan.GetGlobalPosition(id13);
  // G4ThreeVector gPos14=geomMan.GetGlobalPosition(id14);
  // G4ThreeVector gPos15=geomMan.GetGlobalPosition(id15);
  // G4double RotationAngleST2=geomMan.GetRotAngle2(id11);
  // G4RotationMatrix RMST2; RMST2.rotateZ(RotationAngleST2*degree);
  // G4ThreeVector OfsLST2( 0.0, 0.0, 0.0 );

  // Fiber9layers1 *ST2 =
  //   new Fiber9layers1( "ST2", pMother, RMST2, 
  // 		       gPos7, gPos8, gPos9, 
  // 		       gPos10, gPos11, gPos12, 
  // 		       gPos13, gPos14, gPos15,
  // 		       OfsLST2, 
  // 		       id7, id8, id9, 
  // 		       id10, id11, id12, 
  // 		       id13, id14, id15, 
  // 		       matST2Layer, matST2Base, matST2Box );
  
  //  ST2->SetVisAttributes( LayerAttrib, FrameAttrib, &G4VisAttributes::Invisible );


  // ///////////////////// ST3(v,x,u,v,x,u)
  // G4int id51 =geomMan.GetDetectorId("ST3-v-1");
  // G4int id52 =geomMan.GetDetectorId("ST3-x-1");
  // G4int id53 =geomMan.GetDetectorId("ST3-u-1");
  // G4int id54 =geomMan.GetDetectorId("ST3-v-2");
  // G4int id55 =geomMan.GetDetectorId("ST3-x-2");
  // G4int id56 =geomMan.GetDetectorId("ST3-u-2");
  // G4ThreeVector gPos51 =geomMan.GetGlobalPosition(id51);
  // G4ThreeVector gPos52 =geomMan.GetGlobalPosition(id52);
  // G4ThreeVector gPos53 =geomMan.GetGlobalPosition(id53);
  // G4ThreeVector gPos54 =geomMan.GetGlobalPosition(id54);
  // G4ThreeVector gPos55 =geomMan.GetGlobalPosition(id55);
  // G4ThreeVector gPos56 =geomMan.GetGlobalPosition(id56);
  // G4double RotationAngleST3=geomMan.GetRotAngle2(id53);
  // G4RotationMatrix RMST3; RMST3.rotateZ(RotationAngleST3*degree);
  // G4ThreeVector OfsLST3(  0.0, 0.0, 0.0 );

  // STChamber1 *ST3 =
  //   new STChamber1( "ST3", pMother, RMST3, 
  // 		    gPos51, gPos52, gPos53, gPos54, gPos55, gPos56, 
  // 		    OfsLST3, 
  // 		    id51, id52, id53, id54, id55, id56, 
  // 		    matST3Gas, matST3Frame, matST3Box );
  // ST3->SetVisAttributes( &G4VisAttributes::Invisible,
  // 			 FrameAttrib, LayerAttrib,
  // 			 &G4VisAttributes::Invisible );


  // ///////////////////// ST4(v,x,u,v,x,u)
  // G4int id57 =geomMan.GetDetectorId("ST4-v-1");
  // G4int id58 =geomMan.GetDetectorId("ST4-x-1");
  // G4int id59 =geomMan.GetDetectorId("ST4-u-1");
  // G4int id60 =geomMan.GetDetectorId("ST4-v-2");
  // G4int id61 =geomMan.GetDetectorId("ST4-x-2");
  // G4int id62 =geomMan.GetDetectorId("ST4-u-2");
  // G4ThreeVector gPos57 =geomMan.GetGlobalPosition(id57);
  // G4ThreeVector gPos58 =geomMan.GetGlobalPosition(id58);
  // G4ThreeVector gPos59 =geomMan.GetGlobalPosition(id59);
  // G4ThreeVector gPos60 =geomMan.GetGlobalPosition(id60);
  // G4ThreeVector gPos61 =geomMan.GetGlobalPosition(id61);
  // G4ThreeVector gPos62 =geomMan.GetGlobalPosition(id62);
  // G4double RotationAngleST4=geomMan.GetRotAngle2(id59);
  // G4RotationMatrix RMST4; RMST4.rotateZ(RotationAngleST4*degree);
  // G4ThreeVector OfsLST4(  0.0, 0.0, 0.0 );

  // STChamber2 *ST4 =
  //   new STChamber2( "ST4", pMother, RMST4, 
  // 		    gPos57, gPos58, gPos59, gPos60, gPos61, gPos62, 
  // 		    OfsLST4, 
  // 		    id57, id58, id59, id60, id61, id62, 
  // 		    matST4Gas, matST4Frame, matST4Box );
  // ST4->SetVisAttributes( &G4VisAttributes::Invisible,
  // 			 FrameAttrib, LayerAttrib,
  // 			 &G4VisAttributes::Invisible );
  
  // //////// Sensitive Detectors ////////
  // G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  // TrackerSD *trackerSD = new TrackerSD( "/spec/tr" );
  // SDMan->AddNewDetector( trackerSD );  
  // BT1->SetSensitiveDetector( trackerSD );
  // BT2->SetSensitiveDetector( trackerSD );

  // ST1->SetSensitiveDetector( trackerSD );
  // ST2->SetSensitiveDetector( trackerSD );
  // ST3->SetSensitiveDetector( trackerSD );
  // ST4->SetSensitiveDetector( trackerSD );
}


void DetectorConstructionSpec::
MakeTofCounters( G4VPhysicalVolume *pMother )
{
  // const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  // G4VisAttributes *TofScinAttrib =
  //   new G4VisAttributes( G4Colour::G4Colour( 1.0, 0.0, 1.0 ) );
  // G4VisAttributes *TofLGAttrib =
  //   new G4VisAttributes( G4Colour::G4Colour( 0.5, 0.0, 0.5 ) );
  // G4VisAttributes *TofPMTAttrib =
  //   new G4VisAttributes( G4Colour::G4Colour( 1.0, 1.0, 1.0 ) );
  // G4VisAttributes *TofPMTBoxAttrib =
  //   new G4VisAttributes( G4Colour::G4Colour( 1.0, 1.0, 1.0 ) );
  // G4VisAttributes *TofFrameAttrib =
  //   new G4VisAttributes( G4Colour::G4Colour( 0.5, 0.5, 0.5 ) );

  // G4int id1=geomMan.GetDetectorId("T1");
  // G4ThreeVector gPosT1=geomMan.GetGlobalPosition(id1);
  // G4double RotationAngleT1=geomMan.GetRotAngle2(id1);
  // G4RotationMatrix RMT1; RMT1.rotateZ( RotationAngleT1*degree);
  // G4ThreeVector OfsLT1( 0.0, 0.0, 0.0 );

  // T1Wall1 *T1 =
  //   new T1Wall1( "T1", pMother, RMT1, gPosT1, OfsLT1,
  // 		 matT1Scin, matT1LG, matT1PMT, matT1Box );
  // T1->SetVisAttributes( TofScinAttrib, TofLGAttrib, TofPMTAttrib,
  // 			&G4VisAttributes::Invisible );
  
  // G4int id2=geomMan.GetDetectorId("TOF");
  // G4ThreeVector gPosTof1=geomMan.GetGlobalPosition(id2);
  // G4double RotationAngleTof1=geomMan.GetRotAngle2(id2);
  // G4RotationMatrix RMTof1; RMTof1.rotateZ( RotationAngleTof1*degree);
  // G4ThreeVector OfsLTof1( 0.0, 0.0, 0.0 );

  // TofWall1 *Tof1 =
  //   new TofWall1( "TOF", pMother, RMTof1, gPosTof1, OfsLTof1,
  // 		  matTofScin, matTofLG, matTofPMT, matTofPMTBox, 
  // 		  matTofFrame, matTofBox );
  // Tof1->SetVisAttributes( TofScinAttrib, TofLGAttrib, TofPMTAttrib,
  // 			  TofPMTBoxAttrib, TofFrameAttrib,
  // 			  &G4VisAttributes::Invisible );

  // //////// Sensitive Detectors ////////
  // G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  // T1SD *t1SD = new T1SD( "/spec/t1" );
  // SDMan->AddNewDetector( t1SD );
  // T1->SetSensitiveDetector( t1SD );
  // TofSD *tofSD = new TofSD( "/spec/tof" );
  // SDMan->AddNewDetector( tofSD );
  // Tof1->SetSensitiveDetector( tofSD );
}


void DetectorConstructionSpec::
MakeVDetector( G4VPhysicalVolume *pMother )
  
{
  // const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  // G4VisAttributes *VDAttrib =
  //   new G4VisAttributes( G4Colour::G4Colour( 1.0, 1.0, 1.0 ) );
  
  // G4Box *solidVD = 
  //   new G4Box( "VD", 3500./2., 1./2., 2000./2. );

  // G4LogicalVolume *logVD =
  //   new G4LogicalVolume( solidVD, matVD, "VD", 0, 0, 0 );

  // G4int id1=geomMan.GetDetectorId("SpecVp1");
  // G4ThreeVector gPosVD1=geomMan.GetGlobalPosition(id1);
  // G4double RotationAngleVD1=geomMan.GetRotAngle2(id1);
  // G4RotationMatrix RMVD1; RMVD1.rotateZ( RotationAngleVD1*degree);
  // G4ThreeVector OfsLVD1( 0.0, 0.0, 0.0 );

  // // G4VPhysicalVolume *physVD1 =
  // //   new G4PVPlacement( G4Transform3D( RMVD1,
  // //                                     G4ThreeVector( gPosVD1.x()*mm,
  // //                                                    gPosVD1.y()*mm,
  // //                                                    gPosVD1.z()*mm ) ),
  // // 		       "VD1", logVD, pMother, false, id1 );

  // G4int id2=geomMan.GetDetectorId("SpecVp2");
  // G4ThreeVector gPosVD2=geomMan.GetGlobalPosition(id2);
  // G4double RotationAngleVD2=geomMan.GetRotAngle2(id2);
  // G4RotationMatrix RMVD2; RMVD2.rotateZ( RotationAngleVD2*degree);
  // G4ThreeVector OfsLVD2( 0.0, 0.0, 0.0 );

  // G4VPhysicalVolume *physVD2 =
  //   new G4PVPlacement( G4Transform3D( RMVD2,
  //                                     G4ThreeVector( gPosVD2.x()*mm,
  //                                                    gPosVD2.y()*mm,
  //                                                    gPosVD2.z()*mm ) ),
  // 		       "VD2", logVD, pMother, false, id2 );

  // logVD->SetVisAttributes(VDAttrib);

  // //////// Sensitive Detectors ////////
  // G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  // VDSD *vdSD = new VDSD( "/Sks/VD" );
  // SDMan->AddNewDetector( vdSD );
  // logVD->SetSensitiveDetector( vdSD );
}


G4bool DetectorConstructionSpec::
IsVolumeStopper( G4VPhysicalVolume *physVol ) const
{
  G4bool fSpec=DetectorConstruction::IsVolumeStopper( physVol );
  if(fSpec) return true;
  // if other magnets exsits, check them.

  return false;
}
