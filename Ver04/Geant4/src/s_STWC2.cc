/*
  s_STWC2.cc

  Simple Magnet Downstream WC 2
  -> Just plate + Frame

  2015/12  K.Shirotori
*/


#include "s_STWC2.hh"
#include "s_DetectorSize.hh"
#include "ConfMan.hh"
#include <cmath>

const G4double LzLayer1= -120.0;
const G4double LzLayer2=  -90.0;
const G4double LzLayer3=  -60.0;
const G4double LzLayer4=  -30.0;
const G4double LzLayer5=    0.0;
const G4double LzLayer6=   30.0;
const G4double LzLayer7=   60.0;
const G4double LzLayer8=   90.0;
const G4double LzLayer9=  120.0;

const G4double MaxDispl= 0.01; /* 10um */

const G4double Deg2Rad = acos(-1.)/180.;
const G4double Rad2Deg = 180./acos(-1.);

#include "globals.hh"
#include "G4Box.hh"
#include "G4Para.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4ios.hh"

s_STWC2::s_STWC2( const G4String & Cname,
		  G4LogicalVolume *pMother,
		  const G4RotationMatrix &rotMat,
		  const G4ThreeVector &gPos1, 
		  const G4ThreeVector &gPos2,
		  const G4ThreeVector &gPos3, 
		  const G4ThreeVector &gPos4,
		  const G4ThreeVector &gPos5,
		  const G4ThreeVector &gPos6,
		  const G4ThreeVector &gPos7,
		  const G4ThreeVector &gPos8,
		  const G4ThreeVector &gPos9,
		  const G4ThreeVector &OffsetLocal,
		  G4int id1, G4int id2, G4int id3, 
		  G4int id4, G4int id5, G4int id6,
		  G4int id7, G4int id8, G4int id9,
		  G4Material *matWire, 
		  G4Material *matFrame, 
		  G4Material *matGas, 
		  G4Material *matBox )
  : Cname_(Cname)
{
  ConfMan *confMan = ConfMan::GetConfManager();

  G4ThreeVector gPosCent=gPos5;
  G4ThreeVector gOffset=rotMat*OffsetLocal;
  G4ThreeVector DCgPos=gPosCent+gOffset;
  
  G4double LaysizeX     = ST2_LaysizeX;
  G4double LaysizeUV    = ST2_LaysizeUV;
  G4double LaysizeY     = ST2_LaysizeY;
  G4double Laythickness = ST2_LaysizeZ;

  G4double tilted    = ST2_TiltAngle;  
  G4double tilted_BT = ST2_TiltAngle_BT;  
  G4double BaseThick = ST2_BoxThick*9.0 + 20.0*mm;
  
  G4double width = LaysizeUV;//+LaysizeY*tan(tilted);
  G4double FrameSize = width;//+LaysizeY*tan(tilted);
  G4double FrameTH = ST2_FrameThickH;
  G4double FrameTV = ST2_FrameThickV;

  G4double sizeX_BT;
  G4double sizeUV_BT;
  G4double posX_BT;
  G4double posUV_BT;

  if( confMan->WCHoleFlag() ){
    sizeX_BT  = ST2_BTsize;
    sizeUV_BT = sizeX_BT/cos(tilted_BT);
    posX_BT  = ST2_BTpos;
    posUV_BT = ST2_BTpos;
  }
  else{
    sizeX_BT  = 0.0001*mm;
    sizeUV_BT = sizeX_BT/cos(tilted_BT);
    posX_BT  = 0.0*mm;
    posUV_BT = 0.0*mm;
  }

  G4Box *solidBox  = 
    new G4Box( Cname_+"Box", (FrameSize+1.0*mm)/2.+FrameTH, (BaseThick+1.0*mm)/2., (LaysizeY+1.0*mm)/2.+FrameTV );
  G4Box *solidBase  = 
    new G4Box( Cname_+"Base", FrameSize/2.+FrameTH, BaseThick/2., LaysizeY/2.+FrameTV );
  G4Box *solidArea  = 
    new G4Box( Cname_+"Area", FrameSize/2., BaseThick/2.+0.001*mm, LaysizeY/2. );

  G4Box *solidLayerX_0 = 
    new G4Box( Cname_+"LayerX_0", LaysizeX/2., Laythickness/2., LaysizeY/2. );
  G4Box *solidLayerX_BT = 
    new G4Box( Cname_+"LayerX_BT", sizeX_BT/2., Laythickness/2.+0.001*mm, LaysizeY/2. );

  G4Para *solidLayerU_0 = 
    new G4Para( Cname_+"LayerU_0", width/2., Laythickness/2., LaysizeY/2., 
		0.0*degree, (-1.)*tilted, 0.0*degree );
  G4Para *solidLayerU_BT = 
    new G4Para( Cname_+"LayerU", sizeUV_BT/2., Laythickness/2.+0.001*mm, LaysizeY/2., 
		0.0*degree, (-1.)*tilted_BT, 0.0*degree );

  G4Para *solidLayerV_0 = 
    new G4Para( Cname_+"LayerV_0", width/2., Laythickness/2., LaysizeY/2., 
		0.0*degree,  tilted, 0.0*degree );
  G4Para *solidLayerV_BT = 
    new G4Para( Cname_+"LayerV_BY", sizeUV_BT/2., Laythickness/2.+0.001*mm, LaysizeY/2., 
		0.0*degree,  tilted_BT, 0.0*degree );

  G4SubtractionSolid *solidLayerX =
    new G4SubtractionSolid( Cname_+"LayerX", solidLayerX_0, solidLayerX_BT,
			    0, G4ThreeVector( posX_BT, 0.0*mm, 0.0*mm) );
  G4SubtractionSolid *solidLayerU =
    new G4SubtractionSolid( Cname_+"LayerU", solidLayerU_0, solidLayerU_BT,
			    0, G4ThreeVector( posUV_BT, 0.0*mm, 0.0*mm) );
  G4SubtractionSolid *solidLayerV =
    new G4SubtractionSolid( Cname_+"LayerV", solidLayerV_0, solidLayerV_BT,
			    0, G4ThreeVector( posUV_BT, 0.0*mm, 0.0*mm) );

  G4SubtractionSolid *solidFrame =
    new G4SubtractionSolid( Cname_+"Frame", solidBase, solidArea );
  
  logBox   = 
    new G4LogicalVolume( solidBox,  matBox, Cname_+"Box",   0, 0, 0 );
  logFrame   = 
    new G4LogicalVolume( solidFrame,  matFrame, Cname_+"Frame",   0, 0, 0 );
  logArea   = 
    new G4LogicalVolume( solidArea,  matGas, Cname_+"Area",   0, 0, 0 );
  logLayerX = 
    new G4LogicalVolume( solidLayerX, matWire, Cname_+"LayerX", 0, 0, 0 );
  logLayerU = 
    new G4LogicalVolume( solidLayerU, matWire, Cname_+"LayerU", 0, 0, 0 );
  logLayerV = 
    new G4LogicalVolume( solidLayerV, matWire, Cname_+"LayerV", 0, 0, 0 );
  
  G4VPhysicalVolume *physBox =
    new G4PVPlacement( G4Transform3D( rotMat,
				      G4ThreeVector( DCgPos.x()*mm,
						     DCgPos.y()*mm,
						     DCgPos.z()*mm ) ),
		       logBox, Cname_+"Box", pMother, false, 0 );

  G4VPhysicalVolume *physFrame =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       Cname_+"Frame", logFrame, physBox, false, 0 );

  G4VPhysicalVolume *physArea =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       Cname_+"Area", logArea, physBox, false, 0 );
 
  // Consitency should be checked.
  G4RotationMatrix InvMat=rotMat.inverse();
  G4ThreeVector lPos1=InvMat*(gPos1-gPosCent);
  G4ThreeVector lPos2=InvMat*(gPos2-gPosCent);
  G4ThreeVector lPos3=InvMat*(gPos3-gPosCent);
  G4ThreeVector lPos4=InvMat*(gPos4-gPosCent);
  G4ThreeVector lPos5=InvMat*(gPos5-gPosCent);
  G4ThreeVector lPos6=InvMat*(gPos6-gPosCent);
  G4ThreeVector lPos7=InvMat*(gPos7-gPosCent);
  G4ThreeVector lPos8=InvMat*(gPos8-gPosCent);
  G4ThreeVector lPos9=InvMat*(gPos9-gPosCent);
  
  if(fabs(lPos1.x())>MaxDispl || fabs(lPos1.y()+LzLayer1)>MaxDispl ){
    G4cout << Cname_ << "Layer1: Geometry Error!:" << "(0," << -LzLayer1
	   << ",0) --> " << lPos1 << G4endl;
    G4Exception("Geometry Error");
  }
  if(fabs(lPos2.x())>MaxDispl || fabs(lPos2.y()+LzLayer2)>MaxDispl ){
    G4cout << Cname_ << "Layer2: Geometry Error!:" << "(0," << -LzLayer2
	   << ",0) --> " << lPos2 << G4endl;
    G4Exception("Geometry Error");
  }
  if(fabs(lPos3.x())>MaxDispl || fabs(lPos3.y()+LzLayer3)>MaxDispl ){
    G4cout << Cname_ << "Layer3: Geometry Error!:" << "(0," << -LzLayer3
	   << ",0) --> " << lPos3 << G4endl;
    G4Exception("Geometry Error");
  }
  if(fabs(lPos4.x())>MaxDispl || fabs(lPos4.y()+LzLayer4)>MaxDispl ){
    G4cout << Cname_ << "Layer4: Geometry Error!:" << "(0," << -LzLayer4
	   << ",0) --> " << lPos4 << G4endl;
    G4Exception("Geometry Error");
  }
  if(fabs(lPos5.x())>MaxDispl || fabs(lPos5.y()+LzLayer5)>MaxDispl ){
    G4cout << Cname_ << "Layer5: Geometry Error!:" << "(0," << -LzLayer5
	   << ",0) --> " << lPos5 << G4endl;
    G4Exception("Geometry Error");
  }
  if(fabs(lPos6.x())>MaxDispl || fabs(lPos6.y()+LzLayer6)>MaxDispl ){
    G4cout << Cname_ << "Layer6: Geometry Error!:" << "(0," << -LzLayer6
	   << ",0) --> " << lPos6 << G4endl;
    G4Exception("Geometry Error");
  }
  if(fabs(lPos7.x())>MaxDispl || fabs(lPos7.y()+LzLayer7)>MaxDispl ){
    G4cout << Cname_ << "Layer7: Geometry Error!:" << "(0," << -LzLayer7
	   << ",0) --> " << lPos7 << G4endl;
    G4Exception("Geometry Error");
  }
  if(fabs(lPos8.x())>MaxDispl || fabs(lPos8.y()+LzLayer8)>MaxDispl ){
    G4cout << Cname_ << "Layer8: Geometry Error!:" << "(0," << -LzLayer8
	   << ",0) --> " << lPos8 << G4endl;
    G4Exception("Geometry Error");
  }
  if(fabs(lPos9.x())>MaxDispl || fabs(lPos9.y()+LzLayer9)>MaxDispl ){
    G4cout << Cname_ << "Layer9: Geometry Error!:" << "(0," << -LzLayer9
	   << ",0) --> " << lPos9 << G4endl;
    G4Exception("Geometry Error");
  }
    
  // Layer1
  G4VPhysicalVolume *physLayer1 =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, -LzLayer1*mm, 0.0*mm ),
  		       Cname_+"Layer1", logLayerX, physArea, false, id1 );
  // Layer2
  G4VPhysicalVolume *physLayer2 =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, -LzLayer2*mm, 0.0*mm ),
		       Cname_+"Layer2", logLayerU, physArea, false, id2 );
  // Layer3
  G4VPhysicalVolume *physLayer3 =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, -LzLayer3*mm, 0.0*mm ),
  		       Cname_+"Layer3", logLayerV, physArea, false, id3 );
  // Layer4
  G4VPhysicalVolume *physLayer4 =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, -LzLayer4*mm, 0.0*mm ),
  		       Cname_+"Layer4", logLayerU, physArea, false, id4 );
  // Layer5
  G4VPhysicalVolume *physLayer5 =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, -LzLayer5*mm, 0.0*mm ),
  		       Cname_+"Layer5", logLayerV, physArea, false, id5 );
  // Layer6
  G4VPhysicalVolume *physLayer6 =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, -LzLayer6*mm, 0.0*mm ),
  		       Cname_+"Layer6", logLayerX, physArea, false, id6 );
  // Layer7
  G4VPhysicalVolume *physLayer7 =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, -LzLayer7*mm, 0.0*mm ),
  		       Cname_+"Layer7", logLayerU, physArea, false, id7 );
  // Layer8
  G4VPhysicalVolume *physLayer8 =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, -LzLayer8*mm, 0.0*mm ),
  		       Cname_+"Layer8", logLayerV, physArea, false, id8 );
  // Layer9
  G4VPhysicalVolume *physLayer9 =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, -LzLayer9*mm, 0.0*mm ),
  		       Cname_+"Layer9", logLayerX, physArea, false, id9 );
}

void s_STWC2::SetVisAttributes( const G4VisAttributes *attLayer,
				const G4VisAttributes *attFrame,
				const G4VisAttributes *attArea,
				const G4VisAttributes *attBox )
{
  if(attLayer) logLayerX->SetVisAttributes(attLayer); 
  if(attLayer) logLayerU->SetVisAttributes(attLayer); 
  if(attLayer) logLayerV->SetVisAttributes(attLayer); 
  if(attFrame) logFrame->SetVisAttributes(attFrame); 
  if(attArea)  logArea->SetVisAttributes(attArea); 
  if(attBox)   logBox->SetVisAttributes(attBox); 
}

void s_STWC2::SetSensitiveDetector( G4VSensitiveDetector *SD )
{
  logLayerX->SetSensitiveDetector(SD);
  logLayerU->SetSensitiveDetector(SD);
  logLayerV->SetSensitiveDetector(SD);
}
