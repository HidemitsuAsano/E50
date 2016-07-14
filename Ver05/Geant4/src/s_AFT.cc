/*
  s_AFT.cc

  Simple AFT
  -> Just plate + Frame

  2015/12  K.Shirotori
*/


#include "s_AFT.hh"
#include "s_DetectorSize.hh"
#include <cmath>

const G4double LzLayer1= -0.5;
const G4double LzLayer2=  0.5;

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

s_AFT::s_AFT( const G4String & Cname,
	      G4LogicalVolume *pMother,
	      const G4RotationMatrix &rotMat,
	      const G4ThreeVector &gPos1, 
	      const G4ThreeVector &gPos2,
	      const G4ThreeVector &OffsetLocal,
	      G4int id1, 
	      G4int id2, 
	      G4Material *matScin, 
	      G4Material *matFrame,
	      G4Material *matArea,
	      G4Material *matBox )
  : Cname_(Cname)
{
  G4ThreeVector gPosCent=(gPos1+gPos2)/2.;
  G4ThreeVector gOffset=rotMat*OffsetLocal;
  G4ThreeVector DCgPos=gPosCent+gOffset;
  
  G4double LaysizeX     = AFT_LaysizeX;
  G4double LaysizeY     = AFT_LaysizeY;
  G4double Laythickness = AFT_LaysizeZ;

  G4double tilted = AFT_TiltAngle;
  G4double BaseThick = AFT_BoxThick*2.0 + 20.0*mm;

  G4double width = LaysizeX;//+LaysizeY*tan(tilted);
  G4double FrameSize = width;
  G4double FrameT = AFT_FrameThick;

  G4Box *solidBox  = 
    new G4Box( Cname_+"Box", (FrameSize+1.0*mm)/2.+FrameT, (LaysizeY+1.0*mm)/2.+FrameT, (BaseThick+0.5*mm)/2. );
  G4Box *solidBase  = 
    new G4Box( Cname_+"Base", FrameSize/2.+FrameT, LaysizeY/2.+FrameT, BaseThick/2. );
  G4Box *solidArea  = 
    new G4Box( Cname_+"Area", FrameSize/2., LaysizeY/2., BaseThick/2. );
  
  G4Box *solidLayerX = 
    new G4Box( Cname_+"LayerX", LaysizeX/2., LaysizeY/2., Laythickness/2. );
  
  G4SubtractionSolid *solidFrame =
    new G4SubtractionSolid( Cname_+"Frame", solidBase, solidArea );
  
  logBox   = 
    new G4LogicalVolume( solidBox,  matBox, Cname_+"Box",   0, 0, 0 );
  logFrame   = 
    new G4LogicalVolume( solidFrame,  matFrame, Cname_+"Frame",   0, 0, 0 );
  logArea   = 
    new G4LogicalVolume( solidArea,  matArea, Cname_+"Area",   0, 0, 0 );
  logLayerX = 
    new G4LogicalVolume( solidLayerX, matScin, Cname_+"LayerX", 0, 0, 0 );
  
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

  if(fabs(lPos1.x())>MaxDispl || fabs(lPos1.y())>MaxDispl || fabs(lPos1.z()-LzLayer1)>MaxDispl ){
    G4cout << Cname_ << "Layer1: Geometry Error!:" << "(0,0," << LzLayer1
  	   << ") --> " << lPos1 << G4endl;
    G4Exception("Geometry Error");
  }
  if(fabs(lPos2.x())>MaxDispl || fabs(lPos2.y())>MaxDispl || fabs(lPos2.z()-LzLayer2)>MaxDispl ){
    G4cout << Cname_ << "Layer2: Geometry Error!:" << "(0,0," << LzLayer2
  	   << ") --> " << lPos2 << G4endl;
    G4Exception("Geometry Error");
  }
    
  // Layer1
  G4VPhysicalVolume *physLayer1 =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, LzLayer1*mm ),
		       Cname_+"Layer1", logLayerX, physArea, false, id1 );
  // Layer2
  G4VPhysicalVolume *physLayer2 =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, LzLayer2*mm ),
		       Cname_+"Layer2", logLayerX, physArea, false, id2 );
}

void s_AFT::SetVisAttributes( const G4VisAttributes *attLayer,
			      const G4VisAttributes *attFrame,
			      const G4VisAttributes *attArea,
			      const G4VisAttributes *attBox )
{
  if(attLayer) logLayerX->SetVisAttributes(attLayer); 
  if(attFrame) logFrame->SetVisAttributes(attFrame); 
  if(attArea)  logArea->SetVisAttributes(attArea); 
  if(attBox)   logBox->SetVisAttributes(attBox); 
}

void s_AFT::SetSensitiveDetector( G4VSensitiveDetector *SD )
{
  logLayerX->SetSensitiveDetector(SD);
}
