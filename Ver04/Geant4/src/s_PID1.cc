/*
  s_PID1.hh
  
  Simple PID1 -> Just plate and box

  2015/12  K.Shirotori
*/

#include "s_PID1.hh"
#include "s_DetectorSize.hh"

#include "globals.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4ios.hh"

s_PID1::s_PID1( const G4String & Cname, 
		G4LogicalVolume *pMother,
		const G4RotationMatrix & rotMat, 
		const G4ThreeVector & gPos, 
		const G4ThreeVector & OffsetLocal,
		G4Material *matSeg,
		G4Material *matArea,
		G4Material *matFrame,
		G4Material *matBox )
  : Cname_(Cname)
{
  G4ThreeVector gOffset=rotMat*OffsetLocal;
  G4ThreeVector PID1gPos=gPos+gOffset;

  //  G4cout << "OffL " << OffsetLocal << "  gOffset " << gOffset
  //	 << " gPos " << gPos << " PID1gPos " << PID1gPos
  //	 << G4endl;

  //Frame size  
  G4double BoxsizeX = PID1_BoxsizeX;
  G4double BoxsizeY = PID1_BoxsizeY;
  G4double BoxsizeZ = PID1_BoxsizeZ;
  G4double FrameT1 =  PID1_FrameThick1;

  //Radiator size
  G4double SegsizeX =  BoxsizeX-2.*FrameT1;
  G4double SegsizeY =  BoxsizeY-2.*FrameT1;
  G4double SegsizeZ =  PID1_SegsizeZ;

  //Full box size
  G4double BoxsizeX0 = BoxsizeX + 1.0*mm;
  G4double BoxsizeY0 = BoxsizeY + 1.0*mm;
  G4double BoxsizeZ0 = BoxsizeZ + 1.0*mm;

  //Full box
  G4Box *solidBox =
    new G4Box( Cname_+"Box", BoxsizeX0/2., BoxsizeZ0/2., BoxsizeY0/2. );

  //Frame box
  G4Box *solidBase1 =
    new G4Box( Cname_+"Base1", BoxsizeX/2., BoxsizeZ/2., BoxsizeY/2. );
  G4Box *solidBase2 =
    new G4Box( Cname_+"Base2", BoxsizeX/2.-FrameT1, BoxsizeZ/2.+0.001*mm, BoxsizeY/2.-FrameT1 );

  G4SubtractionSolid *solidFrame =
    new G4SubtractionSolid( Cname_+"Frame", solidBase1, solidBase2,
			    0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ) );
  
  G4Box *solidSeg =
    new G4Box( Cname_+"Seg", SegsizeX/2., SegsizeZ/2., SegsizeY/2. );

  logBox = 
    new G4LogicalVolume( solidBox, matBox, Cname_+"Box", 0, 0, 0 );
  logFrame = 
    new G4LogicalVolume( solidFrame, matFrame, Cname_+"Frame", 0, 0, 0 );
  logArea = 
    new G4LogicalVolume( solidBase2, matArea, Cname_+"Area", 0, 0, 0 );
  logSeg = 
    new G4LogicalVolume( solidSeg, matSeg, Cname_+"Seg", 0, 0, 0 );
  
  G4VPhysicalVolume *physBox =
    new G4PVPlacement( G4Transform3D( rotMat,
                                      G4ThreeVector( PID1gPos.x()*mm,
                                                     PID1gPos.y()*mm,
                                                     PID1gPos.z()*mm ) ),
  		       logBox, Cname_+"Box", pMother, false, 0 );
  
  G4VPhysicalVolume *physFrame =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       Cname_+"Frame", logFrame, physBox, false, 0 );
  
  G4VPhysicalVolume *physArea =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       Cname_+"Area", logArea, physBox, false, 0 );
  
  G4VPhysicalVolume *physSeg =
    new G4PVPlacement( 0, G4ThreeVector( (BoxsizeX-SegsizeX)/2.-FrameT1, (BoxsizeZ-SegsizeZ)/2., 0.0*mm ),
                       Cname_+"Seg", logSeg, physArea, false, 1 );
}

void s_PID1::SetVisAttributes( const G4VisAttributes *attSeg,
			       const G4VisAttributes *attFrame,
			       const G4VisAttributes *attArea,
			       const G4VisAttributes *attBox )
{
  if(attSeg)   logSeg->SetVisAttributes(attSeg);
  if(attFrame) logFrame->SetVisAttributes(attFrame);
  if(attArea)  logArea->SetVisAttributes(attArea);
  if(attBox)   logBox->SetVisAttributes(attBox);
}

void s_PID1::SetSensitiveDetector( G4VSensitiveDetector *SD )
{
  logSeg->SetSensitiveDetector(SD);
}
