/*
  s_RICH.hh
  
  Simple RICH -> Just plate and box

  2015/12  K.Shirotori
*/

#include "s_RICH.hh"
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

s_RICH::s_RICH( const G4String & Cname, 
		G4LogicalVolume *pMother,
		const G4RotationMatrix & rotMat, 
		const G4ThreeVector & gPos, 
		const G4ThreeVector & OffsetLocal,
		G4Material *matAerogel,
		G4Material *matGas,
		G4Material *matFrame,
		G4Material *matBox )
  : Cname_(Cname)
{
  G4ThreeVector gOffset=rotMat*OffsetLocal;
  G4ThreeVector RICHgPos=gPos+gOffset;

  //  G4cout << "OffL " << OffsetLocal << "  gOffset " << gOffset
  //	 << " gPos " << gPos << " RICHgPos " << RICHgPos
  //	 << G4endl;

  //AC radiator size
  G4double SegsizeX = RICH_SegsizeX;
  G4double SegsizeY = RICH_SegsizeY;
  G4double SegsizeZ = RICH_SegsizeZ;

  //Frame size  
  G4double BoxsizeX = RICH_BoxsizeX;
  G4double BoxsizeY = RICH_BoxsizeY;
  G4double BoxsizeZ = RICH_BoxsizeZ;
  G4double FrameT1 =  RICH_FrameThick1;
  G4double FrameT2 =  RICH_FrameThick2;

  //Radiator box size
  G4double angle     = RICH_FrameAngle;  
  G4double BoxsizeX1 = SegsizeX + 2.*FrameT2;
  G4double BoxsizeY1 = SegsizeY + 2.*FrameT2;
  G4double BoxsizeZ1 = RICH_BoxsizeZ1;
  G4double BoxsizeX2 = BoxsizeX1 + BoxsizeZ1*tan(angle);
  G4double BoxsizeY2 = BoxsizeY1 + BoxsizeZ1*tan(angle);
  G4double PosRB = (BoxsizeZ + BoxsizeZ1)/2.;
  G4double PosGas = (BoxsizeZ + BoxsizeZ1)/2.;

  //Full box size
  G4double BoxsizeX0 = BoxsizeX + 1.0*mm;
  G4double BoxsizeY0 = BoxsizeY + 1.0*mm;
  G4double BoxsizeZ0 = BoxsizeZ + BoxsizeZ1 + 1.0*mm;

  //Full box
  G4Box *solidBox =
    new G4Box( Cname_+"Box", BoxsizeX0/2., BoxsizeY0/2., BoxsizeZ0/2. );

  //Rectangle box
  G4Box *solidBaseA1 =
    new G4Box( Cname_+"BaseA1", BoxsizeX/2., BoxsizeY/2., BoxsizeZ/2. );
  G4Box *solidBaseA2 =
    new G4Box( Cname_+"BaseA2", BoxsizeX/2.-FrameT1, BoxsizeY/2.-FrameT1, BoxsizeZ/2.-FrameT1 );
  G4Box *solidBaseA3 =
    new G4Box( Cname_+"BaseA3", BoxsizeX2/2.-FrameT2, BoxsizeY2/2.-FrameT2, FrameT1/2. );

  G4SubtractionSolid *solidBaseA4 =
    new G4SubtractionSolid( Cname_+"BaseA4", solidBaseA1, solidBaseA2,
			    0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ) );
  G4SubtractionSolid *solidBaseA =
    new G4SubtractionSolid( Cname_+"BaseA", solidBaseA4, solidBaseA3,
			    0, G4ThreeVector( 0.0*mm, 0.0*mm, (-1.)*(BoxsizeZ-FrameT1)/2.) );
  G4UnionSolid *solidGas1 =
    new G4UnionSolid( Cname_+"Gas1", solidBaseA2, solidBaseA3,
		      0, G4ThreeVector( 0.0*mm, 0.0*mm, (-1.)*(BoxsizeZ-FrameT1)/2.) );

  //Radiator box
  G4Trd *solidBaseB1 =
    new G4Trd( Cname_+"BaseB1", BoxsizeX1/2., BoxsizeX2/2., 
  	       BoxsizeY1/2., BoxsizeY2/2., BoxsizeZ1/2.);
  G4Trd *solidBaseB2 =
    new G4Trd( Cname_+"BaseB2", BoxsizeX1/2.-FrameT2, BoxsizeX2/2.-FrameT2, 
	       BoxsizeY1/2.-FrameT2, BoxsizeY2/2.-FrameT2, BoxsizeZ1/2.-FrameT2 );

  G4Box *solidBaseB3 =
    new G4Box( Cname_+"BaseB3", BoxsizeX2/2.-FrameT2, BoxsizeY2/2.-FrameT2, FrameT2/2. );
  G4Box *solidBaseB4 =
    new G4Box( Cname_+"BaseB4", BoxsizeX1/2.-FrameT2, BoxsizeY1/2.-FrameT2, FrameT2/2. );

  G4SubtractionSolid *solidBaseB5 =
    new G4SubtractionSolid( Cname_+"BaseB5", solidBaseB1, solidBaseB2,
  			    0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ) );
  G4SubtractionSolid *solidBaseB6 =
    new G4SubtractionSolid( Cname_+"BaseB6", solidBaseB5, solidBaseB3,
  			    0, G4ThreeVector( 0.0*mm, 0.0*mm, (BoxsizeZ1-FrameT2)/2.) );
  G4SubtractionSolid *solidBaseB =
    new G4SubtractionSolid( Cname_+"BaseB", solidBaseB6, solidBaseB4,
  			    0, G4ThreeVector( 0.0*mm, 0.0*mm, (-1.)*(BoxsizeZ1-FrameT2)/2.) );

  G4UnionSolid *solidGas2A =
    new G4UnionSolid( Cname_+"Gas2A", solidBaseB2, solidBaseB3,
		      0, G4ThreeVector( 0.0*mm, 0.0*mm, (BoxsizeZ1-FrameT2)/2.) );
  G4UnionSolid *solidGas2 =
    new G4UnionSolid( Cname_+"Gas2", solidGas2A, solidBaseB4,
		      0, G4ThreeVector( 0.0*mm, 0.0*mm, (-1.)*(BoxsizeZ1-FrameT2)/2.) );
			    
  G4RotationMatrix RMx; 
  RMx.rotateX( 0.*degree );
  
  G4UnionSolid *solidFrame =
    new G4UnionSolid( Cname_+"Frame", solidBaseA, solidBaseB,
  		      G4Transform3D( RMx, G4ThreeVector( 0.0*mm, 0.0*mm, (-1.)*PosRB)) );
  G4UnionSolid *solidGas =
    new G4UnionSolid( Cname_+"Gas", solidGas1, solidGas2,
  		      G4Transform3D( RMx, G4ThreeVector( 0.0*mm, 0.0*mm, (-1.)*PosGas)) );

  G4Box *solidAerogel =
    new G4Box( Cname_+"Aerogel", SegsizeX/2., SegsizeY/2., SegsizeZ/2. );

  logBox = 
    new G4LogicalVolume( solidBox, matBox, Cname_+"Box", 0, 0, 0 );
  logFrame = 
    new G4LogicalVolume( solidFrame, matFrame, Cname_+"Frame", 0, 0, 0 );
  logGas = 
    new G4LogicalVolume( solidGas, matGas, Cname_+"Gas", 0, 0, 0 );
  logAerogel = 
    new G4LogicalVolume( solidAerogel, matAerogel, Cname_+"Aerogel", 0, 0, 0 );
  
  G4VPhysicalVolume *physBox =
    new G4PVPlacement( G4Transform3D( rotMat,
                                      G4ThreeVector( RICHgPos.x()*mm,
                                                     RICHgPos.y()*mm,
                                                     RICHgPos.z()*mm ) ),
  		       logBox, Cname_+"Box", pMother, false, 0 );
  
  G4VPhysicalVolume *physFrame =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, BoxsizeZ1/2. ),
                       Cname_+"Frame", logFrame, physBox, false, 0 );

  G4VPhysicalVolume *physGas =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, BoxsizeZ1/2. ),
                       Cname_+"Gas", logGas, physBox, false, 0 );

  G4VPhysicalVolume *physAerogel =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, (-1.)*((BoxsizeZ - SegsizeZ)/2.+ BoxsizeZ1) ),
                       Cname_+"Aerogel", logAerogel, physGas, false, 1 );
}

void s_RICH::SetVisAttributes( const G4VisAttributes *attAerogel,
			       const G4VisAttributes *attFrame,
			       const G4VisAttributes *attGas,
			       const G4VisAttributes *attBox )
{
  if(attAerogel) logAerogel->SetVisAttributes(attAerogel);
  if(attFrame)   logFrame->SetVisAttributes(attFrame);
  if(attGas)     logGas->SetVisAttributes(attGas);
  if(attBox)     logBox->SetVisAttributes(attBox);
}

void s_RICH::SetSensitiveDetector( G4VSensitiveDetector *SD )
{
  logAerogel->SetSensitiveDetector(SD);
}
