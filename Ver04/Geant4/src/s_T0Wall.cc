/*
  s_T0Wall.cc

  Simple T0 -> Just plate
  
  2015/12  K.Shirotori
*/

#include "s_T0Wall.hh"
#include "s_DetectorSize.hh"

#include "globals.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4ios.hh"

s_T0Wall::s_T0Wall( const G4String & Cname, 
		    G4LogicalVolume *pMother,
		    const G4RotationMatrix & rotMat, 
		    const G4ThreeVector & gPos, 
		    const G4ThreeVector & OffsetLocal,
		    G4Material *matScin,
		    G4Material *matBox )
  : Cname_(Cname)
{
  G4ThreeVector gOffset=rotMat*OffsetLocal;
  G4ThreeVector T0gPos=gPos+gOffset;

  //  G4cout << "OffL " << OffsetLocal << "  gOffset " << gOffset
  //	 << " gPos " << gPos << " T0gPos " << T0gPos
  //	 << G4endl;

  G4double SegsizeX = T0_SegsizeX;
  G4double SegsizeY = T0_SegsizeY;
  G4double SegsizeZ = T0_SegsizeZ;

  G4double BoxsizeX = SegsizeX + 1.0*mm;
  G4double BoxsizeY = SegsizeY + 1.0*mm;
  G4double BoxsizeZ = SegsizeZ + 1.0*mm;

  G4Box *solidBox =
    new G4Box( Cname_+"Box", BoxsizeX/2., BoxsizeZ/2., BoxsizeY/2. );
  G4Box *solidSeg =
    new G4Box( Cname_+"Segment", SegsizeX/2., SegsizeZ/2., SegsizeY/2. );

  logBox = 
    new G4LogicalVolume( solidBox, matBox, Cname_+"Box", 0, 0, 0 );
  logScin = 
    new G4LogicalVolume( solidSeg, matScin, Cname_+"Scin", 0, 0, 0 );

  G4VPhysicalVolume *physBox =
    new G4PVPlacement( G4Transform3D( rotMat,
                                      G4ThreeVector( T0gPos.x()*mm,
                                                     T0gPos.y()*mm,
                                                     T0gPos.z()*mm ) ),
		       logBox, Cname_+"Box", pMother, false, 0 );

  for( G4int id=1; id<=T0_SegNum; ++id ){
    G4double ofsScinX=-T0_SegSpacing*(T0_SegNum/2.+0.5-id);
    new G4PVPlacement( 0, G4ThreeVector( ofsScinX, 0.0*mm, 0.0*mm ),
		       Cname_+"Seg", logScin, physBox, false, id );
  }
}

void s_T0Wall::SetVisAttributes( const G4VisAttributes *attScin,
				 const G4VisAttributes *attBox )
{
  if(attScin) logScin->SetVisAttributes(attScin);
  if(attBox) logBox->SetVisAttributes(attBox);
}

void s_T0Wall::SetSensitiveDetector( G4VSensitiveDetector *SD )
{
  logScin->SetSensitiveDetector(SD);
}
