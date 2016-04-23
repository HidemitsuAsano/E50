/*
  s_PAD.cc

  Simple PAD -> Just plate
  
  2015/12  K.Shirotori
*/

#include "s_PAD.hh"
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

s_PAD::s_PAD( const G4String & Cname, 
	      G4LogicalVolume *pMother,
	      const G4RotationMatrix & rotMat, 
	      const G4ThreeVector & gPos, 
	      const G4ThreeVector & OffsetLocal,
	      G4Material *matScin,
	      G4Material *matBox )
  : Cname_(Cname)
{
  G4ThreeVector gOffset=rotMat*OffsetLocal;
  G4ThreeVector PADgPos=gPos+gOffset;

  //  G4cout << "OffL " << OffsetLocal << "  gOffset " << gOffset
  //	 << " gPos " << gPos << " PADgPos " << PADgPos
  //	 << G4endl;

  G4double SegsizeX  = PAD_SegsizeX1;
  G4double SegsizeX2 = PAD_SegsizeX2;
  G4double SegsizeZ  = PAD_SegsizeZ;

  G4double BoxsizeX  = SegsizeX;
  G4double BoxsizeX2 = BoxsizeX/2.;
  G4double BoxsizeZ = SegsizeZ + 1.0*mm;

  G4double SegsizeY_A = PAD_SegsizeY_A;
  G4double SegsizeY_B = PAD_SegsizeY_B;
  G4double SegsizeY_C = PAD_SegsizeY_C;
  G4double SegsizeY_D = PAD_SegsizeY_D;
  G4double SegsizeY_E = PAD_SegsizeY_E;

  G4double OffsetA =  PAD_Offset_A;
  G4double OffsetB =  PAD_Offset_B;
  G4double OffsetC = (SegsizeY_A - SegsizeY_C)/2.;
  G4double OffsetD = (SegsizeY_A - SegsizeY_D)/2.;
  G4double OffsetE = (SegsizeY_A - SegsizeY_E)/2.;

  G4double BoxsizeY_A = SegsizeY_A + 1.0*mm;
  G4double BoxsizeY_B = SegsizeY_B + 1.0*mm;
  G4double BoxsizeY_C = SegsizeY_C + 1.0*mm;
  G4double BoxsizeY_D = SegsizeY_D + 1.0*mm;
  G4double BoxsizeY_E = SegsizeY_E + 1.0*mm;

  //PAD Box
  G4Box *solidBoxA =
    new G4Box( Cname_+"BoxA", BoxsizeX/2., BoxsizeZ/2., BoxsizeY_A/2. );
  G4Box *solidBoxB =
    new G4Box( Cname_+"BoxB", BoxsizeX/2., BoxsizeZ/2., BoxsizeY_B/2. );
  G4Box *solidBoxC =
    new G4Box( Cname_+"BoxC", BoxsizeX/2., BoxsizeZ/2., BoxsizeY_C/2. );
  G4Box *solidBoxD =
    new G4Box( Cname_+"BoxD", BoxsizeX2/2., BoxsizeZ/2., BoxsizeY_D/2. );
  G4Box *solidBoxE =
    new G4Box( Cname_+"BoxE", BoxsizeX2/2., BoxsizeZ/2., BoxsizeY_E/2. );

  G4UnionSolid *solidBox1 =
    new G4UnionSolid( Cname_+"Box1", solidBoxA, solidBoxA,
  		      0, G4ThreeVector( BoxsizeX*(PAD_SegNum/2.+0.5-5.0), 0.0*mm, -OffsetA ) );
  G4UnionSolid *solidBox2 =
    new G4UnionSolid( Cname_+"Box2", solidBox1, solidBoxA,
  		      0, G4ThreeVector( BoxsizeX*(PAD_SegNum/2.+0.5-7.0), 0.0*mm, -OffsetA ) );
  G4UnionSolid *solidBox3 =
    new G4UnionSolid( Cname_+"Box3", solidBox2, solidBoxB,
  		      0, G4ThreeVector( BoxsizeX*(PAD_SegNum/2.+0.5-4.0), 0.0*mm, -OffsetB ) );
  G4UnionSolid *solidBox4 =
    new G4UnionSolid( Cname_+"Box4", solidBox3, solidBoxB,
  		      0, G4ThreeVector( BoxsizeX*(PAD_SegNum/2.+0.5-8.0), 0.0*mm, -OffsetB ) );
  G4UnionSolid *solidBox5 =
    new G4UnionSolid( Cname_+"Box5", solidBox4, solidBoxC,
  		      0, G4ThreeVector( BoxsizeX*(PAD_SegNum/2.+0.5-3.0), 0.0*mm, -OffsetC ) );
  G4UnionSolid *solidBox6 =
    new G4UnionSolid( Cname_+"Box6", solidBox5, solidBoxC,
  		      0, G4ThreeVector( BoxsizeX*(PAD_SegNum/2.+0.5-9.0), 0.0*mm, -OffsetC ) );
  G4UnionSolid *solidBox7 =
    new G4UnionSolid( Cname_+"Box7", solidBox6, solidBoxD,
  		      0, G4ThreeVector( BoxsizeX*(PAD_SegNum/2.+0.5-2.0)-BoxsizeX2/2., 0.0*mm, -OffsetD ) );
  G4UnionSolid *solidBox8 =
    new G4UnionSolid( Cname_+"Box8", solidBox7, solidBoxD,
  		      0, G4ThreeVector( BoxsizeX*(PAD_SegNum/2.+0.5-10.0)+BoxsizeX2/2., 0.0*mm, -OffsetD ) );
  G4UnionSolid *solidBox9 =
    new G4UnionSolid( Cname_+"Box9", solidBox8, solidBoxE,
  		      0, G4ThreeVector( BoxsizeX*(PAD_SegNum/2.+0.5-1.0)-1.5*BoxsizeX2, 0.0*mm, -OffsetE ) );
  G4UnionSolid *solidBox10 =
    new G4UnionSolid( Cname_+"Box10", solidBox9, solidBoxE,
  		      0, G4ThreeVector( BoxsizeX*(PAD_SegNum/2.+0.5-11.0)+1.5*BoxsizeX2, 0.0*mm, -OffsetE ) );

  
  //PAD seg
  G4Box *solidSegA =
    new G4Box( Cname_+"SegmentA", SegsizeX/2., SegsizeZ/2., SegsizeY_A/2. );
  G4Box *solidSegB =
    new G4Box( Cname_+"SegmentB", SegsizeX/2., SegsizeZ/2., SegsizeY_B/2. );
  G4Box *solidSegC =
    new G4Box( Cname_+"SegmentC", SegsizeX/2., SegsizeZ/2., SegsizeY_C/2. );
  G4Box *solidSegD =
    new G4Box( Cname_+"SegmentD", SegsizeX2/2., SegsizeZ/2., SegsizeY_D/2. );
  G4Box *solidSegE =
    new G4Box( Cname_+"SegmentE", SegsizeX2/2., SegsizeZ/2., SegsizeY_E/2. );

  G4UnionSolid *solidSeg1 =
    new G4UnionSolid( Cname_+"Seg1", solidSegA, solidSegA,
  		      0, G4ThreeVector( SegsizeX*(PAD_SegNum/2.+0.5-5.0), 0.0*mm, -OffsetA ) );
  G4UnionSolid *solidSeg2 =
    new G4UnionSolid( Cname_+"Seg2", solidSeg1, solidSegA,
  		      0, G4ThreeVector( SegsizeX*(PAD_SegNum/2.+0.5-7.0), 0.0*mm, -OffsetA ) );
  G4UnionSolid *solidSeg3 =
    new G4UnionSolid( Cname_+"Seg3", solidSeg2, solidSegB,
  		      0, G4ThreeVector( SegsizeX*(PAD_SegNum/2.+0.5-4.0), 0.0*mm, -OffsetB ) );
  G4UnionSolid *solidSeg4 =
    new G4UnionSolid( Cname_+"Seg4", solidSeg3, solidSegB,
  		      0, G4ThreeVector( SegsizeX*(PAD_SegNum/2.+0.5-8.0), 0.0*mm, -OffsetB ) );
  G4UnionSolid *solidSeg5 =
    new G4UnionSolid( Cname_+"Seg5", solidSeg4, solidSegC,
  		      0, G4ThreeVector( SegsizeX*(PAD_SegNum/2.+0.5-3.0), 0.0*mm, -OffsetC ) );
  G4UnionSolid *solidSeg6 =
    new G4UnionSolid( Cname_+"Seg6", solidSeg5, solidSegC,
  		      0, G4ThreeVector( SegsizeX*(PAD_SegNum/2.+0.5-9.0), 0.0*mm, -OffsetC ) );
  G4UnionSolid *solidSeg7 =
    new G4UnionSolid( Cname_+"Seg7", solidSeg6, solidSegD,
  		      0, G4ThreeVector( SegsizeX*(PAD_SegNum/2.+0.5-2.0)-SegsizeX2/2., 0.0*mm, -OffsetD ) );
  G4UnionSolid *solidSeg8 =
    new G4UnionSolid( Cname_+"Seg8", solidSeg7, solidSegD,
  		      0, G4ThreeVector( SegsizeX*(PAD_SegNum/2.+0.5-10.0)+SegsizeX2/2., 0.0*mm, -OffsetD ) );
  G4UnionSolid *solidSeg9 =
    new G4UnionSolid( Cname_+"Seg9", solidSeg8, solidSegE,
  		      0, G4ThreeVector( SegsizeX*(PAD_SegNum/2.+0.5-1.0)-1.5*SegsizeX2, 0.0*mm, -OffsetE ) );
  G4UnionSolid *solidSeg10 =
    new G4UnionSolid( Cname_+"Seg10", solidSeg9, solidSegE,
  		      0, G4ThreeVector( SegsizeX*(PAD_SegNum/2.+0.5-11.0)+1.5*SegsizeX2, 0.0*mm, -OffsetE ) );

  logBox = 
    new G4LogicalVolume( solidBox10, matBox, Cname_+"Box", 0, 0, 0 );
  logScin = 
    new G4LogicalVolume( solidSeg10, matScin, Cname_+"Scin", 0, 0, 0 );

  G4VPhysicalVolume *physBox =
    new G4PVPlacement( G4Transform3D( rotMat,
                                      G4ThreeVector( PADgPos.x()*mm,
                                                     PADgPos.y()*mm,
                                                     PADgPos.z()*mm ) ),
  		       logBox, Cname_+"Box", pMother, false, 0 );

  G4VPhysicalVolume *physScin =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
  		       Cname_+"Seg", logScin, physBox, false, 1 );
}

void s_PAD::SetVisAttributes( const G4VisAttributes *attScin,
			      const G4VisAttributes *attBox )
{ 
  if(attScin) logScin->SetVisAttributes(attScin);
  if(attBox)  logBox->SetVisAttributes(attBox);
}

void s_PAD::SetSensitiveDetector( G4VSensitiveDetector *SD )
{
  logScin->SetSensitiveDetector(SD);
}
