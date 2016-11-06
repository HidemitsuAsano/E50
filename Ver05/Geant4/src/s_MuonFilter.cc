/*
  s_MuonFilter.cc

  Simple Muon Filter -> Just plate
  
  2016/2  K.Shirotori
*/

#include "s_MuonFilter.hh"

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

const double SegSpacing =  3000.0*mm;
const double SegNum     =  1.;

s_MuonFilter::s_MuonFilter( const G4String & Cname, 
			    G4LogicalVolume *pMother,
			    const G4RotationMatrix & rotMat, 
			    const G4ThreeVector & gPosU, 
			    const G4ThreeVector & gPosD, 
			    const G4ThreeVector & OffsetLocal,
			    G4Material *matScin,
			    G4Material *matIron,
			    G4Material *matBox )
  : Cname_(Cname)
{
  G4ThreeVector gPos=(gPosU+gPosD)/2.;
  G4ThreeVector gOffset=rotMat*OffsetLocal;
  G4ThreeVector MFgPos=gPos+gOffset;
  
  //  G4cout << "OffL " << OffsetLocal << "  gOffset " << gOffset
  //	 << " gPos " << gPos << " MFgPos " << MFgPos
  //	 << G4endl;

  double SegsizeX1 =  3000.0*mm;
  double SegsizeY1 =  2500.0*mm;
  double SegsizeZ1 =    30.0*mm;

  double SegsizeX2 =  3000.0*mm;
  double SegsizeY2 =  2500.0*mm;
  double SegsizeZ2 =    30.0*mm;

  double IronsizeX =  3200.0*mm;
  double IronsizeY =  2700.0*mm;
  double IronsizeZ =  1000.0*mm;

  double BoxsizeX = SegsizeX1 + 1.0*mm;
  double BoxsizeY = SegsizeY1 + 1.0*mm;
  double BoxsizeZ = SegsizeZ1 + SegsizeZ2 + IronsizeZ + 2.*30.0*mm + 1.0*mm;

  G4Box *solidBox =
    new G4Box( Cname_+"Box", BoxsizeX/2., BoxsizeY/2., BoxsizeZ/2. );
  G4Box *solidSegU =
    new G4Box( Cname_+"SegmentU", SegsizeX1/2., SegsizeY1/2., SegsizeZ1/2. );
  G4Box *solidSegD =
    new G4Box( Cname_+"SegmentD", SegsizeX2/2., SegsizeY2/2., SegsizeZ2/2. );
  G4Box *solidIron =
    new G4Box( Cname_+"Iron", IronsizeX/2., IronsizeY/2., IronsizeZ/2. );

  logBox = 
    new G4LogicalVolume( solidBox, matBox, Cname_+"Box", 0, 0, 0 );
  logScinU = 
    new G4LogicalVolume( solidSegU, matScin, Cname_+"ScinU", 0, 0, 0 );
  logScinD = 
    new G4LogicalVolume( solidSegD, matScin, Cname_+"ScinD", 0, 0, 0 );
  logIron = 
    new G4LogicalVolume( solidIron, matIron, Cname_+"Iron", 0, 0, 0 );

  G4VPhysicalVolume *physBox =
    new G4PVPlacement( G4Transform3D( rotMat,
                                      G4ThreeVector( MFgPos.x()*mm,
                                                     MFgPos.y()*mm,
                                                     MFgPos.z()*mm ) ),
		       logBox, Cname_+"Box", pMother, false, 0 );

  G4VPhysicalVolume *physIron =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       Cname_+"Iron", logIron, physBox, false, 0 );

  double ofsScinU = SegsizeZ1/2. + IronsizeZ/2. + 30.0*mm;
  double ofsScinD = SegsizeZ2/2. + IronsizeZ/2. + 30.0*mm;

  for( int id=1; id<=SegNum; ++id ){
    G4double ofsScinX=-SegSpacing*(SegNum/2.+0.5-id);
    new G4PVPlacement( 0, G4ThreeVector( ofsScinX, 0.0*mm, (-1.)*ofsScinU ),
		       Cname_+"SegU", logScinU, physBox, false, 1 );
  }
  for( int id=1; id<=SegNum; ++id ){
    G4double ofsScinX=-SegSpacing*(SegNum/2.+0.5-id);
    new G4PVPlacement( 0, G4ThreeVector( ofsScinX, 0.0*mm, ofsScinD ),
		       Cname_+"SegD", logScinD, physBox, false, 2 );
  }

}

void s_MuonFilter::SetVisAttributes( const G4VisAttributes *attScin,
				     const G4VisAttributes *attIron,
				     const G4VisAttributes *attBox )
{
  if(attScin) logScinU->SetVisAttributes(attScin);
  if(attScin) logScinD->SetVisAttributes(attScin);
  if(attIron) logIron->SetVisAttributes(attIron);
  if(attBox) logBox->SetVisAttributes(attBox);
}

void s_MuonFilter::SetSensitiveDetector( G4VSensitiveDetector *SD )
{
  logScinU->SetSensitiveDetector(SD);
  logScinD->SetSensitiveDetector(SD);
}
