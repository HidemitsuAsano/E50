/*
  r_T0Wall.cc

  Real T0

  2015/12  K.Shirotori
*/

#include "r_T0Wall.hh"

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

const double SegSpacing =  3.2*mm;
const double SegNum     = 35.;

r_T0Wall::r_T0Wall( const G4String & Cname, 
		    G4LogicalVolume *pMother,
		    const G4RotationMatrix & rotMat, 
		    const G4ThreeVector & gPos, 
		    const G4ThreeVector & OffsetLocal,
		    G4Material *matScin, 
		    G4Material *matMPPC,
		    G4Material *matPCB, 
		    G4Material *matBox )
  : Cname_(Cname)
{
  G4ThreeVector gOffset=rotMat*OffsetLocal;
  G4ThreeVector T0gPos=gPos+gOffset;

  //  G4cout << "OffL " << OffsetLocal << "  gOffset " << gOffset
  //	 << " gPos " << gPos << " T0gPos " << T0gPos
  //	 << G4endl;

  double BoxsizeX =  150.0*mm;
  double BoxsizeY =  150.0*mm;
  double BoxsizeZ =   10.0*mm;

  double SegsizeX =    3.0*mm;
  double SegsizeY =  100.0*mm;
  double SegsizeZ =    3.0*mm;

  double MPPCsizeX =   3.0*mm;
  double MPPCsizeY =   5.0*mm;
  double MPPCsizeZ =   3.0*mm;

  double PCBsizeX =  110.0*mm;
  double PCBsizeY =   50.0*mm;
  double PCBsizeZ =    5.0*mm;

  G4Box *solidBox =
    new G4Box( Cname_+"Box", BoxsizeX/2., BoxsizeY/2., BoxsizeZ/2. );
  G4Box *solidSeg =
    new G4Box( Cname_+"Segment", SegsizeX/2., SegsizeY/2., SegsizeZ/2. );
  G4Box *solidMPPC =
    new G4Box( Cname_+"MPPC", MPPCsizeX/2., MPPCsizeY/2., MPPCsizeZ/2. );
  G4Box *solidPCB =
    new G4Box( Cname_+"PCB", PCBsizeX/2., PCBsizeY/2., PCBsizeZ/2. );

  logBox = 
    new G4LogicalVolume( solidBox, matBox, Cname_+"Box", 0, 0, 0 );
  logScin = 
    new G4LogicalVolume( solidSeg, matScin, Cname_+"Scin", 0, 0, 0 );
  logMPPC =
    new G4LogicalVolume( solidMPPC, matMPPC, Cname_+"MPPC", 0, 0, 0 );
  logPCB =
    new G4LogicalVolume( solidPCB, matPCB, Cname_+"PCB", 0, 0, 0 );

  
  G4VPhysicalVolume *physBox =
    new G4PVPlacement( G4Transform3D( rotMat,
                                      G4ThreeVector( T0gPos.x()*mm,
                                                     T0gPos.y()*mm,
                                                     T0gPos.z()*mm ) ),
		       logBox, Cname_+"Box", pMother, false, 0 );
  
  for( int id=1; id<=SegNum; ++id ){
    G4double ofsScinX=-SegSpacing*(SegNum/2.+0.5-id);
    new G4PVPlacement( 0, G4ThreeVector( ofsScinX, 0.0*mm, 0.0*mm ),
		       Cname_+"Seg", logScin, physBox, false, id );

    new G4PVPlacement( 0, G4ThreeVector( ofsScinX, (SegsizeY/2.+MPPCsizeY/2.), 0.0*mm ),
		       Cname_+"MPPCu", logMPPC, physBox, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector( ofsScinX, (-1.)*(SegsizeY/2.+MPPCsizeY/2.), 0.0*mm ),
		       Cname_+"MPPCd", logMPPC, physBox, false, 0 );
  }
  {
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, (SegsizeY/2.+MPPCsizeY/2.+PCBsizeY/2.), 0.0*mm ),
		       Cname_+"PCBu", logPCB, physBox, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, (-1.)*(SegsizeY/2.+MPPCsizeY/2.+PCBsizeY/2.), 0.0*mm ),
		       Cname_+"PCBd", logPCB, physBox, false, 0 );
  }
}

void r_T0Wall::SetVisAttributes( const G4VisAttributes *attScin,
				const G4VisAttributes *attMPPC,
				const G4VisAttributes *attPCB,
				const G4VisAttributes *attBox )
{
  if(attScin) logScin->SetVisAttributes(attScin);
  if(attMPPC) logMPPC->SetVisAttributes(attMPPC);
  if(attPCB) logPCB->SetVisAttributes(attPCB);
  if(attBox) logBox->SetVisAttributes(attBox);
}

void r_T0Wall::SetSensitiveDetector( G4VSensitiveDetector *SD )
{
  logScin->SetSensitiveDetector(SD);
}

