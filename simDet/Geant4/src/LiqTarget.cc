/*
  LigTarget.cc

  2012/5  K.Shirotori
*/

#include "LiqTarget.hh"

#include "globals.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4ios.hh"

const double TargetLength = 120.0*mm;
const double TargetDiameter = 67.2*mm;

LiqTarget::LiqTarget( const G4String & Cname, G4VPhysicalVolume *pMother,
		      const G4RotationMatrix & rotMat,
		      const G4ThreeVector & gPos, 
		      const G4ThreeVector & OffsetLocal,
		      G4Material *matTarget, G4Material *matTargCell,
		      G4Material *matVessel, G4Material *matWindow,
		      G4Material *matVacuum, G4Material *matBox )
  : Cname_(Cname), length_(TargetLength), diameter_(TargetDiameter)
{
  G4ThreeVector gOffset=rotMat*OffsetLocal;
  G4ThreeVector TargPos=gPos+gOffset;

  G4RotationMatrix RMCell; RMCell.rotateX(90.*degree);

  G4Box *solidBox =
    new G4Box( Cname_+"Box", 150.0*mm, 170.0*mm, 450.0*mm ); 

  //  G4VSolid *solidVessel = solidVessel3;

  G4Tubs *solidVacuum0 =
    new G4Tubs( Cname_+"Vacuum0", 0.*mm, 135.*mm, 375.*mm,
		0.0*degree, 360.0*degree );
  G4Tubs *solidVacuum1 =
    new G4Tubs( Cname_+"Vacuum1", 0.*mm, 132.*mm, 375.*mm,
		0.0*degree, 360.0*degree );
  G4Tubs *solidVacuum2 =
    new G4Tubs( Cname_+"Vacuum2", 0.*mm,  50.*mm, 67.5*mm,
		0.0*degree, 360.0*degree );
  G4Tubs *solidVacuum3 =
    new G4Tubs( Cname_+"Vacuum3", 0.*mm,  45.*mm, 81.0*mm,
		0.0*degree, 360.0*degree );
  G4UnionSolid *solidVacuum4 =
    new G4UnionSolid( Cname_+"Vacuum4", solidVacuum1, solidVacuum2,
		      G4Transform3D( RMCell,
				     G4ThreeVector( 0., -67.5*mm, 0. ) ) );
  G4IntersectionSolid *solidVacuum6 =
    new G4IntersectionSolid( Cname_+"Vacuum6", solidVacuum4, solidVacuum0,
			     0, G4ThreeVector( 0., 0., 0. ) );
  G4UnionSolid *solidVacuum5 =
    new G4UnionSolid( Cname_+"Vacuum5", solidVacuum6, solidVacuum3,
		      G4Transform3D( RMCell,
				     G4ThreeVector( 0., 81.*mm, 0. ) ) );

  G4VSolid *solidVacuum = solidVacuum5;

  G4Tubs *solidVessel1 =
    new G4Tubs( Cname_+"Vessel1", 132.*mm, 135.*mm, 375.*mm,
		0.0*degree, 360.0*degree );
  G4Tubs *solidVessel2 =
    new G4Tubs( Cname_+"Vessel2", 0.*mm, 50.*mm, 136.0*mm,
		0.0*degree, 360.0*degree );

  G4Tubs *solidVessel4 =
    new G4Tubs( Cname_+"Vessel4", 45.*mm, 50.*mm, 19.*mm,
		0.0*degree, 360.0*degree );
  G4Tubs *solidVessel5 =
    new G4Tubs( Cname_+"Vessel5", 50.*mm, 75.*mm, 10.*mm, 
		0.0*degree, 360.0*degree );

  G4UnionSolid *solidVesselU1 =
    new G4UnionSolid( Cname_+"VesselU1", solidVessel1, solidVessel4,
		      G4Transform3D( RMCell,
				     G4ThreeVector( 0., 143.*mm, 0. ) ) );
  G4UnionSolid *solidVesselU2 =
    new G4UnionSolid( Cname_+"VesselU2", solidVesselU1, solidVessel5,
		      G4Transform3D( RMCell,
				     G4ThreeVector( 0., 152.*mm, 0. ) ) );
  //  G4VSolid *solidVessel = solidVesselU2;

  G4SubtractionSolid *solidVessel = 
   new G4SubtractionSolid( Cname_+"Vessel", solidVesselU2, solidVacuum5,
		    0, G4ThreeVector( 0., 0., 0. ) );

  G4Tubs *solidCell1 =
    new G4Tubs( Cname_+"Cell1", diameter_/2., diameter_/2.+0.3*mm, 
		length_/2., 0.0*degree, 360.0*degree );
  G4Tubs *solidCell2 =
    new G4Tubs( Cname_+"Cell2",  0.0*mm, diameter_/2.+0.3*mm, 0.125*mm, 
		0.0*degree, 360.0*degree );
  G4UnionSolid *solidCell3 =
    new G4UnionSolid( Cname_+"Cell3", solidCell1, solidCell2,
		      0, G4ThreeVector( 0.0, 0.0, length_/2.+0.125*mm ) );
  G4UnionSolid *solidCell =
    new G4UnionSolid( Cname_+"Cell", solidCell3, solidCell2,
		      0, G4ThreeVector( 0.0, 0.0, -length_/2.-0.125*mm ) );

  G4Tubs *solidTarg =
    new G4Tubs( Cname_+"Targ", 0.0, diameter_/2., length_/2.,
		0.0*degree, 360.0*degree );

  G4Tubs *solidWinIn =
    new G4Tubs( Cname_+"WinIn", 0.0*mm, 45.0*mm, 0.1*mm,
                0.0*degree, 360.0*degree );
  G4Tubs *solidWinOut1 =
    new G4Tubs( Cname_+"WinOut1", 134.8*mm, 135.0*mm, 100.0*mm, 
                0.0*degree, 360.0*degree );
  G4IntersectionSolid *solidWinOut =
    new G4IntersectionSolid( Cname_+"WinOut", solidWinOut1, solidVacuum2,
			     G4Transform3D( RMCell,
					    G4ThreeVector( 0., -67.5*mm, 0. ) )
			     );



  logBox =
    new G4LogicalVolume( solidBox, matBox, Cname_+"Box", 0, 0, 0 );

  logVessel =
    new G4LogicalVolume( solidVessel, matVessel, Cname_+"Vessel", 0, 0, 0 );
  logVacuum =
    new G4LogicalVolume( solidVacuum, matVacuum, Cname_+"Vacuum", 0, 0, 0 );

  logCell = 
    new G4LogicalVolume( solidCell, matTargCell, Cname_+"Cell", 0, 0, 0 );
  logTarg =
    new G4LogicalVolume( solidTarg, matTarget, Cname_+"Targ", 0, 0, 0 );

  logWindowIn =
    new G4LogicalVolume( solidWinIn, matWindow, 
			 Cname_+"WindowIn", 0, 0, 0 );
  logWindowOut =
    new G4LogicalVolume( solidWinOut, matWindow, 
			 Cname_+"WindowOut", 0, 0, 0 );

  G4VPhysicalVolume *physBox =
    new G4PVPlacement( G4Transform3D( rotMat,
				      G4ThreeVector( TargPos.x()*mm,
						     TargPos.y()*mm,
						     TargPos.z()*mm ) ),
		       Cname_+"Box", logBox, pMother, false, 0 );

  G4VPhysicalVolume *physVessel =
    new G4PVPlacement( 0, G4ThreeVector( 0., 0., 0. ),
		       Cname_+"Vessel", logVessel, physBox, false, 0 );

  G4VPhysicalVolume *physVacuum =
    new G4PVPlacement( 0, G4ThreeVector( 0., 0., 0. ),
		       Cname_+"Vacuum", logVacuum, physBox, false, 0 );
  G4VPhysicalVolume *physCell =
    new G4PVPlacement( G4Transform3D( RMCell, G4ThreeVector( 0., 0., 0. ) ),
		       Cname_+"Cell", logCell, physVacuum, false, 0 );
  G4VPhysicalVolume *physTarg =
    new G4PVPlacement( 0, G4ThreeVector( 0., 0., 0. ),
		       Cname_+"Targ", logTarg, physCell, false, 0 );
  G4VPhysicalVolume *physWindowIn =
    new G4PVPlacement( G4Transform3D( RMCell, 
				      G4ThreeVector( 0., 161.9*mm, 0. ) ),
		       Cname_+"WindowIn", logWindowIn, physVacuum, 
		       false, 0 );
  G4VPhysicalVolume *physWindowOut =
    new G4PVPlacement( 0, G4ThreeVector( 0., 0., 0. ),
		       Cname_+"WindowOut", logWindowOut, physVacuum, 
		       false, 0 );
}

void LiqTarget::SetVisAttributes( const G4VisAttributes *attTarg,
				  const G4VisAttributes *attCell,
				  const G4VisAttributes *attVessel,
				  const G4VisAttributes *attVacuum,
				  const G4VisAttributes *attWindow,
				  const G4VisAttributes *attBox )
{
  if(attTarg) logTarg->SetVisAttributes(attTarg);
  if(attCell) logCell->SetVisAttributes(attCell);
  if(attVessel) logVessel->SetVisAttributes(attVessel);
  if(attVacuum) logVacuum->SetVisAttributes(attVacuum);
  if(attWindow){
    logWindowIn->SetVisAttributes(attWindow);
    logWindowOut->SetVisAttributes(attWindow);
  }
  if(attBox) logBox->SetVisAttributes(attBox);
}


void LiqTarget::SetSensitiveDetector( G4VSensitiveDetector *SD )
{

}

