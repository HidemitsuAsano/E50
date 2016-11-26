/*
  s_SFT.cc

  Simple SFT
  -> Just plate + Frame

  2016/4  K.Shirotori
  2016/11  Updated by H.Asano
*/


#include "s_SFT.hh"
#include "s_DetectorSize.hh"
#include <cmath>


static const G4int sSFT_nLayer = 12;

const G4double LzLayer[sSFT_nLayer]={
-165.0*mm,
-135.0*mm,
-105.0*mm,
 -75.0*mm,
 -45.0*mm,
 -15.0*mm,
  15.0*mm,
  45.0*mm,
  75.0*mm,
 105.0*mm,
 135.0*mm,
 165.0*mm
};


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

s_SFT::s_SFT( const G4String & Cname,
	      G4VPhysicalVolume *pMother,
	      const G4RotationMatrix &rotMat,
	      const G4ThreeVector *gPos, 
	      const G4ThreeVector &OffsetLocal,
	      const G4int *detid, 
	      G4Material *matScin, 
	      G4Material *matFrame,
	      G4Material *matArea,
	      G4Material *matBox )
  : Cname_(Cname)
{
  G4ThreeVector gPosCent=(gPos[5]+gPos[6])/2.;
  G4ThreeVector gOffset=rotMat*OffsetLocal;
  G4ThreeVector DCgPos=gPosCent+gOffset;
  
  G4double LaysizeX     = sSFT_LaysizeX;
  G4double LaysizeY     = sSFT_LaysizeY;
  G4double Laythickness = sSFT_LaysizeZ;

  G4double tilted = sSFT_TiltAngle;
  G4double BaseThick = sSFT_BoxThick*12.0 + 20.0*mm;

  G4double width = LaysizeX;//+LaysizeY*tan(tilted);
  G4double FrameSize = sSFT_FramesizeX;
  G4double FrameT = sSFT_FrameThick;

  G4Box *solidBox  = 
    new G4Box( Cname_+"Box", (FrameSize+1.0*mm)/2.+FrameT, (LaysizeY+1.0*mm)/2.+FrameT, (BaseThick+0.5*mm)/2. );
  G4Box *solidBase  = 
    new G4Box( Cname_+"Base", FrameSize/2.+FrameT, LaysizeY/2.+FrameT, BaseThick/2. );
  G4Box *solidArea  = 
    new G4Box( Cname_+"Area", FrameSize/2., LaysizeY/2., BaseThick/2.+0.001*mm );
  
  G4Box *solidLayerX = 
    new G4Box( Cname_+"LayerX", LaysizeX/2., LaysizeY/2.,Laythickness/2. );
  G4Para *solidLayerU = 
    new G4Para( Cname_+"LayerU", width/2., LaysizeY/2., Laythickness/2.,
		(-1.)*tilted, 0.0*degree, 0.0*degree  );
  G4Para *solidLayerV = 
    new G4Para( Cname_+"LayerV", width/2., LaysizeY/2., Laythickness/2.,
		tilted, 0.0*degree, 0.0*degree );
  
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
  logLayerU = 
    new G4LogicalVolume( solidLayerU, matScin, Cname_+"LayerU", 0, 0, 0 );
  logLayerV = 
    new G4LogicalVolume( solidLayerV, matScin, Cname_+"LayerV", 0, 0, 0 );
  
  G4VPhysicalVolume *physBox =
    new G4PVPlacement( G4Transform3D( rotMat,
				      G4ThreeVector( DCgPos.x()*mm,
						     DCgPos.y()*mm,
						     DCgPos.z()*mm ) ),
		       Cname_+"Box", logBox, pMother, false, 0 );

//  G4VPhysicalVolume *physFrame =
//    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
//                       Cname_+"Frame", logFrame, physBox, false, 0 );

  G4VPhysicalVolume *physArea =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       Cname_+"Area", logArea, physBox, false, 0 );
 
  // Consitency should be checked.
  G4RotationMatrix InvMat=rotMat.inverse();
  G4ThreeVector lPos[12];
  for(G4int ilr = 0; ilr < sSFT_nLayer;ilr++){
    lPos[ilr] = InvMat*(gPos[ilr]-gPosCent);
    if( fabs(lPos[ilr].z()-LzLayer[ilr])>MaxDispl){
      G4cout << Cname_ << "Layer" << ilr+1 << ": Geometry Error!:" 
        << "(0,0," << LzLayer[ilr] << ") --> " << lPos[ilr] << G4endl;
      G4Exception("Geometry Error");
    }
  }
  
  //x u v x u v u v x u v x 
  //Here
  //x: 0 
  //u: 1
  //v: 2
  const G4int layer_conf[sSFT_nLayer]={0, 1, 2, 0, 1, 2, 1, 2, 0, 1, 2, 0}; 
  G4LogicalVolume *logLayer[3]={logLayerX,logLayerU,logLayerV};

  G4VPhysicalVolume *physLayer[sSFT_nLayer];
  for(int ilr=0;ilr<sSFT_nLayer;ilr++){
    physLayer[ilr] = new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, -LzLayer[ilr]*mm ),
        Cname_+"Layer1", logLayer[layer_conf[ilr]], physArea, false, detid[ilr] );
  }

}

void s_SFT::SetVisAttributes( const G4VisAttributes *attLayer,
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

void s_SFT::SetSensitiveDetector( G4VSensitiveDetector *SD )
{
  logLayerX->SetSensitiveDetector(SD);
  logLayerU->SetSensitiveDetector(SD);
  logLayerV->SetSensitiveDetector(SD);
}
