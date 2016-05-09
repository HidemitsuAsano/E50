/*
  r_SFT_C.cc

  Segmented SFT (0.5mm)
  -> square Fiber + Frame

  May 2016 H.Asano

*/

#include "r_SFT_C.hh"
#include "r_DetectorSize.hh"
#include <cmath>
#include <sstream>
#include <string>

const G4double LzLayer[rSFT_nLayer]={
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

#include "globals.hh"
#include "G4Box.hh"
#include "G4Para.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4ios.hh"

r_SFT_C::r_SFT_C( const G4String & Cname,
		  G4VPhysicalVolume *pMother,
		  const G4RotationMatrix &rotMat,
		  const G4ThreeVector *gPos, 
		  const G4ThreeVector &OffsetLocal,
		  G4int *detid, 
		  G4Material *matScin, 
		  G4Material *matPMMA, 
		  G4Material *matFrame,
		  G4Material *matArea,
		  G4Material *matBox )
  : Cname_(Cname)
{
  G4ThreeVector gPosCent=(gPos[5]+gPos[6])/2.;
  G4ThreeVector gOffset=rotMat*OffsetLocal;
  G4ThreeVector DCgPos=gPosCent+gOffset;
  
  G4double LaysizeX     = rSFT_C_LaysizeX;
  G4double LaysizeY     = rSFT_C_LaysizeY;
  G4double Laythickness = rSFT_C_LaysizeZ;

  G4double tilted = rSFT_TiltAngle;
  G4double BaseThick = rSFT_BoxThick*12.0 + 20.0*mm; // 12 number of layer , 20 mm ?

  G4double width = LaysizeX;//
  G4double livewidth = LaysizeX*rSFT_livefraction;
  G4double thickness = Laythickness;
  G4double livethickness = Laythickness*rSFT_livefraction;
  G4double FrameSize = rSFT_FramesizeX;
  G4double FrameT = rSFT_FrameThick;

  G4Box *solidBox  = 
    new G4Box( Cname_+"Box", (FrameSize+1.0*mm)/2.+FrameT, (LaysizeY+1.0*mm)/2.+FrameT, (BaseThick+0.5*mm)/2. );
  G4Box *solidBase  = 
    new G4Box( Cname_+"Base", FrameSize/2.+FrameT, LaysizeY/2.+FrameT, BaseThick/2. );
  G4Box *solidArea  = 
    new G4Box( Cname_+"Area", FrameSize/2., LaysizeY/2., BaseThick/2.+0.001*mm );
  
  G4Box *solidLayerX = 
    new G4Box( Cname_+"LayerX", width/2., LaysizeY/2.,thickness/2. );
  G4Para *solidLayerU = 
    new G4Para( Cname_+"LayerU", width/2., LaysizeY/2., thickness/2.,
		(-1.)*tilted, 0.0*degree, 0.0*degree  );
  G4Para *solidLayerV = 
    new G4Para( Cname_+"LayerV", width/2., LaysizeY/2., thickness/2.,
		tilted, 0.0*degree, 0.0*degree );
  
  G4Box *solidLayerX_core = 
    new G4Box( Cname_+"LayerX_core", livewidth/2., LaysizeY/2.,livethickness/2. );
  G4Para *solidLayerU_core = 
    new G4Para( Cname_+"LayerU_core", livewidth/2., LaysizeY/2., livethickness/2.,
		(-1.)*tilted, 0.0*degree, 0.0*degree  );
  G4Para *solidLayerV_core = 
    new G4Para( Cname_+"LayerV_core", livewidth/2., LaysizeY/2., livethickness/2.,
		tilted, 0.0*degree, 0.0*degree );
  
  G4SubtractionSolid *solidFrame =
    new G4SubtractionSolid( Cname_+"Frame", solidBase, solidArea );

   
  G4SubtractionSolid *solidLayerX_clad =
    new G4SubtractionSolid( Cname_+"LayerX_clad", solidLayerX, solidLayerX_core );

  G4SubtractionSolid *solidLayerU_clad =
    new G4SubtractionSolid( Cname_+"LayerU_clad", solidLayerU, solidLayerU_core );
  
  G4SubtractionSolid *solidLayerV_clad =
    new G4SubtractionSolid( Cname_+"LayerV_clad", solidLayerV, solidLayerV_core );
  
  logBox   = 
    new G4LogicalVolume( solidBox,  matBox, Cname_+"Box",   0, 0, 0 );
  logFrame   = 
    new G4LogicalVolume( solidFrame,  matFrame, Cname_+"Frame",   0, 0, 0 );
  logArea   = 
    new G4LogicalVolume( solidArea,  matArea, Cname_+"Area",   0, 0, 0 );
  
  logLayerX_core = 
    new G4LogicalVolume( solidLayerX_core, matScin, Cname_+"LayerX_core", 0, 0, 0 );
  logLayerU_core = 
    new G4LogicalVolume( solidLayerU_core, matScin, Cname_+"LayerU_core", 0, 0, 0 );
  logLayerV_core = 
    new G4LogicalVolume( solidLayerV_core, matScin, Cname_+"LayerV_core", 0, 0, 0 );
  
  logLayerX_clad = 
    new G4LogicalVolume( solidLayerX_clad, matPMMA, Cname_+"LayerX_clad", 0, 0, 0 );
  logLayerU_clad = 
    new G4LogicalVolume( solidLayerU_clad, matPMMA, Cname_+"LayerU_clad", 0, 0, 0 );
  logLayerV_clad = 
    new G4LogicalVolume( solidLayerV_clad, matPMMA, Cname_+"LayerV_clad", 0, 0, 0 );
  
  G4VPhysicalVolume *physBox =
    new G4PVPlacement( G4Transform3D( rotMat,
				      G4ThreeVector( DCgPos.x()*mm,
						     DCgPos.y()*mm,
						     DCgPos.z()*mm ) ),
		       Cname_+"Box", logBox, pMother, false, 0 );

  G4VPhysicalVolume *physFrame =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       Cname_+"Frame", logFrame, physBox, false, 0 );

  G4VPhysicalVolume *physArea =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       Cname_+"Area", logArea, physBox, false, 0 );
 
  // Consitency should be checked.
  G4RotationMatrix InvMat=rotMat.inverse();
  G4ThreeVector lPos[12];
  for(G4int ilr = 0; ilr < rSFT_nLayer;ilr++){
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
  const G4int layer_conf[rSFT_nLayer]={0, 1, 2, 0, 1,2 ,1,2,0,1,2,0}; 

  for( int id=1; id<=rSFT_C_SegNum; ++id ){ 
    G4double ofsScinX=-rSFT_C_SegSpacing*(rSFT_SegNum/2.+0.5-id);
    
    for(int ilr=0; ilr<rSFT_nLayer; ilr++){
      G4int conf = layer_conf[ilr];

      //set name
      std::ostringstream os;
      os << ilr+1;
      G4String Lname = os.str();;//+"_"+1;

      //
      if(conf==0){
        // first sub layer
        new G4PVPlacement( 0, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer[ilr] ),
            Cname_+"Layer"+Lname+"_1_clad", logLayerX_clad, physArea, false, 1000*detid[ilr]+id );

        new G4PVPlacement( 0, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer[ilr] ),
            Cname_+"Layer"+Lname+"_1_core", logLayerX_core, physArea, false, 1000*detid[ilr]+id );

        // second sub layer
        new G4PVPlacement( 0, G4ThreeVector(  ofsScinX+rSFT_C_LaysizeX/2.0, 0.0*mm, -LzLayer[ilr]+rSFT_C_LaysizeZ ),
            Cname_+"Layer"+Lname+"_2_clad", logLayerX_clad, physArea, false, 1000*detid[ilr]+id+rSFT_C_SegNum );

        new G4PVPlacement( 0, G4ThreeVector(  ofsScinX+rSFT_C_LaysizeX/2.0, 0.0*mm, -LzLayer[ilr]+rSFT_C_LaysizeZ ),
            Cname_+"Layer"+Lname+"_2_core", logLayerX_core, physArea, false, 1000*detid[ilr]+id+rSFT_C_SegNum );

      }else if(conf==1){
        // first sub layer
        new G4PVPlacement( 0, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer[ilr] ),
            Cname_+"Layer"+Lname+"_1_clad", logLayerU_clad, physArea, false, 1000*detid[ilr]+id );

        new G4PVPlacement( 0, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer[ilr] ),
            Cname_+"Layer"+Lname+"_1_core", logLayerU_core, physArea, false, 1000*detid[ilr]+id );

        // second sub layer
        new G4PVPlacement( 0, G4ThreeVector(  ofsScinX+rSFT_C_LaysizeX/2.0, 0.0*mm, -LzLayer[ilr]+rSFT_C_LaysizeZ ),
            Cname_+"Layer"+Lname+"_2_clad", logLayerU_clad, physArea, false, 1000*detid[ilr]+id+rSFT_C_SegNum );

        new G4PVPlacement( 0, G4ThreeVector(  ofsScinX+rSFT_C_LaysizeX/2.0, 0.0*mm, -LzLayer[ilr]+rSFT_C_LaysizeZ ),
            Cname_+"Layer"+Lname+"_2_core", logLayerU_core, physArea, false, 1000*detid[ilr]+id+rSFT_C_SegNum );

      }else if(conf==2){
        // first sub layer
        new G4PVPlacement( 0, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer[ilr] ),
            Cname_+"Layer"+Lname+"_1_clad", logLayerV_clad, physArea, false, 1000*detid[ilr]+id );

        new G4PVPlacement( 0, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer[ilr] ),
            Cname_+"Layer"+Lname+"_1_core", logLayerV_core, physArea, false, 1000*detid[ilr]+id );

        // second sub layer
        new G4PVPlacement( 0, G4ThreeVector(  ofsScinX+rSFT_C_LaysizeX/2.0, 0.0*mm, -LzLayer[ilr]+rSFT_C_LaysizeZ ),
            Cname_+"Layer"+Lname+"_2_clad", logLayerV_clad, physArea, false, 1000*detid[ilr]+id+rSFT_C_SegNum );

        new G4PVPlacement( 0, G4ThreeVector(  ofsScinX+rSFT_C_LaysizeX/2.0, 0.0*mm, -LzLayer[ilr]+rSFT_C_LaysizeZ ),
            Cname_+"Layer"+Lname+"_2_core", logLayerV_core, physArea, false, 1000*detid[ilr]+id+rSFT_C_SegNum );
      }
    }//ilr
  }//iseg
}

void r_SFT_C::SetVisAttributes( const G4VisAttributes *attLayer,
        const G4VisAttributes *attClad,
				const G4VisAttributes *attFrame,
				const G4VisAttributes *attArea,
				const G4VisAttributes *attBox )
{
  if(attLayer) logLayerX_core->SetVisAttributes(attLayer); 
  if(attLayer) logLayerU_core->SetVisAttributes(attLayer); 
  if(attLayer) logLayerV_core->SetVisAttributes(attLayer); 
  if(attLayer) logLayerX_clad->SetVisAttributes(attClad); 
  if(attLayer) logLayerU_clad->SetVisAttributes(attClad); 
  if(attLayer) logLayerV_clad->SetVisAttributes(attClad); 
  if(attFrame) logFrame->SetVisAttributes(attFrame); 
  if(attArea)  logArea->SetVisAttributes(attArea); 
  if(attBox)   logBox->SetVisAttributes(attBox); 
}

void r_SFT_C::SetSensitiveDetector( G4VSensitiveDetector *SD )
{
  logLayerX_core->SetSensitiveDetector(SD);
  logLayerU_core->SetSensitiveDetector(SD);
  logLayerV_core->SetSensitiveDetector(SD);
}
