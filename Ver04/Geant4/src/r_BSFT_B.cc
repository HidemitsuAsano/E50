/*
  r_BSFT_B.cc

  code for placement of Segmented SFT
  -Multi Caldding round Fiber + Frame
  -Kuraray's Scintillation fiber SCSF-78M is implemented
  
  10 fibers (r=0.5 mm) x 2 sublayers x 12 layers (x,u,v x,u,v, u,v,x, u,v,x)

  Apr. 27th H.Asano
*/

#include "r_BSFT_B.hh"
#include "r_DetectorSize.hh"
#include <cmath>
#include <sstream>
#include <string>

const G4double LzLayer[rBSFT_nLayer]={
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
//#include "G4Para.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4ios.hh"

r_BSFT_B::r_BSFT_B( const G4String & Cname,
		  G4VPhysicalVolume *pMother,
		  const G4RotationMatrix &rotMat,
		  const G4ThreeVector *gPos, 
		  const G4ThreeVector &OffsetLocal,
		  G4int *detid,
		  G4Material *matScin, 
		  G4Material *matPMMA, 
		  G4Material *matFP, 
		  G4Material *matFrame,
		  G4Material *matArea,
		  G4Material *matBox )
  : Cname_(Cname)
{
  G4ThreeVector gPosCent=(gPos[5]+gPos[6])/2.;
  G4ThreeVector gOffset=rotMat*OffsetLocal;
  G4ThreeVector DCgPos=gPosCent+gOffset;
  
  G4double Radius     = rBSFT_B_Radius;

  G4double BaseThick = rBSFT_BoxThick*12.0 + 20.0*mm;

  G4double radius_core_in = 0*mm;
  G4double radius_core_out = Radius*rBSFT_B_Corefraction;
  G4double radius_innerclad_in = radius_core_out ;
  G4double radius_innerclad_out = radius_innerclad_in+Radius*rBSFT_B_InnerCladding;
  G4double radius_outerclad_in = radius_innerclad_out ;
  G4double radius_outerclad_out =  radius_outerclad_in+Radius*rBSFT_B_OuterCladding;
  G4double FrameSize = rBSFT_B_FramesizeX;
  G4double FrameT = rBSFT_B_FrameThick;

  G4Box *solidBox  = 
    new G4Box( Cname_+"Box", (FrameSize+1.0*mm)/2.+FrameT, (rBSFT_B_Length+1.0*mm)/2.+FrameT, (BaseThick+0.5*mm)/2. );
  G4Box *solidBase  = 
    new G4Box( Cname_+"Base", FrameSize/2.+FrameT, rBSFT_B_Length/2.+FrameT, BaseThick/2. );
  G4Box *solidArea  = 
    new G4Box( Cname_+"Area", FrameSize/2., rBSFT_B_Length/2., BaseThick/2.+0.001*mm );
  
  
  //Note the local coordinate system is different from type A, C (square fiber) because of the symmetry.
  G4Tubs *solidLayer_core = 
    new G4Tubs( Cname_+"Layer", radius_core_in, radius_core_out, rBSFT_B_Length/2.,0.*deg, 360.*deg);
  
  G4Tubs *solidLayer_innerclad = 
    new G4Tubs( Cname_+"Layer_innerclad", radius_innerclad_in, radius_innerclad_out, rBSFT_B_Length/2.,0.*deg, 360.*deg);
  
  G4Tubs *solidLayer_outerclad = 
    new G4Tubs( Cname_+"Layer_outerclad", radius_outerclad_in, radius_outerclad_out, rBSFT_B_Length/2.,0.*deg, 360.*deg);
  
  G4SubtractionSolid *solidFrame =
    new G4SubtractionSolid( Cname_+"Frame", solidBase, solidArea );

   
  logBox   = 
    new G4LogicalVolume( solidBox,  matBox, Cname_+"Box",   0, 0, 0 );
  logFrame   = 
    new G4LogicalVolume( solidFrame,  matFrame, Cname_+"Frame",   0, 0, 0 );
  logArea   = 
    new G4LogicalVolume( solidArea,  matArea, Cname_+"Area",   0, 0, 0 );
  
  logLayer_core = 
    new G4LogicalVolume( solidLayer_core, matScin, Cname_+"Layer_core", 0, 0, 0 );
  
  logLayer_innerclad = 
    new G4LogicalVolume( solidLayer_innerclad, matPMMA, Cname_+"Layer_innerclad", 0, 0, 0 );
  
  logLayer_outerclad = 
    new G4LogicalVolume( solidLayer_outerclad, matFP, Cname_+"Layer_outerclad", 0, 0, 0 );



  G4VPhysicalVolume *physBox =
    new G4PVPlacement( G4Transform3D( rotMat,
				      G4ThreeVector( DCgPos.x()*mm,
						     DCgPos.y()*mm,
						     DCgPos.z()*mm ) ),
		       Cname_+"Box", logBox, pMother, false, 0 );//name, logical volume, mother volume, no booliean operations its copy number
  G4VPhysicalVolume *physFrame =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       Cname_+"Frame", logFrame, physBox, false, 0 );

  G4VPhysicalVolume *physArea =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       Cname_+"Area", logArea, physBox, false, 0 );
 
  // Consitency should be checked.
  G4RotationMatrix InvMat=rotMat.inverse();
  G4ThreeVector lPos[12];
  for(G4int ilr = 0; ilr < rBSFT_nLayer;ilr++){
    lPos[ilr] = InvMat*(gPos[ilr]-gPosCent);
  
    if( fabs(lPos[ilr].z()-LzLayer[ilr])>MaxDispl){
      G4cout << Cname_ << "Layer" << ilr+1 << ": Geometry Error!:" 
        << "(0,0," << LzLayer[ilr] << ") --> " << lPos[ilr] << G4endl;
      G4Exception("Geometry Error");
    }
  }
 

  const G4double offset_2ndsublayer = 0.5*mm;
  //x u v x u v u v x u v x 
  //Here
  //x: 0 
  //u: 1
  //v: 2
  const G4int layer_conf[rBSFT_nLayer]={0, 1, 2, 0, 1, 2, 1, 2, 0, 1, 2, 0}; 
  
  G4RotationMatrix *layrot[3];
  //x layer local -> global rotation
  layrot[0] = new G4RotationMatrix();
  layrot[0]->rotateX(90.*deg);
  
  //u layer local -> global rotation
  layrot[1] = new G4RotationMatrix();
  layrot[1]->rotateZ( (-1.0)*rBSFT_TiltAngle);
  layrot[1]->rotateX( 90.*deg);
  
  //v layer local -> global rotation
  layrot[2] = new G4RotationMatrix();
  layrot[2]->rotateZ( (1.0)*rBSFT_TiltAngle);
  layrot[2]->rotateX( 90.*deg);

    
  for( G4int iseg=1; iseg<=rBSFT_SegNum; ++iseg ){ 
    G4double ofsScinX=-rBSFT_SegSpacing*(rBSFT_SegNum/2.+0.5-iseg);
    for( G4int ilr=0; ilr<rBSFT_nLayer;ilr++){
      G4int conf = layer_conf[ilr];

      //set name
      std::ostringstream os;
      os << ilr+1;
      G4String Lname = os.str();;//+"_"+1;

      //first sublayer
      new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer[ilr] ),
          Cname_+"Layer"+Lname+"_1_core", logLayer_core, physArea, false, 1000*detid[ilr]+2*iseg-1 );

      new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer[ilr] ),
          Cname_+"Layer"+Lname+"_1_innerclad", logLayer_innerclad, physArea, false, 1000*detid[ilr]+2*iseg-1 );

      new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer[ilr] ),
          Cname_+"Layer"+Lname+"_1_outerclad", logLayer_outerclad, physArea, false, 1000*detid[ilr]+2*iseg-1 );

      //second sublayer
      new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX+offset_2ndsublayer, 0.0*mm, -LzLayer[ilr]+rBSFT_B_Radius*1.732051 ),
          Cname_+"Layer"+Lname+"_2_core", logLayer_core, physArea, false, 1000*detid[ilr]+2*iseg );

      new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX+offset_2ndsublayer, 0.0*mm, -LzLayer[ilr]+rBSFT_B_Radius*1.732051 ),
          Cname_+"Layer"+Lname+"_2_innerclad", logLayer_innerclad, physArea, false, 1000*detid[ilr]+2*iseg );

      new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX+offset_2ndsublayer, 0.0*mm, -LzLayer[ilr]+rBSFT_B_Radius*1.732051 ),
          Cname_+"Layer"+Lname+"_2_outerclad", logLayer_outerclad, physArea, false, 1000*detid[ilr]+2*iseg );
    }
    
  }
}

void r_BSFT_B::SetVisAttributes( const G4VisAttributes *attLayer,
        const G4VisAttributes *attInnerClad,
        const G4VisAttributes *attOuterClad,
				const G4VisAttributes *attFrame,
				const G4VisAttributes *attArea,
				const G4VisAttributes *attBox )
{
  if(attLayer) logLayer_core->SetVisAttributes(attLayer); 
  if(attLayer) logLayer_innerclad->SetVisAttributes(attInnerClad); 
  if(attLayer) logLayer_outerclad->SetVisAttributes(attOuterClad); 
  if(attFrame) logFrame->SetVisAttributes(attFrame); 
  if(attArea)  logArea->SetVisAttributes(attArea); 
  if(attBox)   logBox->SetVisAttributes(attBox); 
}

void r_BSFT_B::SetSensitiveDetector( G4VSensitiveDetector *SD )
{
  logLayer_core->SetSensitiveDetector(SD);
}
