/*
  r_SFT_Round.cc

  code for placement of a layer of Segmented SFT
  -Multi Caldding round Fiber + Frame
  -Kuraray's Scintillation fiber SCSF-78M is implemented
  
  10 fibers (r=0.5 mm) x 2 sublayers x 12 layers (x,u,v x,u,v, u,v,x, u,v,x)

  Apr. 27th H.Asano

  
  description of fiber geometry
  
  LzLayer = layer position with respect to the center of physArea.
  physArea is G4VPhysicalVolume , which represents the cage of of fiber tracker.
  The center position of the physArea is the intermediate point between layer 5 and 6.
  The global position of each layer is defined in "param_DetGeom001", however in this GEANT simulation, only the position of layer 5 and 6 is used.
  The position of other layers is used in analysis stage.
  
  global z position of the layer i = LzLayer[i] + (global position of 5 + global position of 6)/2.0

*/

#include "r_SFT_Round.hh"
#include "r_DetectorSize.hh"
#include <cmath>
#include <sstream>
#include <string>


/*-->defined in DetectorConstructionSpec.cc
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
};*/

static const G4double rSFT_TiltAngle  =   30.0*degree;
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

r_SFT_Round::r_SFT_Round( const G4String & Cname,
		  G4VPhysicalVolume *pMother,
		  const G4RotationMatrix &rotMat,
		  const G4ThreeVector *gPos, 
		  const G4ThreeVector &OffsetLocal,
		  const G4int *detid,
      const G4double *LzLayer,
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
  
  G4double Radius     = rSFT_Round_Radius;

  G4double BaseThick = rSFT_BoxThick*12.0 + 20.0*mm;

  G4double radius_core_in = 0*mm;
  G4double radius_core_out = Radius*rSFT_Round_Corefraction;
  G4double radius_innerclad_in = radius_core_out ;
  G4double radius_innerclad_out = radius_innerclad_in+Radius*rSFT_Round_InnerCladding;
  G4double radius_outerclad_in = radius_innerclad_out ;
  G4double radius_outerclad_out =  radius_outerclad_in+Radius*rSFT_Round_OuterCladding;
  G4double FrameSize = rSFT_Round_FramesizeX;
  G4double FrameT = rSFT_Round_FrameThick;

  G4Box *solidBox  = 
    new G4Box( Cname_+"Box", (FrameSize+1.0*mm)/2.+FrameT, (rSFT_Round_Length+1.0*mm)/2.+FrameT, (BaseThick+0.5*mm)/2. );
  G4Box *solidBase  = 
    new G4Box( Cname_+"Base", FrameSize/2.+FrameT, rSFT_Round_Length/2.+FrameT, BaseThick/2. );
  G4Box *solidArea  = 
    new G4Box( Cname_+"Area", FrameSize/2., rSFT_Round_Length/2., BaseThick/2.+0.001*mm );
  
  
  //Note the local coordinate system is different from type A, C (square fiber) because of the symmetry.
  G4Tubs *solidLayer_core = 
    new G4Tubs( Cname_+"Layer", radius_core_in, radius_core_out, rSFT_Round_Length/2.,0.*deg, 360.*deg);
  
  G4Tubs *solidLayer_innerclad = 
    new G4Tubs( Cname_+"Layer_innerclad", radius_innerclad_in, radius_innerclad_out, rSFT_Round_Length/2.,0.*deg, 360.*deg);
  
  G4Tubs *solidLayer_outerclad = 
    new G4Tubs( Cname_+"Layer_outerclad", radius_outerclad_in, radius_outerclad_out, rSFT_Round_Length/2.,0.*deg, 360.*deg);
  
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
//  G4VPhysicalVolume *physFrame =
//    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
//                       Cname_+"Frame", logFrame, physBox, false, 0 );

  G4VPhysicalVolume *physArea =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       Cname_+"Area", logArea, physBox, false, 0 );
 
  // Consitency should be checked.
  // Here, the position located in the GEANT world and parameter file (param/*** ) should be same
  G4RotationMatrix InvMat=rotMat.inverse();
  G4ThreeVector lPos[12];
  for(G4int ilr = 0; ilr < rSFT_nLayer;ilr++){
    lPos[ilr] = InvMat*(gPos[ilr]-gPosCent);
    if( fabs(lPos[ilr].z()-LzLayer[ilr])>MaxDispl){
      G4cout << Cname_ << "Layer" << ilr << ": Geometry Error!:" 
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
  const G4int layer_conf[rSFT_nLayer]={0, 1, 2, 0, 1, 2, 1, 2, 0, 1, 2, 0}; 
  
  G4RotationMatrix *layrot[3];
  //x layer rotation local -> global 
  layrot[0] = new G4RotationMatrix();
  layrot[0]->rotateX(90.*deg);
  
  //u layer rotation local -> global 
  layrot[1] = new G4RotationMatrix();
  layrot[1]->rotateZ( (-1.0)*rSFT_TiltAngle);
  layrot[1]->rotateX( 90.*deg);
  
  //v layer rotation local -> global 
  layrot[2] = new G4RotationMatrix();
  layrot[2]->rotateZ( (1.0)*rSFT_TiltAngle);
  layrot[2]->rotateX( 90.*deg);

  
  //numbering convention
  //layer , segment (= a fiber) and sublayer number starts from 0 , not 1
  for( G4int iseg=0; iseg<rSFT_SegNum; ++iseg ){ 
    G4double ofsScinX=-rSFT_SegSpacing*(rSFT_SegNum/2.+0.5-iseg);
    for( G4int ilr=0; ilr<rSFT_nLayer;ilr++){
      G4int conf = layer_conf[ilr];

      //set name
      std::ostringstream os;
      os << ilr;
      G4String Lname = os.str();;//+"_"+1;

      //first sublayer
      new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX, 0.0*mm, LzLayer[ilr] ),
          Cname_+"Layer"+Lname+"_0_core", logLayer_core, physArea, false, 1000*detid[ilr]+2*iseg );

      new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX, 0.0*mm, LzLayer[ilr] ),
          Cname_+"Layer"+Lname+"_0_innerclad", logLayer_innerclad, physArea, false, 1000*detid[ilr]+2*iseg );

      new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX, 0.0*mm, LzLayer[ilr] ),
          Cname_+"Layer"+Lname+"_0_outerclad", logLayer_outerclad, physArea, false, 1000*detid[ilr]+2*iseg );

      //second sublayer
      new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX+offset_2ndsublayer, 0.0*mm, LzLayer[ilr]+rSFT_Round_Radius*1.732051 ),
          Cname_+"Layer"+Lname+"_1_core", logLayer_core, physArea, false, 1000*detid[ilr]+2*iseg+1 );

      new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX+offset_2ndsublayer, 0.0*mm, LzLayer[ilr]+rSFT_Round_Radius*1.732051 ),
          Cname_+"Layer"+Lname+"_1_innerclad", logLayer_innerclad, physArea, false, 1000*detid[ilr]+2*iseg+1 );

      new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX+offset_2ndsublayer, 0.0*mm, LzLayer[ilr]+rSFT_Round_Radius*1.732051 ),
          Cname_+"Layer"+Lname+"_1_outerclad", logLayer_outerclad, physArea, false, 1000*detid[ilr]+2*iseg+1 );
    }
  }

}

void r_SFT_Round::SetVisAttributes( const G4VisAttributes *attLayer,
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

void r_SFT_Round::SetSensitiveDetector( G4VSensitiveDetector *SD )
{
  logLayer_core->SetSensitiveDetector(SD);
}
